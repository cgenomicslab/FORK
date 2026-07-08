"""
Visualization helper script for the FORK.

All functions return a BytesIO image buffer that can be passed directly
to st.image() in Streamlit, or saved to disk with open(..., "wb").write(buf.getvalue()).

Dependencies:
    matplotlib, seaborn, biopython, numpy

Usage example (in the GUI):
    import viz_utils as viz
    buf = viz.draw_domain_architecture(domain_records)
    st.image(buf, use_container_width=True)
"""

import hashlib
import io
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")  # Non-interactive backend — must come before pyplot import.
# "Agg" renders to an in-memory buffer instead of a screen
# window. This is required in server/headless environments
# like the lab server running Streamlit.
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _domain_color(name: str) -> tuple:
    """
    Derive a consistent RGBA color from a domain name string.

    How it works:
    - Take the first 6 hex characters of the MD5 hash of the name.
    - Interpret them as an integer, map to a hue in [0, 1].
    - Build an HSV color with fixed saturation (0.65) and value (0.88)
      so all domain colors are equally vivid — no domain gets an
      accidentally washed-out or near-black color.
    - Convert to RGBA with alpha=0.88.

    The result is deterministic: the same domain name always produces
    the same color, across different function calls and different sessions.
    """
    h = int(hashlib.md5(name.encode()).hexdigest()[:6], 16)
    hue = (h % 360) / 360.0
    # HSV → RGB, then add alpha
    import colorsys

    r, g, b = colorsys.hsv_to_rgb(hue, 0.65, 0.88)
    return (r, g, b, 0.88)


def _fig_to_buffer(fig) -> io.BytesIO:
    """
    Render a matplotlib figure to a PNG BytesIO buffer and close the figure.

    Always call this instead of plt.savefig() to avoid accumulating open
    figure objects, which would slowly eat server memory over a long session.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    plt.close(fig)
    return buf


# ---------------------------------------------------------------------------
# 1. Domain architecture diagram
# ---------------------------------------------------------------------------


def draw_domain_architecture(
    domain_records: list,
    protein_lengths: dict = None,
    title: str = None,
    max_proteins: int = 40,
    min_label_width_fraction: float = 0.04,
) -> io.BytesIO:
    """
    Draw a domain architecture diagram: one horizontal bar per protein,
    colored domain blocks overlaid at their alignment positions.


    Example
    -------
    domains = uni.fetch_domains_by_accession("2026_01", acc_list, db_config=config)
    buf = viz.draw_domain_architecture(domains, title="Homeodomain proteins")
    st.image(buf, use_container_width=True)
    """
    if not domain_records:
        # Return a minimal figure with a message rather than crashing.
        fig, ax = plt.subplots(figsize=(8, 2))
        ax.text(
            0.5,
            0.5,
            "No domain records to display.",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
        )
        ax.axis("off")
        return _fig_to_buffer(fig)

    # ---------------- Group domain records by accession ----------------
    # We want: protein_order (stable, insertion order), and for each protein
    # a list of (hmm_name, ali_from, ali_to) tuples sorted by ali_from.
    domains_by_acc = defaultdict(list)
    meta_by_acc = {}  # accession → display label

    for rec in domain_records:
        acc = rec["accession"]
        domains_by_acc[acc].append(rec)
        if acc not in meta_by_acc:
            name = rec.get("protein_name") or ""
            org = rec.get("organism") or ""
            meta_by_acc[acc] = (
                f"{acc}"
                + (f" | {name}" if name else "")
                + (f" [{org[:25]}]" if org else "")
            )

    # Sort domains within each protein by start position
    for acc in domains_by_acc:
        domains_by_acc[acc].sort(key=lambda r: r.get("ali_from", 0))

    protein_list = list(domains_by_acc.keys())
    truncated = len(protein_list) > max_proteins
    protein_list = protein_list[:max_proteins]

    # ---------------- Resolve protein lengths ----------------
    # Priority: (1) caller-supplied dict, (2) "length" field in record,
    # (3) max(env_to) fallback.
    lengths = {}
    for acc in protein_list:
        if protein_lengths and acc in protein_lengths:
            lengths[acc] = protein_lengths[acc]
        else:
            rec_len = domains_by_acc[acc][0].get("length")
            if rec_len:
                lengths[acc] = rec_len
            else:
                lengths[acc] = max(
                    r.get("env_to") or r.get("ali_to") or 1 for r in domains_by_acc[acc]
                )

    max_len = max(lengths.values()) if lengths else 1

    # ---------------- Collect all unique domain names for the legend ----------------
    all_domain_names = sorted(
        {r["hmm_name"] for acc in protein_list for r in domains_by_acc[acc]}
    )
    color_map = {name: _domain_color(name) for name in all_domain_names}

    # ---------------- Figure layout ----------------
    # Row height: 0.55 inches per protein + 0.4 inches top margin
    # + space for legend at the bottom (0.35 per legend row, ~4 entries/row)
    n_proteins = len(protein_list)
    legend_rows = max(1, (len(all_domain_names) + 3) // 4)
    fig_height = 0.4 + n_proteins * 0.55 + legend_rows * 0.45 + 0.8

    # Label column width: fixed fraction of figure width
    label_col_fraction = 0.30  # 30% of figure width for protein labels

    fig, ax = plt.subplots(figsize=(14, fig_height))
    ax.set_xlim(0, max_len)
    ax.set_ylim(-0.5, n_proteins - 0.5)
    ax.invert_yaxis()  # First protein at the top, matching alignment conventions

    ax.set_xlabel("Residue position", fontsize=9)
    ax.set_title(title or "Domain Architecture", fontsize=11, fontweight="bold", pad=10)

    # Hide y-axis ticks -- we'll draw labels manually as text inside the axes
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    for i, acc in enumerate(protein_list):
        y = i  # vertical centre of this protein's row
        prot_len = lengths[acc]

        # ---------------- Backbone: thin grey rectangle spanning protein length ----------------
        ax.broken_barh(
            [(0, prot_len)],
            (y - 0.06, 0.12),
            facecolors="#cccccc",
            edgecolors="none",
            zorder=1,
        )

        # Keep track of text boundaries in this row to prevent overlapping text blobs
        occupied_text_regions = []

        # ---------------- Domain blocks ----------------
        for rec in domains_by_acc[acc]:
            d_start = rec.get("ali_from", 0)
            d_end = rec.get("ali_to", 0)
            d_width = max(d_end - d_start, 1)

            # Make color slightly transparent (alpha 0.7) so stacked domains are visible
            r, g, b, _ = color_map[rec["hmm_name"]]
            color_with_alpha = (r, g, b, 0.7)

            ax.broken_barh(
                [(d_start, d_width)],
                (y - 0.22, 0.44),
                facecolors=[color_with_alpha],
                edgecolors="white",
                linewidth=0.5,
                zorder=2,
            )

            # --- Text Label Logic ---
            hmm_name = rec["hmm_name"]

            # Heuristic: Estimate text width in data coordinates
            # ~0.6% of the max axis length per character at font size 6.5
            est_text_width = len(hmm_name) * max_len * 0.006

            if d_width >= est_text_width and (
                d_width / max_len >= min_label_width_fraction
            ):

                label_x = d_start + d_width / 2
                label_start = label_x - (est_text_width / 2)
                label_end = label_x + (est_text_width / 2)

                is_overlapping = any(
                    not (label_end < occ_start or label_start > occ_end)
                    for occ_start, occ_end in occupied_text_regions
                )

                if not is_overlapping:
                    ax.text(
                        label_x,
                        y,
                        hmm_name,
                        ha="center",
                        va="center",
                        fontsize=6.5,
                        fontweight="bold",
                        color="white",
                        zorder=3,
                        clip_on=True,
                    )
                    # Record this space as occupied
                    occupied_text_regions.append((label_start, label_end))

        # ---------------- Protein label to the left of the backbone ----------------
        label_x = -max_len * 0.01
        ax.text(
            label_x,
            y,
            meta_by_acc[acc],
            ha="right",
            va="center",
            fontsize=7,
            color="#333333",
            transform=ax.transData,
        )

    # ---------------- Legend ----------------
    legend_handles = [
        mpatches.Patch(color=color_map[name], label=name) for name in all_domain_names
    ]

    y_offset = -(0.8 / max(1, n_proteins))
    ax.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, y_offset),  # below the plot
        ncol=4,
        fontsize=8,
        frameon=False,
        # title="Domain",
        title_fontsize=8,
    )

    if truncated:
        fig.text(
            0.5,
            0.01,
            f"Showing first {max_proteins} of {len(domain_records)} proteins.",
            ha="center",
            fontsize=8,
            color="grey",
            style="italic",
        )

    plt.tight_layout(pad=1.5)
    return _fig_to_buffer(fig)


# ---------------------------------------------------------------------------
# 2. Inline phylogenetic tree rendering
# ---------------------------------------------------------------------------


def render_tree(
    newick_string: str,
    ladderize: bool = True,
    label_font_size: int = 7,
    branch_color: str = "#2a5f8f",
    max_leaves_before_compact: int = 80,
) -> io.BytesIO:
    """
    Render a Newick tree string as a PNG image using Bio.Phylo + matplotlib.


    Example
    -------
    with open("myrun.mft.gt01.lg.fasttree") as f:
        nwk = f.read()
    buf = viz.render_tree(nwk)
    st.image(buf, use_container_width=True)
    """
    try:
        from Bio import Phylo
    except ImportError as e:
        raise ImportError(
            "BioPython is required for render_tree(). "
            "Install it with: pip install biopython"
        ) from e

    tree = Phylo.read(io.StringIO(newick_string), "newick")

    if ladderize:
        tree.ladderize()

    n_leaves = len(tree.get_terminals())
    compact = (
        n_leaves > max_leaves_before_compact
    )  # Check if number of leaves is more than compact limit

    # Adapt font size and figure height to leaf count
    if n_leaves > 60:
        label_font_size = max(4, label_font_size - 2)
    fig_height = max(5, min(n_leaves * 0.22, 60))  # cap at 60 inches

    fig, ax = plt.subplots(figsize=(12, fig_height))

    # Bio.Phylo.draw() uses matplotlib internally.
    # We override the default grey branch color with blue.
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=(lambda x: x.name if (x.name and not compact) else ""),
    )

    # Re-color all branch lines after drawing
    for line in ax.get_lines():
        line.set_color(branch_color)
        line.set_linewidth(0.9)

    # Re-apply font size to all text elements (leaf labels)
    for text in ax.texts:
        text.set_fontsize(label_font_size)

    ax.set_xlabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(
        f"Phylogenetic Tree  ({n_leaves} leaves)", fontsize=11, fontweight="bold"
    )

    if compact:
        ax.set_title(
            f"Phylogenetic Tree  ({n_leaves} leaves — labels hidden, too many to display)",
            fontsize=10,
            fontweight="bold",
        )

    plt.tight_layout()
    return _fig_to_buffer(fig)


# ---------------------------------------------------------------------------
# 3. Clustered presence/absence heatmap
# ---------------------------------------------------------------------------


def draw_presence_absence_heatmap(
    matrix_df,
    title: str = "Presence / Absence Heatmap",
    cmap: str = "viridis",
    cluster: bool = True,
) -> io.BytesIO:
    """
    Draw a clustered heatmap of the presence/absence matrix.

    seaborn.clustermap() reorders both rows (organisms) and columns (profiles)
    by hierarchical clustering, grouping organisms with similar profile
    complements together. This reveals evolutionary patterns — e.g. all
    Bacteria clustering together because they share a set of profiles that
    Eukaryotes lack.

    """
    try:
        import seaborn as sns
    except ImportError as e:
        raise ImportError(
            "seaborn is required for draw_presence_absence_heatmap(). "
            "Install it with: pip install seaborn"
        ) from e

    # Clustermap needs at least 2 rows and 2 columns
    if cluster and (matrix_df.shape[0] < 2 or matrix_df.shape[1] < 2):
        cluster = False

    n_rows, n_cols = matrix_df.shape
    # Scale figure so cells are approximately square regardless of matrix size
    fig_width = max(6, n_cols * 1.1 + 3)
    fig_height = max(4, n_rows * 0.5 + 2)

    if cluster:
        # clustermap manages its own figure — get it back from the returned object
        g = sns.clustermap(
            matrix_df.astype(float),
            cmap=cmap,
            figsize=(fig_width, fig_height),
            linewidths=0.4,
            linecolor="#e0e0e0",
            annot=(n_rows <= 30 and n_cols <= 20),  # show count numbers if small
            fmt=".0f",
            cbar_kws={"label": "Protein count"},
            # Leave room for the left-facing text
            cbar_pos=(0.03, 0.15, 0.02, 0.15),
            dendrogram_ratio=(0.15, 0.1),
            xticklabels=True,
            yticklabels=True,
        )

        # Flip the colorbar ticks and label to the LEFT side
        # so they point away from the heatmap and dendrogram
        g.cax.yaxis.set_ticks_position("left")
        g.cax.yaxis.set_label_position("left")

        g.ax_heatmap.set_xlabel("HMM Profile", fontsize=9)
        g.ax_heatmap.set_ylabel("Organism", fontsize=9)
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=8
        )
        g.ax_heatmap.set_yticklabels(
            g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=7
        )
        g.fig.suptitle(title, y=1.01, fontsize=11, fontweight="bold")
        buf = io.BytesIO()

        # bbox_inches="tight" will automatically expand the left edge if the text needs more room
        g.fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        buf.seek(0)
        plt.close(g.fig)
        return buf

    else:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        sns.heatmap(
            matrix_df.astype(float),
            cmap=cmap,
            ax=ax,
            linewidths=0.6,
            linecolor="#e0e0e0",
            annot=(n_rows <= 30 and n_cols <= 20),
            fmt=".0f",
            # 'pad' to push the colorbar away from the heatmap
            cbar_kws={"label": "Protein count", "shrink": 0.6, "pad": 0.04},
        )
        ax.set_xlabel("HMM Profile", fontsize=9)
        ax.set_ylabel("Organism", fontsize=9)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)
        ax.set_title(title, fontsize=11, fontweight="bold")
        plt.tight_layout()
        return _fig_to_buffer(fig)


# ---------------------------------------------------------------------------
# 4. High-resolution phylogenetic profile heatmap
# ---------------------------------------------------------------------------


def draw_highres_profile_heatmap(
    matrix_df,
    column_origin: dict = None,
    taxon_names: dict = None,
    missing_accessions=None,
    title: str = "High-Resolution Phylogenetic Profile",
    binary: bool = False,
    log_scale: bool = False,
    cluster_rows: bool = True,
    cluster_cols: bool = False,
    cmap: str = None,
    figsize=None,
) -> io.BytesIO:
    """
    Draw a clustered heatmap of the high-resolution phylogenetic profile.

    Rows are taxa, columns are Pfam-subclade pairs (e.g. "PF00041-A").
    A colored stripe above the columns groups subclades by their parent Pfam,
    making paralog-group structure visible at a glance.



    Example
    -------
    out = fetch_highres_profile("2026_01", pfam_subclade_map, binary=False)
    buf = viz.draw_highres_profile_heatmap(
        out["matrix"],
        column_origin       = out["column_origin"],
        taxon_names         = out["taxon_names"],
        missing_accessions  = out["missing_accessions"],
        log_scale           = True,
    )
    st.image(buf, use_container_width=True)
    """
    try:
        import seaborn as sns
    except ImportError as e:
        raise ImportError(
            "seaborn is required for draw_highres_profile_heatmap()."
        ) from e

    if matrix_df.empty:
        fig, ax = plt.subplots(figsize=(8, 2))
        ax.text(
            0.5,
            0.5,
            "No data to display.",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
        )
        ax.axis("off")
        return _fig_to_buffer(fig)

    n_rows, n_cols = matrix_df.shape

    # ---------------- Auto-select colormap ----------------
    if cmap is None:
        cmap = "Greys" if binary else "viridis"

    # ---------------- Build display matrix (log-transform if requested) ----------------
    if log_scale and not binary:
        display_df = np.log1p(matrix_df.astype(float))
        cbar_label = "log1p(protein count)"
    else:
        display_df = matrix_df.astype(float)
        cbar_label = "Presence (0/1)" if binary else "Protein count"

    # ---------------- Row labels: "9606  Homo sapiens" if names provided ----------------
    if taxon_names:
        new_index = [
            f"{tx}  {taxon_names.get(tx, '?')}" if taxon_names.get(tx) else str(tx)
            for tx in display_df.index
        ]
        display_df = display_df.copy()
        display_df.index = new_index

    # ---------------- Build col_colors stripe grouping subclades by their Pfam ----------------
    col_colors = None
    pfam_color_map = {}
    if column_origin:
        # One color per distinct Pfam
        unique_pfams = []
        for col in display_df.columns:
            origin = column_origin.get(col)
            if origin is None:
                continue
            pfam = origin[0]
            if pfam not in unique_pfams:
                unique_pfams.append(pfam)

        # Use a categorical palette
        palette_name = "tab10" if len(unique_pfams) <= 10 else "tab20"
        palette = sns.color_palette(palette_name, n_colors=max(len(unique_pfams), 1))
        pfam_color_map = dict(zip(unique_pfams, palette))

        col_colors = [
            pfam_color_map.get(column_origin.get(c, (None,))[0], (1, 1, 1))
            for c in display_df.columns
        ]

    # ---------------- Figure size ----------------

    max_label_len = max(
        (
            len(str(tx)) + len(str(taxon_names.get(tx, ""))) + 2
            for tx in matrix_df.index
        ),
        default=10,
    )
    label_extra = max(0, (max_label_len - 20) * 0.07)
    if figsize is None:
        fig_width = max(9, n_cols * 0.55 + 4 + label_extra)
        fig_height = max(4, n_rows * 0.45 + 2)
    else:
        fig_width, fig_height = figsize

    # ---------------- clustermap ----------------
    # Disable clustering on an axis if there's only one row/column on it.
    rc = cluster_rows and n_rows >= 2
    cc = cluster_cols and n_cols >= 2

    annot_show = n_rows <= 30 and n_cols <= 20
    # When log-scaling, color uses log1p but annotations should show
    # the original counts
    if annot_show and log_scale and not binary:
        annot = matrix_df.astype(int).copy()
        if taxon_names:
            annot.index = display_df.index  # match the relabeled rows
        annot_fmt = "d"
    elif annot_show:
        annot = True
        annot_fmt = ".0f"
    else:
        annot = False
        annot_fmt = ".0f"

    g = sns.clustermap(
        display_df,
        cmap=cmap,
        figsize=(fig_width, fig_height),
        row_cluster=rc,
        col_cluster=cc,
        col_colors=col_colors,
        linewidths=0.3,
        linecolor="#e0e0e0",
        annot=annot,
        fmt=annot_fmt,
        cbar_kws={"label": cbar_label},
        cbar_pos=(0.03, 0.15, 0.02, 0.15),
        dendrogram_ratio=(0.12, 0.08 if cc else 0.02),
        xticklabels=True,
        yticklabels=True,
    )

    # Colorbar layout adjustment (Outward vertical style)
    g.cax.yaxis.set_ticks_position("left")
    g.cax.yaxis.set_label_position("left")
    g.cax.set_ylabel(cbar_label, fontsize=8, labelpad=12, rotation=90)
    g.cax.tick_params(axis="y", labelsize=8)

    # For binary, force integer 0 / 1 ticks
    if binary:
        g.cax.set_yticks([0, 1])

    # Axis labels & tick rotations
    g.ax_heatmap.set_xlabel("Pfam · Subclade", fontsize=9)
    g.ax_heatmap.set_ylabel("Taxon", fontsize=9)
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=8
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)

    # Title
    g.fig.suptitle(title, y=1.01, fontsize=11, fontweight="bold")

    # Pfam legend
    if pfam_color_map:
        import matplotlib.patches as mpatches

        handles = [
            mpatches.Patch(color=color, label=pfam)
            for pfam, color in pfam_color_map.items()
        ]

        # Placing the legend safely below the heatmap's x-axis labels
        g.ax_heatmap.legend(
            handles=handles,
            title="Pfam",
            loc="upper center",
            bbox_to_anchor=(
                0.5,
                -0.22,
            ),  # Pushes the legend safely below the 45-degree rotated labels
            ncol=len(pfam_color_map),
            fontsize=8,
            title_fontsize=9,
            frameon=False,
        )

    # missing accessions note
    if missing_accessions:
        n_missing = len(missing_accessions)
        g.fig.text(
            0.5,
            -0.02,
            f"Note: {n_missing} accession(s) from the input trees were not "
            f"found in the database for this version and were excluded "
            f"from counts.",
            ha="center",
            va="top",
            fontsize=8,
            style="italic",
            color="#666",
        )

    buf = io.BytesIO()
    g.fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    plt.close(g.fig)
    return buf
