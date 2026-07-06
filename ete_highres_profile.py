"""ETE4 interactive viewer for a High-Res Profile matrix on the NCBI taxonomy tree.

Reads the matrix CSV produced by the 'Assemble Profile' step, extracts taxids
from row labels, builds the NCBI species tree, and renders each matrix column
(subclade) as a coloured aligned rectangle — matching the style of ete_profile.py.

Called by Flask via subprocess:
    python ete_highres_profile.py -m matrix.csv [-c colormap.txt] [-p PORT]
"""

import argparse
import csv
import sys
import seaborn as sns
from ete4 import NCBITaxa
from ete4.smartview import Layout, TextFace, RectFace
import utils

parser = argparse.ArgumentParser()
parser.add_argument(
    "-m",
    "--matrix",
    help="Profile matrix CSV from High-Res Profile step 3",
    required=True,
)
parser.add_argument(
    "-c", "--colormap", help="Branch-colour map: taxid <TAB> hex-colour, one per line"
)
parser.add_argument(
    "-p", "--port", help="Port for ete4 explore server", type=int, default=5003
)
args = parser.parse_args()

# ── Load matrix ─────────────────────────────────────────────────────────────

matrix = {}  # {taxid_str: {col: float}}
columns = []
taxid_to_label = {}

with open(args.matrix, newline="") as fh:
    reader = csv.DictReader(fh)
    fieldnames = reader.fieldnames or []
    if not fieldnames:
        print("Empty or invalid matrix CSV.")
        sys.exit(1)
    row_col = fieldnames[0]
    columns = fieldnames[1:]

    for row in reader:
        label = row[row_col]
        # Row labels are "taxid  sci_name" (two spaces) — extract taxid
        taxid = label.strip().split()[0]
        if not taxid.isdigit():
            continue
        matrix[taxid] = {}
        for col in columns:
            try:
                matrix[taxid][col] = float(row[col])
            except (ValueError, TypeError):
                matrix[taxid][col] = 0.0
        taxid_to_label[taxid] = label.strip()

if len(matrix) < 2:
    print(f"Need at least 2 taxids in the matrix (found {len(matrix)}).")
    sys.exit(1)

taxids = list(matrix.keys())

# ── Optional colormap ────────────────────────────────────────────────────────

taxid2color: dict[str, str] = {}
if args.colormap:
    with open(args.colormap) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                taxid2color[parts[0]] = parts[1]

# ── NCBI taxonomy tree ───────────────────────────────────────────────────────

ncbi = NCBITaxa()
tree = ncbi.get_topology(taxids, intermediate_nodes=False)
tree.to_ultrametric()
tax2names, tax2lineages, tax2rank = tree.annotate_ncbi_taxa(taxid_attr="taxid")

# ── Colour palette per column ────────────────────────────────────────────────

col_palette = sns.color_palette("viridis", len(columns))
col_colors = dict(zip(columns, col_palette))
col_max = {
    col: max((matrix.get(t, {}).get(col, 0) for t in taxids), default=1) or 1
    for col in columns
}

# ── Tree style ───────────────────────────────────────────────────────────────

my_tree_style = {
    "shape": "rectangular",
    "aligned-leaves": True,
    "show-popup-props": None,
    # No hz-line/vt-line here — branch_col handles all branch styling per node
}

# ── Pre-assign branch_color prop ───────────
# ── Branch colour pre-computation (same pattern as tree_from_db.py) ──────────

_PALETTE = [
    "#e6194B",
    "#3cb44b",
    "#ffe119",
    "#4363d8",
    "#f58231",
    "#911eb4",
    "#42d4f4",
    "#f032e6",
    "#bfef45",
    "#469990",
    "#dcbeff",
    "#9A6324",
    "#800000",
    "#aaffc3",
    "#808000",
    "#ffd8b1",
    "#000075",
    "#a9a9a9",
]

if taxid2color:
    # User provided a colormap: group-level taxids → color.
    # Walk each node's NCBI lineage; first matching colormap entry wins.
    # Put more-specific groups before broader ones in your colormap file.
    priority = list(taxid2color.keys())
    leaf_color = {}
    for _tid_str in taxids:
        try:
            _lin = ncbi.get_lineage(int(_tid_str)) or []
        except Exception:
            continue
        for _p in priority:
            try:
                if int(_p) in _lin:
                    leaf_color[_tid_str] = taxid2color[_p]
                    break
            except ValueError:
                continue
else:
    # No colormap: auto-assign one distinct color per leaf taxid (like tree_from_db.py)
    leaf_color = {
        tid: _PALETTE[i % len(_PALETTE)] for i, tid in enumerate(sorted(taxids))
    }

# Assign branch_color to every node; internal nodes get the shared color only
# when all their leaf descendants belong to the same color group.
for _node in tree.traverse():
    if _node.is_leaf:
        _tid_str = str(_node.props.get("taxid", ""))
        _col = leaf_color.get(_tid_str)
    else:
        _cols = {
            leaf_color.get(str(_l.props.get("taxid", ""))) for _l in _node.leaves()
        }
        _cols.discard(None)
        _col = next(iter(_cols)) if len(_cols) == 1 else None
    if _col:
        _node.add_prop("branch_color", _col)


def branch_col(node):
    # Always return a style — same pattern as tree_from_db.py.
    # Returning None lets the tree-level style override the per-node color.
    col = node.props.get("branch_color", "#333333")
    return {
        "hz-line": {"stroke": col, "stroke-width": 3},
        "vt-line": {"stroke": col, "stroke-width": 3},
    }


def name_layout(node):
    if "sci_name" in node.props:
        sci = str(node.props["sci_name"])
        return TextFace(sci, position="right", fs_min=6, fs_max=25)
    else:
        return TextFace(node.name, position="right", fs_min=6, fs_max=25)


def header_layout(tree):
    faces = []
    for n, col in enumerate(columns):
        short = col[:22] + ("…" if len(col) > 22 else "")
        faces.append(
            TextFace(
                short,
                fs_min=5,
                fs_max=15,
                rotation=-45,
                anchor=(-1, 1),
                position="header",
                column=n,
            )
        )
    return faces


def count_layout(node):
    if not node.is_leaf:
        return
    taxid = str(node.props.get("taxid", ""))
    faces = []
    for n, col in enumerate(columns):
        val = matrix.get(taxid, {}).get(col, 0)
        cmax = col_max[col]
        if val > 0:
            ratio = min(val / cmax, 1.0)
            colour = utils.color_gradient("white", col_colors[col], mix=ratio)
            text = str(int(val)) if val == int(val) else f"{val:.2f}"
        else:
            colour = "#f0f0f0"
            text = "–"
        faces.append(
            RectFace(
                wmax=160,
                hmax=80,
                style={"fill": colour},
                position="aligned",
                column=n,
                text=text,
            )
        )
    return faces


# ── Launch ───────────────────────────────────────────────────────────────────

base_layout = Layout(name="base", draw_tree=my_tree_style)
branch_layout = Layout(name="branch colours", draw_node=branch_col)
label_layout = Layout(name="scientific names", draw_node=name_layout)
hdr_layout = Layout(name="column headers", draw_tree=header_layout)
cnt_layout = Layout(name="subclade counts", draw_node=count_layout)

print(
    f"Launching ETE profile viewer on port {args.port} with {len(taxids)} taxa, {len(columns)} columns."
)

tree.explore(
    layouts=[base_layout, branch_layout, label_layout, hdr_layout, cnt_layout],
    port=args.port,
    open_browser=False,
    keep_server=True,
)
