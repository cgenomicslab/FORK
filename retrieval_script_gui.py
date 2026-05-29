import streamlit as st
import get_reference_uniprot_set_lib as uni
import subprocess
import pandas as pd
import os
import re

import viz_utils as viz
import streamlit.components.v1 as components
import interactive_tree_component as itc

st.set_page_config(page_title="UniProt Lab Manager", layout="wide")

st.title("UniProt Reference Set Manager")
st.markdown("---")

# Sidebar for Database Connection Info
with st.sidebar.expander("DB Config"):
    default_config = uni.get_db_config()
    host = st.text_input("Host", value=default_config["host"])
    user = st.text_input("User", value=default_config["user"])
    db_name = st.text_input("Database", value=default_config["database"])
    config = uni.get_db_config(host=host, user=user, database=db_name)

# ==========================================================================================================

#                                               MAIN NAVIGATION

# ==========================================================================================================


menu = [
    "Standard Retrieval",
    "HMM Search",
    "Accession Lookup",
    "Domain Coordinate Lookup",
    "Database Info",
    "Phylogenetic Tree",
    "GO → Domain Profiles",
    "Presence/Absence & Drill-down",
]
choice = st.selectbox("Select Functionality", menu)



# ==========================================================================================================

#                                               STANDARD RETRIEVAL TAB

# ==========================================================================================================

if choice == "Standard Retrieval":
    st.header("Filtered Sequence Retrieval")
    col1, col2 = st.columns(2)
    with col1:
        ver = st.text_input("UniProt Version", value="2026_01")
        tax = st.text_input("Taxonomy IDs (comma separated, e.g. 9606, 10090)")
        tax_file = st.file_uploader("Or upload txt (one ID per line)", type=["txt"], key="sr_tax_file")
    with col2:
        proteome = st.text_input("Proteome ID (Optional)")
        go_id = st.text_input("GO ID (Optional)")
        pfam_id = st.text_input("Pfam ID (Optional)")

    if st.button("Fetch Sequences"):
        tax_ids = [int(l.strip()) for l in tax_file.read().decode().splitlines() if l.strip()] if tax_file else ([int(t.strip()) for t in tax.split(",") if t.strip()] if tax else None)
        records = uni.fetch_sequences(
            ver, tax_ids, proteome, go_id, pfam_id, db_config=config
        )
        if records:
            st.success(f"Retrieved {len(records)} sequences.")
            st.dataframe(pd.DataFrame(records).drop(columns=["sequence"]))
            fasta_str = uni.fetch_fasta_string(
                ver, tax_ids, proteome, go_id, pfam_id, db_config=config
            )
            st.download_button(
                "Download FASTA", fasta_str, file_name=f"uniprot_{ver}.fasta"
            )
        else:
            st.warning("No records found.")


# ==========================================================================================================

#                                               HMM SEARCH TAB

# ==========================================================================================================

elif choice == "HMM Search":
    st.header("HMM Hit Retrieval")
    ver = st.text_input("UniProt Version", value="2026_01")
    hmm_query = st.text_input("HMM Name or Accession (e.g. Homeodomain or PF00046)")
    use_evalue = st.checkbox("Apply E-value cutoff")
    eval_cutoff = (
        st.number_input("E-value Cutoff", value=1e-5, format="%.1e")
        if use_evalue
        else None
    )
    tax = st.text_input("Taxonomy Filter (Optional, comma separated)")
    tax_file = st.file_uploader("Or upload TXT (one ID per line)", type=["txt"], key="hmm_tax_file")

    if st.button("Run HMM Search"):
        tax_ids = [int(l.strip()) for l in tax_file.read().decode().splitlines() if l.strip()] if tax_file else ([int(t.strip()) for t in tax.split(",") if t.strip()] if tax else None)
        records = uni.fetch_sequences_by_hmm_hit(
            ver, hmm_query, eval_cutoff, tax_ids, db_config=config
        )
        if records:
            st.success(f"Found {len(records)} HMM hits.")
            st.dataframe(pd.DataFrame(records).drop(columns=["sequence"]))
            fasta_str = uni.fetch_fasta_string_by_hmm_hit(
                ver, hmm_query, eval_cutoff, tax_ids, db_config=config
            )
            st.download_button(
                "Download FASTA", fasta_str, file_name=f"hmm_{hmm_query}.fasta"
            )
        else:
            st.warning("No records found.")


# ==========================================================================================================

#                                               ACCESSION LOOKUP TAB

# ==========================================================================================================

elif choice == "Accession Lookup":
    st.header("Batch Accession or Protein name Retrieval")
    ver = st.text_input("UniProt Version", value="2026_01")
    acc_input = st.text_area("Paste Accessions or Protein names (one per line or space separated)")

    if st.button("Get Sequences"):
        acc_list = acc_input.replace(",", " ").split()
        records = uni.fetch_sequences_by_accession(ver, acc_list, db_config=config)
        if records:
            st.dataframe(pd.DataFrame(records))
            with uni.UniProtRetriever(config) as db:
                fasta_str = db.to_fasta_string(records)
            st.download_button(
                "Download FASTA", fasta_str, file_name="accessions.fasta"
            )
        else:
            st.warning("No records found.")

# ==========================================================================================================

#                                              DOMAIN COORDINATES LOOKUP TAB

# ==========================================================================================================

elif choice == "Domain Coordinate Lookup":
    st.header("Protein Domain Architecture")
    ver = st.text_input("UniProt Version", value="2026_01")
    acc_input = st.text_area("Enter Accessions or Protein names")
    use_evalue = st.checkbox("Apply E-value cutoff")
    eval_cutoff = (
        st.number_input("E-value Cutoff", value=1e-5, format="%.1e")
        if use_evalue
        else None
    )

    # 1. Fetch data and save it to session_state
    if st.button("Analyze Domains"):
        acc_list = acc_input.replace(",", " ").split()
        if not acc_list:
            st.warning("Please enter at least one accession.")
        else:
            with st.spinner("Fetching domain data..."):
                domains = uni.fetch_domains_by_accession(
                    ver, acc_list, eval_cutoff, db_config=config
                )

            if domains:
                # Store results so they survive
                st.session_state["dc_domains"] = domains
                st.session_state["dc_acc_list"] = acc_list
                # Reset the image flag so it hides until explicitly drawn
                st.session_state["show_dc_arch"] = False
            else:
                st.warning("No domains found.")
                st.session_state.pop("dc_domains", None)

    # 2. Show UI if data exists in memory
    if st.session_state.get("dc_domains"):
        domains = st.session_state["dc_domains"]
        acc_list = st.session_state["dc_acc_list"]

        df = pd.DataFrame(domains)
        st.dataframe(df)
        st.download_button(
            "Download CSV", df.to_csv(index=False), "domains.csv", key="dc_csv_dl"
        )

        # ──----------------- Architecture diagram ────────────────────────────────────
        # The domain records already in memory are passed directly to
        # viz_utils — no additional DB call needed.
        # We wrap it in a button so the figure is only rendered when
        # the user explicitly asks for it (rendering can be slow for
        # large protein sets).
        st.markdown("---")

        # 3. Use the button just to flip a toggle in session_state!
        if st.button("Draw Domain Architectures", key="dc_draw_btn"):
            st.session_state["show_dc_arch"] = True

        # 4.
        # This keeps the image alive so the Download button works.
        if st.session_state.get("show_dc_arch"):
            with st.spinner("Rendering architecture diagram..."):
                buf = viz.draw_domain_architecture(
                    domains,
                    title=f"Domain Architecture — {', '.join(acc_list[:3])}"
                    + (" …" if len(acc_list) > 3 else ""),
                )

            st.image(buf, use_container_width=True)

            st.download_button(
                "Download Architecture PNG",
                buf,
                file_name="domain_architecture.png",
                mime="image/png",
                key="dc_arch_png",
            )


# ==========================================================================================================

#                                               DATABASE INFO TAB

# ==========================================================================================================
elif choice == "Database Info":
    st.header("Database Status")
    if st.button("List Versions & Stats"):
        with uni.UniProtRetriever(config) as db:
            versions = db.list_available_versions()
            st.table(versions)
            

# ==========================================================================================================

#                                              PHYLOGENETIC TREE TAB

# ==========================================================================================================

elif choice == "Phylogenetic Tree":
    st.header("Phylogenetic Tree Builder")
    col1, col2 = st.columns(2)
    with col1:
        ver = st.text_input("UniProt Version", value="2026_01")
        pfam = st.text_input(
            "Pfam IDs or HMM names (comma separated)",
            placeholder="e.g. Homeodomain, PF00001,PF00002",
        )
        tax = st.text_input("Taxonomy IDs (comma separated, optional)")
        tax_file = st.file_uploader("Or upload txt (one ID per line)", type=["txt"], key="tree_tax_file")
        exclude_tax = st.text_input("Exclude Taxonomy IDs (comma separated, optional)")
        exclude_tax_file = st.file_uploader("Or upload exclude txt", type=["txt"], key="tree_excl_file")
        prefix = st.text_input("Output Prefix", placeholder="e.g. myrun")
        output_dir = st.text_input("Output Directory", placeholder="e.g. /home/user/results (leave empty for current dir)")

    with col2:
        aln = st.selectbox("Alignment Tool", ["mafft", "einsi", "clustalo"])
        ml = st.selectbox("Tree Method", ["fasttree", "iqtree"])
        cpu = st.text_input("Threads", value="32")
        use_evalue = st.checkbox("Apply E-value cutoff")
        evalue = (
            st.number_input("E-value Cutoff", value=1e-5, format="%.1e")
            if use_evalue
            else None
        )
        no_ncbi = st.checkbox("Skip NCBI annotation (faster)", value=True)

        # ETE4 Server
        use_ete4 = st.checkbox(
            "Start ETE4 Interactive Server (Alternative Viewer)", value=False
        )
        ete4_port = st.number_input("ETE4 Port", value=5001) if use_ete4 else 5001

    if st.button("Run Full Tree Pipeline"):
        if not pfam or not prefix:
            st.warning(
                "Please fill in at least one Pfam ID or HMM name and Output Prefix."
            )
        else:
            cmd = [
                "python",
                "tree_from_db.py",
                "--pfam",
                pfam.replace(" ", ""),
                "--version",
                ver,
                "--prefix",
                prefix,
                "--aln",
                aln,
                "--ml",
                ml,
                "--cpu",
                cpu,
                
            ]
            
            if output_dir:
                cmd += ["--output_dir", output_dir.strip()]
                
            # tax
            if tax_file:
                tax_str = ",".join(l.strip() for l in tax_file.read().decode().splitlines() if l.strip())
            elif tax:
                tax_str = tax.replace(" ", "")
            else:
                tax_str = None
            if tax_str:
                cmd += ["--taxids", tax_str]

            # exclude_tax
            if exclude_tax_file:
                excl_str = ",".join(l.strip() for l in exclude_tax_file.read().decode().splitlines() if l.strip())
            elif exclude_tax:
                excl_str = exclude_tax.replace(" ", "")
            else:
                excl_str = None
            if excl_str:
                cmd += ["--exclude_taxids", excl_str]
                
            if evalue is not None:
                cmd += ["--evalue", str(evalue)]
            if no_ncbi:
                cmd += ["--no_ncbi"]

            # ── PATH A: USE ETE4 BACKGROUND SERVER ──
            if use_ete4:
                cmd += ["--port", str(int(ete4_port))]
                with st.spinner(
                    f"Preparing ETE4 on port {ete4_port}. Alignment & Tree building may take 1-2 minutes..."
                ):

                    # 1. Automatically kill any old zombie server holding this port

                    os.system(f"fuser -k {ete4_port}/tcp >/dev/null 2>&1")

                    # 2. Start the pipeline in the background
                    proc = subprocess.Popen(cmd)

                    # 3. Give it 30 seconds to run MAFFT/FastTree before showing the frame
                    import time

                    time.sleep(30)

                st.success(f"ETE4 Server launched! Bound to port {ete4_port}.")
                st.info(
                    "If the frame below says 'Cannot connect', the alignment is still running! Just wait 30 seconds and click the blue (Try Again) button in the frame."
                )

                # Tell session state to show ETE4, not D3
                st.session_state["tree_ready"] = True
                st.session_state["tree_prefix"] = prefix
                st.session_state["viewer_mode"] = "ete4"
                st.session_state["ete4_port"] = ete4_port

            # ── PATH B: USE D3 VIEWER (DEFAULT, NO SERVER) ──
            else:
                cmd += ["--no_explore"]  # Skip the ETE4 server block

                with st.spinner(
                    "Building alignment and tree... Check your terminal for live progress!"
                ):
                    try:
                        subprocess.run(cmd, check=True)
                    except subprocess.CalledProcessError as e:
                        st.error(
                            f"Pipeline crashed (Exit status {e.returncode}). Check your terminal for the exact error message."
                        )
                        st.stop()

                if ml == "fasttree":
                    tree_file = f"{prefix}.mft.gt01.lg.fasttree"
                elif ml == "iqtree":
                    tree_file = f"{prefix}.mft.gt01.treefile"

                if os.path.isfile(tree_file):
                    with open(tree_file) as f:
                        st.session_state["tree_data"] = f.read()
                    st.session_state["tree_prefix"] = prefix

                    itol_colors_path = f"{tree_file}.itol_colors.txt"
                    if os.path.isfile(itol_colors_path):
                        with open(itol_colors_path) as f:
                            st.session_state["itol_data"] = f.read()

                    # Tell session state to show D3, not ETE4
                    st.session_state["tree_ready"] = True
                    st.session_state["viewer_mode"] = "d3"
                    st.session_state["tree_method"] = ml
                else:
                    st.error("Pipeline finished but the tree file was not found.")

    # ==========================================================
    # Display the correct viewer based on what was just run
    # ==========================================================
    if st.session_state.get("tree_ready"):
        p = st.session_state["tree_prefix"]
        st.markdown("---")

        # --- SHOW ETE4 VIEWER ---
        if st.session_state.get("viewer_mode") == "ete4":
            st.subheader("ETE4 Interactive Explorer")
            port = st.session_state["ete4_port"]
            st.caption(f"Connected to ETE4 server on port {port}")
            components.iframe(
                f"http://localhost:{port}", width=1200, height=800, scrolling=True
            )

        # --- SHOW D3 VIEWER ---
        elif st.session_state.get("viewer_mode") == "d3":
            st.subheader("Interactive Tree Preview")
            st.caption(
                "Scroll to zoom, click nodes to collapse/expand. Hover over leaves for details."
            )

            if st.session_state.get("tree_method") == "fasttree":
                itol_path = f"{p}.mft.gt01.lg.fasttree.itol_colors.txt"
            else:
                itol_path = f"{p}.mft.gt01.treefile.itol_colors.txt"

            leaf_colors = None
            if os.path.isfile(itol_path):
                leaf_colors = itc.parse_itol_colors(itol_path)

            tree_html = itc.build_tree_html(
                newick_str=st.session_state["tree_data"],
                leaf_colors=leaf_colors,
                title=f"Phylogeny: {p}",
                height=800,
            )

            # Revert back to components.html which safely sandboxes the JS
            components.html(tree_html, height=800, scrolling=True)

            # ── Static PNG download ──────────────────────────────────
            if st.button("Render static PNG", key="tree_png_btn"):
                with st.spinner("Rendering..."):
                    png_buf = viz.render_tree(st.session_state["tree_data"])
                st.image(png_buf, use_container_width=True)
                png_buf.seek(0)
                st.download_button(
                    "Download Tree PNG",
                    png_buf,
                    file_name=f"{p}_tree.png",
                    mime="image/png",
                    key="tree_png_dl",
                )

        st.markdown("---")

        # ── Newick and iTOL downloads ────────────────────────────────────
        st.caption("For publication figures, upload the Newick to iTOL:")
        dl_col1, dl_col2 = st.columns(2)
        with dl_col1:
            if "tree_data" in st.session_state:
                st.download_button(
                    "Download Tree (Newick)",
                    st.session_state["tree_data"],
                    file_name=f"{p}.nwk",
                    key="tree_nwk_dl",
                )
        with dl_col2:
            if "itol_data" in st.session_state:
                st.download_button(
                    "Download iTOL Colors",
                    st.session_state["itol_data"],
                    file_name=f"{p}_itol_styles.txt",
                    key="tree_itol_dl",
                )
        st.success("Done! Upload the .nwk file to https://itol.embl.de")


# ==========================================================================================================

#                                          GO -> DOMAIN PROFILES TAB

# ==========================================================================================================

elif choice == "GO → Domain Profiles":
    st.header("GO Term → HMM Domain Profiles")
    ver = st.text_input("UniProt Version", value="2026_01")
    go_input = st.text_input("GO Term (e.g. GO:0005634)")
    use_evalue = st.checkbox("Apply E-value cutoff")
    eval_cutoff = (
        st.number_input("E-value Cutoff", value=1e-5, format="%.1e")
        if use_evalue
        else None
    )

    if st.button("Find Domain Profiles"):
        results = uni.fetch_domains_by_go(ver, go_input, eval_cutoff, db_config=config)
        if results:
            df = pd.DataFrame(results)
            st.success(f"Found {len(results)} domain profiles in GO term {go_input}")
            st.dataframe(df)
            st.download_button(
                "Download CSV",
                df.to_csv(index=False),
                file_name=f"go_{go_input.replace(':','_')}_domains.csv",
            )
        else:
            st.warning("No domain profiles found for this GO term.")


# ==========================================================================================================

#                                       PRESENCE ABSENCE & DRILL DOWN TAB

# ==========================================================================================================

elif choice == "Presence/Absence & Drill-down":
    st.header("Presence / Absence Matrix + Functional Drill-down")

    # -----------------------------------------------------------------
    # EXPLAINER
    # -----------------------------------------------------------------
    with st.expander("How this works", expanded=False):
        st.markdown("""
        **Step 1 — Build the matrix**
        Enter one or more Pfam names or accessions and a set of taxonomy IDs.
        The result is a table: rows = organisms, columns = profiles, cells = how many
        proteins in that organism carry that profile.

        **Step 2 — Drill down into a cell**
        Pick a cell (organism + profile combination) and click *Drill Down*.
        Two views open:
        - **Sub-profiles:** which more specific HMMs (TIGRFAM, sub-Pfam, …) are
        enriched in *this* subset of proteins? This surfaces enzyme-level resolution
        (e.g. Histamine DC vs Tyramine DC within a decarboxylase family).
        - **Domain architectures:** what domain combinations do these proteins carry?
        Useful when enzyme identity is defined by multi-domain context rather than
        a single specific HMM.
        """)

    st.markdown("---")

    # =================================================================
    # SECTION 1 — Inputs for the matrix query
    # =================================================================
    # We use two columns to keep the form compact, matching the style
    # already used in the Phylogenetic Tree tab.
    col1, col2 = st.columns(2)

    with col1:
        pa_ver = st.text_input("UniProt Version", value="2026_01", key="pa_ver")
        pa_pfam_input = st.text_input(
            "Pfam Names/Accessions (comma separated)",
            placeholder="e.g. Homeodomain, PF00046, PBC",
            key="pa_pfam_input",
        )
        pa_tax_input = st.text_input(
            "Taxonomy IDs (comma separated, leave empty for all)",
            placeholder="e.g. 9606, 10090, 7227",
            key="pa_tax_input",
        )
        pa_tax_file = st.file_uploader("Or upload TXT (one ID per line)", type=["txt"], key="pa_tax_file")

    with col2:
        pa_use_evalue = st.checkbox("Apply E-value cutoff", key="pa_use_evalue")
        pa_evalue = (
            st.number_input(
                "E-value Cutoff", value=1e-5, format="%.1e", key="pa_evalue"
            )
            if pa_use_evalue
            else None
        )

    # -----------------------------------------------------------------
    #                       Build Matrix button
    # -----------------------------------------------------------------
    # When clicked, we run the matrix query and store the result in
    # session_state so it survives subsequent interactions (drill-down
    # clicks, selectbox changes, etc.) without re-querying the DB.
    # -----------------------------------------------------------------
    if st.button("Build Matrix", key="pa_build"):
        if not pa_pfam_input.strip():
            st.warning("Please enter at least one Pfam name or accession.")
        else:
            pfam_queries = [q.strip() for q in pa_pfam_input.split(",") if q.strip()]
            if pa_tax_file:
                tax_ids = [int(l.strip()) for l in pa_tax_file.read().decode().splitlines() if l.strip()]
            elif pa_tax_input.strip():
                tax_ids = [int(t.strip()) for t in pa_tax_input.split(",") if t.strip()]
            else:
                tax_ids = None


            with st.spinner("Querying database..."):
                rows = uni.fetch_presence_absence_matrix(
                    pa_ver,
                    pfam_queries,
                    tax_ids,
                    pa_evalue,
                    db_config=config,
                )

            # Store raw rows — we need them to populate the drill-down
            # selectboxes later, even after the matrix is built.
            st.session_state["pa_rows"] = rows
            # Clear any previous drill-down results so stale data from a
            # previous query doesn't appear below the new matrix.
            st.session_state.pop("pa_drill_results", None)
            st.session_state.pop("pa_arch_results", None)

    # =================================================================
    # SECTION 2 — Display the matrix (if it exists in session_state)
    # =================================================================
    # We check session_state, not a local variable, because Streamlit
    # reruns from the top on every interaction. The matrix must already
    # be stored to still be visible when the user changes a selectbox.
    # =================================================================
    if st.session_state.get("pa_rows"):
        rows = st.session_state["pa_rows"]

        if not rows:
            st.warning("No hits found for these Pfam queries and taxa.")
        else:
            st.success(
                f"Matrix built — {len(rows)} (taxon × profile) combinations found."
            )

            # ── Build a human-readable taxon label ──────────────────
            # Prefer "9606 · Homo sapiens" over bare "9606" where the
            # proteomes table provided a scientific name.
            def taxon_label(taxon_id, scientific_name):
                if scientific_name:
                    return f"{taxon_id} · {scientific_name}"
                return str(taxon_id)

            df_flat = pd.DataFrame(rows)

            # Add the label column for display; keep taxon_id for queries.
            df_flat["taxon_label"] = df_flat.apply(
                lambda r: taxon_label(r["taxon_id"], r.get("scientific_name")), axis=1
            )

            # ─────────────────── Pivot into matrix  ───────────────────────────────
            # Rows = organisms, columns = HMM profiles, values = protein counts.
            # fill_value=0 makes absent (taxon, profile) pairs explicit zeros
            # rather than NaN — important for the color gradient.
            matrix_df = df_flat.pivot_table(
                index="taxon_label",
                columns="hmm_name",
                values="protein_count",
                aggfunc="sum",  # sums across multiple hmm_accession versions
                fill_value=0,
            )
            matrix_df.index.name = "Organism (taxon_id · name)"
            matrix_df.columns.name = None

            # ───────────────────── Color the matrix ────────────────────────────────────
            # background_gradient applies a colour scale across the
            # entire table (axis=None), so we compare all cells together.
            # Zeros stay near-white; high counts go dark orange/red.
            # The format call removes decimal places (counts are integers).
            styled_matrix = matrix_df.style.background_gradient(
                cmap="viridis", axis=None
            ).format("{:.0f}")

            st.subheader("Step 1 — Presence / Absence Matrix")
            st.caption(
                "Cell values = distinct proteins. "
                "Colour intensity reflects count (white = 0, dark = high)."
            )
            st.dataframe(styled_matrix, width="stretch")

            # Download the raw (unstyled) matrix as CSV
            st.download_button(
                "Download Matrix CSV",
                matrix_df.reset_index().to_csv(index=False),
                file_name="presence_absence_matrix.csv",
            )

            # ──--------------- Clustered heatmap ────────────────────────────────────────
            # The clustermap below reorders rows and columns by similarity,
            #  useful when you have many organisms/profiles and
            # want to see which groups of organisms share a profile complement.
            # on-demand button (not automatic) because
            # seaborn clustermap takes a moment on large matrices and we
            # don't want it blocking the page on every rerun.
            if matrix_df.shape[0] >= 2 and matrix_df.shape[1] >= 1:
                if st.button("Draw Clustered Heatmap", key="pa_heatmap"):
                    with st.spinner("Clustering and rendering heatmap..."):
                        heatmap_buf = viz.draw_presence_absence_heatmap(
                            matrix_df,
                            title="Presence / Absence — Clustered Heatmap",
                            cluster=(matrix_df.shape[1] >= 2),
                            cmap="viridis",
                        )
                    st.image(heatmap_buf, use_container_width=True)
                    heatmap_buf.seek(0)  # ← reset after st.image() consumed it
                    st.download_button(
                        "Download Heatmap PNG",
                        heatmap_buf,
                        file_name="presence_absence_heatmap.png",
                        mime="image/png",
                        key="pa_heatmap_dl",
                    )

            # =============================================================
            # SECTION 3 — Cell selector for drill-down
            # =============================================================
            # Two selectboxes (taxon + profile) populated from
            # the matrix. The user picks the cell they want to investigate,
            # then clicks "Drill Down".
            # =============================================================
            st.markdown("---")
            st.subheader("Step 2 — Drill Down into a Cell")

            # Build (label → taxon_id) and (label → hmm_name) maps
            # for the two selectboxes.
            taxon_options = (
                df_flat[["taxon_label", "taxon_id"]]
                .drop_duplicates()
                .sort_values("taxon_label")
            )
            profile_options = sorted(df_flat["hmm_name"].unique().tolist())

            col3, col4 = st.columns(2)
            with col3:
                selected_taxon_label = st.selectbox(
                    "Select Organism",
                    taxon_options["taxon_label"].tolist(),
                    key="pa_sel_taxon",
                )
            with col4:
                selected_profile = st.selectbox(
                    "Select Profile",
                    profile_options,
                    key="pa_sel_profile",
                )

            # Resolve the label back to an integer taxon_id for the DB call
            selected_taxon_id = int(
                taxon_options.loc[
                    taxon_options["taxon_label"] == selected_taxon_label, "taxon_id"
                ].values[0]
            )

            # Show the count for the selected cell as a quick sanity check
            cell_count = int(
                matrix_df.loc[selected_taxon_label, selected_profile]
                if selected_profile in matrix_df.columns
                and selected_taxon_label in matrix_df.index
                else 0
            )
            st.caption(
                f"Selected cell: **{selected_taxon_label}** × **{selected_profile}** "
                f"→ {cell_count} protein(s)"
            )

            if cell_count == 0:
                st.info("This cell is zero — no proteins to drill into.")
            else:
                if st.button("Drill Down", key="pa_drill"):
                    with st.spinner("Fetching accessions and computing drill-down..."):

                        # ── Bridge query: get accessions for this cell ──
           
                        cell_records = uni.fetch_accessions_for_cell(
                            pa_ver,
                            selected_profile,
                            selected_taxon_id,
                            pa_evalue,
                            db_config=config,
                        )
                        accessions = [r["accession"] for r in cell_records]

                        # ───────────────── Path A: sub-profile hits ────────────────────
                        # We exclude the original query profile from the
                        # results so it doesn't dominate the view 
                        subprofile_rows = uni.fetch_subprofile_hits(
                            pa_ver,
                            accessions,
                            pa_evalue,
                            exclude_queries=[selected_profile],
                            db_config=config,
                        )

                        # ───────────────────Path B: domain architectures ────────────────
                        arch_rows = uni.fetch_domain_architectures(
                            pa_ver,
                            accessions,
                            pa_evalue,
                            collapse_repeats=True,
                            db_config=config,
                        )

                        # Store everything in session_state.
                        # Also store context strings so the display section
                        # below can show which cell was drilled.
                        st.session_state["pa_drill_results"] = subprofile_rows
                        st.session_state["pa_arch_results"] = arch_rows
                        st.session_state["pa_drill_context"] = {
                            "taxon_label": selected_taxon_label,
                            "profile": selected_profile,
                            "n_proteins": len(accessions),
                        }

    # =================================================================
    # SECTION 4 — Display drill-down results (if they exist)
    # =================================================================
    # Again read from session_state, not local variables, so the
    # results stay visible even if the user later changes a selectbox
    # without clicking "Drill Down" again.
    # =================================================================
    if st.session_state.get("pa_drill_results") is not None:
        ctx = st.session_state["pa_drill_context"]
        st.markdown("---")
        st.subheader(
            f"Drill-down results: {ctx['profile']} in {ctx['taxon_label']} "
            f"({ctx['n_proteins']} proteins)"
        )

        # Use st.tabs so both views are co-equal and immediately switchable.
        # Tab names are self-explanatory to a lab member who hasn't read docs.
        tab_a, tab_b = st.tabs(
            [
                "Sub-profiles (deeper HMM resolution)",
                "Domain Architectures (co-occurrence)",
            ]
        )

        # ───────────────────── Tab A: sub-profile enrichment ───────────────────────────────
        with tab_a:
            st.caption(
                "All HMM profiles found on these proteins, ranked by how many "
                "proteins carry them. The original query profile has been removed. "
                "TIGRFAM and SUPERFAMILY entries here indicate enzyme-level specificity."
            )
            subprofile_rows = st.session_state["pa_drill_results"]
            if subprofile_rows:
                df_sub = pd.DataFrame(subprofile_rows)

                # Reorder columns for readability: most informative first.
                col_order = [
                    "hmm_name",
                    "hmm_type",
                    "protein_count",
                    "coverage",
                    "best_evalue",
                    "best_score",
                    "hmm_accession",
                ]
                df_sub = df_sub[[c for c in col_order if c in df_sub.columns]]

                # Format coverage as percentage for readability.
                df_sub["coverage"] = (df_sub["coverage"] * 100).round(1).astype(
                    str
                ) + "%"

                # Highlight by protein_count so the dominant sub-profiles
                # stand out visually — same cmap as the matrix for consistency.
                styled_sub = df_sub.style.background_gradient(
                    subset=["protein_count"], cmap="viridis"
                ).format({"best_evalue": "{:.2e}", "best_score": "{:.1f}"})
                st.dataframe(styled_sub, width="stretch")
                st.download_button(
                    "Download Sub-profiles CSV",
                    df_sub.to_csv(index=False),
                    file_name=f"subprofiles_{ctx['profile']}_{ctx['taxon_label'].split(' ')[0]}.csv",
                )
            else:
                st.info("No additional profiles found on these proteins.")

        # ───────────────────── Tab B: domain architectures ─────────────────────────────────
        with tab_b:
            st.caption(
                "Domain architecture patterns across these proteins, ranked by "
                "frequency. Architecture = ordered domain names (left to right "
                "on the protein). Repeated domains are collapsed by default."
            )
            arch_rows = st.session_state["pa_arch_results"]
            if arch_rows:
                df_arch = pd.DataFrame(arch_rows)

                # example_accessions is a list; show it as a comma-joined
                # string so it displays cleanly in the dataframe.
                if "example_accessions" in df_arch.columns:
                    df_arch["example_accessions"] = df_arch["example_accessions"].apply(
                        lambda x: ", ".join(x) if isinstance(x, list) else x
                    )

                col_order = [
                    "architecture",
                    "protein_count",
                    "coverage",
                    "arch_accessions",
                    "example_accessions",
                ]
                df_arch = df_arch[[c for c in col_order if c in df_arch.columns]]
                df_arch["coverage"] = (df_arch["coverage"] * 100).round(1).astype(
                    str
                ) + "%"

                styled_arch = df_arch.style.background_gradient(
                    subset=["protein_count"], cmap="viridis"
                )
                st.dataframe(styled_arch, width="stretch")
                st.download_button(
                    "Download Architectures CSV",
                    df_arch.to_csv(index=False),
                    file_name=f"architectures_{ctx['profile']}_{ctx['taxon_label'].split(' ')[0]}.csv",
                )

                # ───────────────────── Architecture diagram for a selected pattern ─────────────
            # Tab B already shows a ranked table of architecture patterns.
            # Functionality for the user to pick one pattern and visualize the
            # example proteins from that pattern as a domain diagram.
            # We can only draw proteins for which we have domain records,
            # so we do a small targeted DB call using the example_accessions
            # stored in the architecture result.
            if arch_rows:
                arch_pattern_options = [r["architecture"] for r in arch_rows]
                selected_arch = st.selectbox(
                    "Draw architecture diagram for pattern:",
                    arch_pattern_options,
                    key="pa_arch_sel",
                )

                if st.button("Draw Selected Architecture", key="pa_arch_draw"):
                    # Find the example accessions for the chosen pattern
                    selected_row = next(
                        r for r in arch_rows if r["architecture"] == selected_arch
                    )
                    example_accs = selected_row.get("example_accessions", [])

                    if not example_accs:
                        st.info("No example accessions available for this pattern.")
                    else:
                        with st.spinner("Fetching domain data and rendering..."):
                            # Fetch full domain records for these example proteins
                            draw_domains = uni.fetch_domains_by_accession(
                                pa_ver,
                                example_accs,
                                pa_evalue,
                                db_config=config,
                            )
                            buf = viz.draw_domain_architecture(
                                draw_domains,
                                title=f"Architecture: {selected_arch}",
                            )
                        st.image(buf, width="stretch")
                        st.download_button(
                            "Download Architecture PNG",
                            buf,
                            file_name=f"arch_{selected_arch[:40].replace('+','_')}.png",
                            mime="image/png",
                            key="pa_arch_png",
                        )

            else:
                st.info("No domain architecture data found for these proteins.")
