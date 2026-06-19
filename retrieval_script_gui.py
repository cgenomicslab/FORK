import streamlit as st
import get_reference_uniprot_set_lib as uni
import subprocess
import pandas as pd
import os
import signal
import socket
import time
import io
import sys
from pathlib import Path


import viz_utils as viz
import streamlit.components.v1 as components
import interactive_tree_component as itc

import subclade_partition as sp
import tree_builder as tb
from get_reference_uniprot_set_lib import fetch_highres_profile

from PIL import Image

Image.MAX_IMAGE_PIXELS = (
    None  # to avoid too many pixels error from the ete tree display
)

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
    "Extract Downloaded Branch",
    "High-Resolution Phylogenetic Profile",
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
        tax_file = st.file_uploader(
            "Or upload txt (one ID per line)", type=["txt"], key="sr_tax_file"
        )
    with col2:
        proteome = st.text_input("Proteome ID (Optional)")
        go_id = st.text_input("GO ID (Optional)")
        pfam_id = st.text_input("Pfam ID (Optional)")

    if st.button("Fetch Sequences", type="primary"):
        tax_ids = (
            [
                int(l.strip())
                for l in tax_file.getvalue().decode().splitlines()
                if l.strip()
            ]
            if tax_file
            else (
                [int(t.strip()) for t in tax.split(",") if t.strip()] if tax else None
            )
        )
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
    tax_file = st.file_uploader(
        "Or upload txt (one ID per line)", type=["txt"], key="hmm_tax_file"
    )
    if tax_file is not None:
        st.session_state["hmm_tax_content"] = tax_file.getvalue().decode()

    if st.button("Run HMM Search", type="primary"):
        if (
            "hmm_tax_content" in st.session_state
            and st.session_state["hmm_tax_content"].strip()
        ):
            tax_ids = ",".join(
                l.strip()
                for l in st.session_state["hmm_tax_content"].splitlines()
                if l.strip()
            )
        elif tax:
            tax_ids = tax.replace(" ", "")
        else:
            tax_ids = None
        tax_ids = [int(t) for t in tax_ids.split(",") if t.strip()] if tax_ids else None

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
    acc_input = st.text_area(
        "Paste Accessions or Protein names (one per line or space separated)"
    )

    if st.button("Get Sequences", type="primary"):
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
    if st.button("Analyze Domains", type="primary"):
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
        if st.button("Draw Domain Architectures", type="primary", key="dc_draw_btn"):
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
    if st.button("List Versions & Stats", type="primary"):
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
        tax_file = st.file_uploader(
            "Or upload txt (one ID per line)", type=["txt"], key="tree_tax_file"
        )
        if tax_file is not None:
            st.session_state["tree_tax_content"] = tax_file.getvalue().decode()
        exclude_tax = st.text_input("Exclude Taxonomy IDs (comma separated, optional)")
        exclude_tax_file = st.file_uploader(
            "Or upload exclude txt", type=["txt"], key="tree_excl_file"
        )
        if exclude_tax_file is not None:
            st.session_state["tree_excl_content"] = exclude_tax_file.getvalue().decode()
        prefix = st.text_input("Output Prefix", placeholder="e.g. myrun")
        output_dir = st.text_input(
            "Output Directory",
            placeholder="e.g. /home/user/results (leave empty for current dir)",
        )

    with col2:
        aln = st.selectbox("Alignment Tool", ["mafft", "einsi", "clustalo"])
        ml = st.selectbox("Tree Method", ["fasttree", "iqtree"])
        trimal_th = st.text_input("Set trimming threshold", value="0.01")
        cpu = st.text_input("Threads", value="4")
        use_evalue = st.checkbox("Apply E-value cutoff")
        evalue = (
            st.number_input("E-value Cutoff", value=1e-5, format="%.1e")
            if use_evalue
            else None
        )
        no_ncbi = st.checkbox("Skip NCBI annotation (faster)", value=True)

        pfam_source = st.selectbox(
            "Pfam source",
            ["hmmsearch", "uniprot"],
            help="hmmsearch = local Pfam-A HMM results; uniprot = UniProt's own "
            "protein_pfam assignments (no e-value/exclude filter).",
        )
        pfam_logic = st.selectbox(
            "Multi-Pfam logic",
            ["or", "and"],
            help="With >1 Pfam: or = any Pfam; and = only proteins with all Pfams.",
        )
        color_by = st.selectbox(
            "Colour branches by",
            ["taxon", "pfam"],
            help="taxon = by lineage; pfam = by which queried Pfam(s) each protein carries.",
        )
        local_fasta = st.text_input(
            "Local FASTA path (optional)",
            placeholder="/path/to/seqs.fa (headers as taxid.accession)",
        )
        attach_msa = st.checkbox(
            "Attach MSA to leaves (ETE4 viewer)",
            value=False,
            help="Show the aligned sequences next to leaves in the ETE4 interactive explorer.",
        )
        msa_range = st.text_input(
            "MSA display range (optional)",
            placeholder="e.g. 30:60 (alignment columns, 0-based)",
            help="Show only alignment columns start:end in the MSA. Applies to "
            "both the ETE4 viewer and the static image. Leave blank for the full alignment.",
        )

        # --- INDEPENDENT CHECKBOXES ---
        use_ete4 = st.checkbox("Start ETE4 Interactive Server", value=False)
        ete4_port = st.number_input("ETE4 Port", value=5001) if use_ete4 else 5001

        render_static_ete = st.checkbox(
            "Generate Static ETE4 Image (Custom Domains)", value=False
        )
        if render_static_ete:
            st.markdown("**Static image layers:**")
            show_names = st.checkbox("Leaf names", value=True, key="sl_names")
            show_domains = st.checkbox("Domain shapes", value=True, key="sl_domains")
            show_colors = st.checkbox("Branch colouring", value=True, key="sl_colors")
            show_gene = st.checkbox("Gene name labels", value=True, key="sl_gene")
            show_msa = st.checkbox(
                "Multiple sequence alignment", value=False, key="sl_msa"
            )

    if st.button("Run Full Tree Pipeline", type="primary"):

        st.session_state.pop("tree_ready", None)
        st.session_state.pop("viewer_mode", None)
        st.session_state.pop("ete4_port", None)
        st.session_state.pop("output_dir", None)

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
                "--gt",
                trimal_th,
            ]

            if output_dir:
                cmd += ["--output_dir", output_dir.strip()]
                st.session_state["output_dir"] = output_dir.strip()

            # tax
            if (
                "tree_tax_content" in st.session_state
                and st.session_state["tree_tax_content"].strip()
            ):
                tax_str = ",".join(
                    l.strip()
                    for l in st.session_state["tree_tax_content"].splitlines()
                    if l.strip()
                )
            elif tax:
                tax_str = tax.replace(" ", "")
            else:
                tax_str = None
            if tax_str:
                cmd += ["--taxids", tax_str]
                st.caption(
                    f"Taxids loaded: {tax_str[:80]}{'...' if len(tax_str) > 80 else ''}"
                )
            else:
                fa_path = (
                    os.path.join(output_dir.strip(), f"{prefix}.fa")
                    if output_dir
                    else f"{prefix}.fa"
                )
                if os.path.isfile(fa_path):
                    st.info(
                        "No taxids provided — reusing cached sequence file from previous run."
                    )
                else:
                    st.warning(
                        "No taxids detected — fetching ALL taxa. This may be very slow."
                    )

            # exclude_tax
            if (
                "tree_excl_content" in st.session_state
                and st.session_state["tree_excl_content"].strip()
            ):
                excl_str = ",".join(
                    l.strip()
                    for l in st.session_state["tree_excl_content"].splitlines()
                    if l.strip()
                )
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

            cmd += ["--pfam_source", pfam_source]
            cmd += ["--pfam_logic", pfam_logic]
            cmd += ["--color_by", color_by]
            if local_fasta and local_fasta.strip():
                cmd += ["--local_fasta", local_fasta.strip()]
            if msa_range and msa_range.strip():
                cmd += ["--positions", msa_range.strip()]

            # Prevent Qt Crash
            run_env = os.environ.copy()
            run_env["QT_QPA_PLATFORM"] = "offscreen"

            # -------- A: USE ETE4 BACKGROUND SERVER ---------
            if use_ete4:
                cmd += ["--port", str(int(ete4_port))]
                if attach_msa:
                    cmd += ["--MSA"]
                with st.spinner(
                    f"Preparing ETE4 on port {ete4_port}. Alignment & Tree building may take 1-2 minutes..."
                ):
                    os.system(f"fuser -k {ete4_port}/tcp >/dev/null 2>&1")
                    proc = subprocess.Popen(cmd, env=run_env)

                    deadline = time.time() + 600
                    interval = 3
                    connected = False
                    while time.time() < deadline:
                        try:
                            with socket.create_connection(
                                ("localhost", int(ete4_port)), timeout=2
                            ):
                                break
                        except OSError:
                            time.sleep(interval)
                    else:
                        st.error(
                            "ETE4 server did not start within 10 minutes. Check your terminal."
                        )
                        st.stop()

                st.success(f"ETE4 Server launched. Bound to port {ete4_port}.")
                st.info(
                    "If the frame below says 'Cannot connect', the alignment is still running."
                )

                st.session_state["tree_ready"] = True
                st.session_state["tree_prefix"] = prefix
                st.session_state["viewer_mode"] = "ete4"
                st.session_state["ete4_port"] = ete4_port

            # -------- B: USE ETE4 TO RENDER STATIC IMAGE ---------
            elif render_static_ete:

                layers = []
                if st.session_state.get("sl_names", True):
                    layers.append("names")
                if st.session_state.get("sl_domains", True):
                    layers.append("domains")
                if st.session_state.get("sl_colors", True):
                    layers.append("colors")
                if st.session_state.get("sl_gene", True):
                    layers.append("gene")
                if st.session_state.get("sl_msa", False):
                    layers.append("msa")

                cmd += ["--port", str(int(ete4_port))]
                cmd += ["--render_ete_static"]
                cmd += ["--no_explore"]

                if layers:
                    cmd += ["--static_layers", ",".join(layers)]

                with st.spinner(
                    "Generating static ETE image with domains displayed, please wait."
                ):
                    try:
                        subprocess.run(cmd, check=True, env=run_env)
                    except subprocess.CalledProcessError as e:
                        st.error(
                            f"Pipeline crashed (Exit status {e.returncode}). Check your terminal."
                        )
                        st.stop()

                st.session_state["tree_ready"] = True
                st.session_state["tree_prefix"] = prefix
                st.session_state["viewer_mode"] = "ete4_static"

            # -------- C: USE D3 VIEWER (DEFAULT, NO SERVER) --------
            else:
                cmd += ["--no_explore"]
                with st.spinner(
                    "Building alignment and tree... Check your terminal for live progress."
                ):
                    try:
                        subprocess.run(cmd, check=True, env=run_env)
                    except subprocess.CalledProcessError as e:
                        st.error(
                            f"Pipeline crashed (Exit status {e.returncode}). Check your terminal."
                        )
                        st.stop()

                # --- DYNAMIC FILENAME RECONSTRUCTION BASED ON TRIMAL THRESHOLD ---
                aln_ext = ".mft" if aln == "mafft" else f".{aln}"
                trim_ext = f".gt{trimal_th.replace('.', '')}"

                if ml == "fasttree":
                    tree_file_base = f"{prefix}{aln_ext}{trim_ext}.lg.fasttree"
                elif ml == "iqtree":
                    tree_file_base = f"{prefix}{aln_ext}{trim_ext}.treefile"

                tree_file_path = tree_file_base
                if output_dir:
                    tree_file_path = os.path.join(output_dir.strip(), tree_file_base)

                if os.path.isfile(tree_file_path):
                    with open(tree_file_path) as f:
                        st.session_state["tree_data"] = f.read()

                    st.session_state["tree_prefix"] = prefix
                    st.session_state["tree_method"] = ml
                    st.session_state["final_tree_basename"] = (
                        tree_file_base  # Save exact name
                    )

                    itol_colors_path = f"{tree_file_path}.itol_colors.txt"
                    if os.path.isfile(itol_colors_path):
                        with open(itol_colors_path) as f:
                            st.session_state["itol_data"] = f.read()

                    itol_domains_path = f"{tree_file_path}.itol_domains.txt"
                    if os.path.isfile(itol_domains_path):
                        with open(itol_domains_path) as f:
                            st.session_state["itol_domains"] = f.read()

                    st.session_state["tree_ready"] = True
                    st.session_state["viewer_mode"] = "d3"
                else:
                    st.error(
                        f"Pipeline finished but expected file not found:\n`{tree_file_path}`"
                    )

    # ==========================================================
    # Display the correct viewer based on what was just run
    # ==========================================================
    if st.session_state.get("tree_ready"):
        p = st.session_state["tree_prefix"]
        out_dir = st.session_state.get("output_dir", "")
        st.markdown("---")

        # -------- SHOW ETE4 INTERACTIVE SERVER --------
        if st.session_state.get("viewer_mode") == "ete4":
            st.subheader("ETE4 Interactive Explorer")
            port = st.session_state["ete4_port"]
            st.caption(f"Connected to ETE4 server on port {port}")
            st.iframe(
                f"http://localhost:{port}/static/gui.html?tree=tree-1",
                width=1800,
                height=1200,
            )

        # -------- SHOW ETE4 STATIC IMAGE --------
        elif st.session_state.get("viewer_mode") == "ete4_static":
            st.subheader("Phylogeny & Domain Architecture (Static High-Res)")

            expected_img_path = f"{p}_tree_domains.png"
            if out_dir:
                expected_img_path = os.path.join(out_dir, expected_img_path)

            if os.path.isfile(expected_img_path):
                st.image(expected_img_path, width="stretch")
                with open(expected_img_path, "rb") as img_file:
                    st.download_button(
                        label="Download Motif Tree PNG",
                        data=img_file.read(),
                        file_name=os.path.basename(expected_img_path),
                        mime="image/png",
                        key="ete4_png_dl",
                    )
            else:
                st.error(
                    f"Image not found at {expected_img_path}. Check terminal for ETE4 rendering errors."
                )

        # -------- SHOW D3 VIEWER --------
        elif st.session_state.get("viewer_mode") == "d3":
            st.subheader("Interactive Tree Preview")
            st.caption(
                "Scroll to zoom, click nodes to collapse/expand. Hover over leaves for details."
            )

            # Look up the exact dynamic name we calculated earlier
            tree_file_base = st.session_state.get("final_tree_basename")
            if out_dir:
                tree_file_base = os.path.join(out_dir, tree_file_base)

            itol_path = f"{tree_file_base}.itol_colors.txt"

            leaf_colors = None
            if os.path.isfile(itol_path):
                leaf_colors = itc.parse_itol_colors(itol_path)

            tree_html = itc.build_tree_html(
                newick_str=st.session_state["tree_data"],
                leaf_colors=leaf_colors,
                title=f"Phylogeny: {p}",
                height=800,
            )

            components.html(tree_html, height=800, scrolling=True)

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

        # -------------- Newick and iTOL downloads --------------
        if st.session_state.get("viewer_mode") in ["d3", "ete4_static"]:
            st.caption("For publication figures, upload the Newick to iTOL:")
            dl_col1, dl_col2, dl_col3 = st.columns(3)
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
            with dl_col3:
                if "itol_domains" in st.session_state:
                    st.download_button(
                        "Download iTOL Domains",
                        st.session_state["itol_domains"],
                        file_name=f"{p}_itol_domains.txt",
                        key="domains_itol_dl",
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

    if st.button("Find Domain Profiles", type="primary"):
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
        pa_tax_file = st.file_uploader(
            "Or upload TXT (one ID per line)", type=["txt"], key="pa_tax_file"
        )
        if pa_tax_file is not None:
            st.session_state["pa_tax_content"] = pa_tax_file.getvalue().decode()

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
    if st.button("Build Matrix", type="primary", key="pa_build"):
        if not pa_pfam_input.strip():
            st.warning("Please enter at least one Pfam name or accession.")
        else:
            pfam_queries = [q.strip() for q in pa_pfam_input.split(",") if q.strip()]
            if (
                "pa_tax_content" in st.session_state
                and st.session_state["pa_tax_content"].strip()
            ):
                tax_ids = [
                    int(l.strip())
                    for l in st.session_state["pa_tax_content"].splitlines()
                    if l.strip()
                ]
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

            # -------------- Build a human-readable taxon label --------------

            def taxon_label(taxon_id, scientific_name):
                if scientific_name:
                    return f"{taxon_id} · {scientific_name}"
                return str(taxon_id)

            df_flat = pd.DataFrame(rows)

            # Add the label column for display; keep taxon_id for queries.
            df_flat["taxon_label"] = df_flat.apply(
                lambda r: taxon_label(r["taxon_id"], r.get("scientific_name")), axis=1
            )

            # -------------- Pivot into matrix  --------------
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

            # -------------- Color the matrix --------------
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

            # ──--------------- Clustered heatmap --------------
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
                if st.button("Drill Down", type="primary", key="pa_drill"):
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

                        # -------------- Path A: sub-profile hits --------------
                        # We exclude the original query profile from the
                        # results so it doesn't dominate the view
                        subprofile_rows = uni.fetch_subprofile_hits(
                            pa_ver,
                            accessions,
                            pa_evalue,
                            exclude_queries=[selected_profile],
                            db_config=config,
                        )

                        # -------------- Path B: domain architectures --------------
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

        # -------------- Tab A: sub-profile enrichment --------------
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

        # -------------- Tab B: domain architectures --------------
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

                # -------------- Architecture diagram for a selected pattern --------------
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

# ==========================================================================================================

#                                    EXTRACT DOWNLOADED BRANCH TAB

# ==========================================================================================================

elif choice == "Extract Downloaded Branch":
    st.header("Extract Accessions & Sequences from Downloaded Branch")
    st.markdown("""
    1. In the ETE4 Interactive Viewer, right-click the root of the clade you want.
    2. Click **"Download branch as newick"**.
    3. Upload that specific file below to extract its accessions and fetch their sequences.
    """)

    col1, col2 = st.columns(2)
    with col1:
        uploaded_branch = st.file_uploader(
            "Upload Branch (.nwk, .nw, .tree)", type=["nwk", "nw", "txt", "tree"]
        )
    with col2:
        # Version to query the database for the sequences
        ver = st.text_input("UniProt Version", value="2026_01", key="branch_ver")

    if uploaded_branch and st.button("Extract Data", type="primary"):
        import tempfile
        import os
        from ete4 import PhyloTree

        # Streamlit file uploaders return bytes -> save to a temp file for ETE4
        with tempfile.NamedTemporaryFile(delete=False, suffix=".nwk") as tmp:
            tmp.write(uploaded_branch.getvalue())
            tmp_path = tmp.name

        try:
            t = PhyloTree(tmp_path)
            extracted_accs = []

            # Grab every leaf from the branch
            for leaf in t.leaves():
                # Clean up the name ("taxid.accession" format)
                name_parts = leaf.name.split(".")
                acc = name_parts[1] if len(name_parts) > 1 else leaf.name
                extracted_accs.append(acc)

            if extracted_accs:
                st.success(
                    f"Successfully extracted {len(extracted_accs)} accessions from the tree!"
                )

                with st.spinner("Fetching sequences from the database..."):
                    # Query the DB for the full protein records
                    records = uni.fetch_sequences_by_accession(
                        ver, extracted_accs, db_config=config
                    )

                if records:
                    st.success(f"Retrieved {len(records)} sequences from the database!")

                    # Convert records to a FASTA string
                    with uni.UniProtRetriever(config) as db:
                        fasta_str = db.to_fasta_string(records)

                    txt_output = "\n".join(sorted(extracted_accs))

                    # Provide download buttons for both formats
                    dl_col1, dl_col2 = st.columns(2)
                    with dl_col1:
                        st.download_button(
                            label="Download Accessions (.txt)",
                            data=txt_output,
                            file_name="branch_accessions.txt",
                        )
                    with dl_col2:
                        st.download_button(
                            label="Download Sequences (.fasta)",
                            data=fasta_str,
                            file_name="branch_sequences.fasta",
                        )

                    with st.expander("Preview FASTA"):
                        # Just show the first 1500 characters so it doesn't lag the browser
                        st.text(
                            fasta_str[:1500]
                            + ("\n... [truncated]" if len(fasta_str) > 1500 else "")
                        )
                else:
                    st.warning(
                        "Could not find sequences for these accessions in the database. Double-check your UniProt Version."
                    )
                    # Still allow them to download the accessions text file even if the DB lookup fails
                    txt_output = "\n".join(sorted(extracted_accs))
                    st.download_button(
                        label="Download Accessions (.txt)",
                        data=txt_output,
                        file_name="branch_accessions.txt",
                    )

            else:
                st.warning("No leaves found in the uploaded branch.")

        except Exception as e:
            st.error(f"Error parsing branch file: {e}")
        finally:
            # Clean up the temporary file
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

# ==========================================================================================================

#                                    HIGH-RESOLUTION PHYLOGENETIC PROFILE TAB

# ==========================================================================================================

elif choice == "High-Resolution Phylogenetic Profile":
    import subclade_partition as sp
    import tree_builder as tb

    st.header("High-Resolution Phylogenetic Profile")
    st.markdown("""
    Build a profile where columns are **subclades within each Pfam's gene tree**
    (paralog groups), not just Pfam families. Reveals lineage-specific paralog
    expansions and contractions that a family-level matrix would average out.

    **Workflow**
    1. Enter Pfams + taxa → build one tree per Pfam.
    2. For each tree, partition into subclades (depth slider or manual MRCA picks).
    3. Lock in each Pfam's partition, then assemble the profile matrix.
    """)

    # -------------------------------------------------------------
    # SECTION 1: Input form
    # -------------------------------------------------------------
    with st.expander(
        "**1. Input — Pfams, taxa, and tree-build parameters**",
        expanded=not st.session_state.get("hrp_trees_built", False),
    ):
        col1, col2 = st.columns(2)
        with col1:
            version = st.text_input(
                "UniProt Version", value="2026_01", key="hrp_version"
            )
            pfam_text = st.text_area(
                "Pfam IDs or HMM names (comma/newline separated)",
                placeholder="PF00041, PF00069, PF00100",
                key="hrp_pfam_text",
                height=80,
            )
            evalue = st.number_input(
                "E-value cutoff (optional)",
                value=None,
                format="%.1e",
                key="hrp_evalue",
            )
        with col2:
            taxids_text = st.text_input(
                "Taxonomy IDs (comma separated, optional)", key="hrp_taxids_text"
            )
            tax_file = st.file_uploader(
                "Or upload txt (one ID per line)", type=["txt"], key="hrp_tax_file"
            )
            if tax_file is not None:
                st.session_state["hrp_tax_content"] = tax_file.getvalue().decode()

            exclude_tax = st.text_input(
                "Exclude Taxonomy IDs (comma separated, optional)",
                key="hrp_exclude_text",
            )
            exclude_tax_file = st.file_uploader(
                "Or upload exclude txt", type=["txt"], key="hrp_excl_file"
            )
            if exclude_tax_file is not None:
                st.session_state["hrp_excl_content"] = (
                    exclude_tax_file.getvalue().decode()
                )

        tree_output_root = st.text_input(
            "Tree output / cache directory",
            value="/home/user/highres_runs",
            key="hrp_output_root",
            help="Per-config subdirectories are created here. Re-runs with the same parameters reuse cached trees.",
        )

        with st.expander("Advanced tree-build options", expanded=False):
            ca, cb, cc, cd = st.columns(4)
            with ca:
                aln = st.selectbox(
                    "Aligner", ["mafft", "einsi", "clustalo"], index=0, key="hrp_aln"
                )
            with cb:
                ml = st.selectbox(
                    "Tree method", ["fasttree", "iqtree"], index=0, key="hrp_ml"
                )
            with cc:
                gt = st.text_input("trimAl -gt", value="0.01", key="hrp_gt")
            with cd:
                cpu = st.number_input(
                    "CPU threads", min_value=1, value=8, key="hrp_cpu"
                )

        if st.button("Build trees", type="primary", key="hrp_build_btn"):
            # Parse Pfams: split on commas/whitespace/newlines, dedupe
            seen, pfams = set(), []
            for tok in pfam_text.replace(",", " ").split():
                if tok and tok not in seen:
                    seen.add(tok)
                    pfams.append(tok)

            # Combine taxid text input + uploaded file content
            tax_combined = (
                taxids_text + "\n" + st.session_state.get("hrp_tax_content", "")
            )
            taxids = []
            for tok in tax_combined.replace(",", " ").split():
                try:
                    taxids.append(int(tok))
                except ValueError:
                    pass

            # Combine exclude taxid text input + uploaded file content
            excl_combined = (
                exclude_tax + "\n" + st.session_state.get("hrp_excl_content", "")
            )
            exclude_taxids = []
            for tok in excl_combined.replace(",", " ").split():
                try:
                    exclude_taxids.append(int(tok))
                except ValueError:
                    pass

            if not pfams:
                st.error("Please supply at least one Pfam ID or HMM name.")
            else:
                # Stream progress as trees build
                prog_area = st.empty()
                log_lines = []

                def cb(msg):
                    log_lines.append(msg)
                    prog_area.code("\n".join(log_lines[-15:]))

                with st.spinner(f"Building {len(pfams)} tree(s)..."):
                    results = tb.build_trees(
                        pfams=pfams,
                        output_root=tree_output_root,
                        version=version,
                        taxids=taxids or None,
                        exclude_taxids=exclude_taxids or None,
                        evalue=float(evalue) if evalue else None,
                        aln=aln,
                        ml=ml,
                        gt=gt,
                        cpu=int(cpu),
                        tree_from_db_path="tree_from_db.py",
                        progress_callback=cb,
                    )

                st.session_state["hrp_tree_results"] = results
                st.session_state["hrp_trees_built"] = True
                st.session_state["hrp_partitions"] = {}
                st.session_state["hrp_profile_output"] = None
                st.session_state["hrp_build_version"] = version
                st.session_state["hrp_build_taxids"] = taxids

                st.rerun()

    # -------------------------------------------------------------
    # SECTION 2: Tree-build summary
    # -------------------------------------------------------------
    if not st.session_state.get("hrp_trees_built", False):
        st.info("Configure inputs above and click **Build trees** to proceed.")

    else:
        results = st.session_state.get("hrp_tree_results", {})
        partitions = st.session_state.get("hrp_partitions", {})

        st.subheader("2. Tree-build summary")
        summary_rows = []
        for pfam, r in results.items():
            summary_rows.append(
                {
                    "Pfam": pfam,
                    "Status": (
                        "OK" if r["error"] is None else f"ERROR: {r['error'][:60]}"
                    ),
                    "Cached": "yes" if r["cached"] else "no",
                    "Leaves": len(r["leaves"]),
                    "Tree path": r["tree_path"] or "(none)",
                }
            )
        st.dataframe(
            pd.DataFrame(summary_rows), use_container_width=True, hide_index=True
        )

        failed = [p for p, r in results.items() if r["error"]]
        if failed:
            with st.expander(
                f"{len(failed)} Pfam(s) failed — show stderr", expanded=False
            ):
                for p in failed:
                    st.markdown(f"**{p}**")
                    st.code(
                        results[p]["stderr"] or "(no stderr captured)", language="text"
                    )

        # -------------------------------------------------------------
        # SECTION 3: Per-Pfam partition selection
        # -------------------------------------------------------------
        st.subheader("3. Partition each tree into subclades")
        st.markdown(
            "For each Pfam, pick a partition mode. **Depth slider** = cut the "
            "tree at a chosen root-to-node distance. **Manual MRCA** = paste "
            "groups of leaf names (one group per line, comma-separated); the "
            "MRCA of each group becomes a subclade."
        )

        successful_pfams = [p for p, r in results.items() if r["error"] is None]

        for pfam in successful_pfams:
            r = results[pfam]
            tree = r["tree"]
            leaves = r["leaves"]
            n_leaves = len(leaves)
            max_d = sp.get_max_root_distance(tree)

            with st.expander(
                f"**{pfam}** — {n_leaves} leaves, max root depth {max_d:.3f}",
                expanded=(pfam not in partitions),
            ):
                mode = st.radio(
                    "Partition mode",
                    ["Depth slider", "Manual MRCA", "Node path"],
                    horizontal=True,
                    key=f"hrp_mode_{pfam}",
                )

                new_parts = None

                if mode == "Depth slider":
                    default_d = max(0.001, max_d * 0.5)
                    d = st.slider(
                        "Cut depth (root-to-node distance)",
                        min_value=float(0.0),
                        max_value=float(max_d * 1.01),
                        value=float(default_d),
                        step=float(max_d / 200),
                        key=f"hrp_depth_{pfam}",
                    )
                    if d > 0:
                        try:
                            new_parts = sp.partition_by_depth(tree, threshold=d)
                        except Exception as e:
                            st.error(f"Partition failed: {e}")

                elif mode == "Manual MRCA":
                    mrca_text = st.text_area(
                        "Leaf-name groups (one group per line, comma-separated)",
                        placeholder=(
                            "# Each line = one subclade defined by its MRCA.\n"
                            f"{leaves[0] if leaves else ''}, "
                            f"{leaves[1] if len(leaves)>1 else ''}\n"
                            f"{leaves[2] if len(leaves)>2 else ''}"
                        ),
                        height=120,
                        key=f"hrp_mrca_text_{pfam}",
                    )
                    include_unassigned = st.checkbox(
                        "Include 'unassigned' bucket for leaves not under any picked MRCA",
                        value=True,
                        key=f"hrp_unassigned_{pfam}",
                    )

                    # Parse groups: one group per non-blank, non-comment line
                    groups = []
                    for line in mrca_text.splitlines():
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue
                        members = [m.strip() for m in line.split(",") if m.strip()]
                        if members:
                            groups.append(members)

                    if groups:
                        leaf_set = set(leaves)
                        bad = [m for grp in groups for m in grp if m not in leaf_set]
                        if bad:
                            st.warning(
                                f"{len(bad)} leaf name(s) not in this tree — "
                                f"will be ignored. First few: {bad[:5]}"
                            )
                        try:
                            new_parts = sp.partition_by_mrca(
                                tree,
                                groups,
                                include_unassigned=include_unassigned,
                            )
                        except Exception as e:
                            st.error(f"Partition failed: {e}")

                elif mode == "Node path":
                    # Optional node browser
                    if st.checkbox(
                        "Show node list (path + leaf count)", key=f"hrp_nodelist_{pfam}"
                    ):
                        node_rows = sp.list_internal_nodes(tree)
                        st.dataframe(
                            pd.DataFrame(node_rows),
                            use_container_width=True,
                            hide_index=True,
                        )

                    path_text = st.text_area(
                        "Node paths (one per line, comma-separated child indices)",
                        placeholder="# Each line = one subclade, e.g.\n0,1\n1,0\n1,1,0",
                        height=120,
                        key=f"hrp_pathtext_{pfam}",
                    )

                    # Parse: one path per non-blank, non-comment line.
                    paths = []
                    for line in path_text.splitlines():
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue
                        try:
                            paths.append(
                                [int(x) for x in line.split(",") if x.strip() != ""]
                            )
                        except ValueError:
                            st.warning(f"Skipping unparseable line: {line}")

                    if paths:
                        try:
                            new_parts = sp.partition_by_node_path(tree, paths)
                        except Exception as e:
                            st.error(f"Partition failed: {e}")

                # Partition preview
                if new_parts:
                    summary = pd.DataFrame(
                        [
                            {
                                "Subclade": label,
                                "# leaves": len(members),
                                "Sample leaves": ", ".join(sorted(members)[:3])
                                + (" ..." if len(members) > 3 else ""),
                            }
                            for label, members in new_parts.items()
                        ]
                    )
                    st.dataframe(summary, use_container_width=True, hide_index=True)

                    lock_col, status_col = st.columns([1, 4])
                    with lock_col:
                        if st.button("Use this partition", key=f"hrp_lock_{pfam}"):
                            partitions[pfam] = tb.strip_leaf_prefix_in_subclades(
                                new_parts
                            )
                            st.session_state["hrp_partitions"] = partitions
                            st.session_state["hrp_profile_output"] = None
                            st.rerun()
                    with status_col:
                        if partitions.get(pfam):
                            st.success(
                                f"Partition locked in ({len(partitions[pfam])} subclades)."
                            )

                # Tree preview + Newick download
                with st.expander("Show tree preview / download Newick"):
                    try:
                        with open(r["tree_path"]) as fh:
                            nwk_str = fh.read()

                        preview_mode = st.radio(
                            "Display mode",
                            [
                                "Static PNG (Bio.Phylo)",
                                "Interactive D3 (inline)",
                                "ETE4 smartview (server)",
                            ],
                            horizontal=True,
                            key=f"hrp_preview_mode_{pfam}",
                        )

                        # -------------- Mode 1: static PNG, current behavior --------------
                        if preview_mode == "Static PNG (Bio.Phylo)":
                            if n_leaves <= 100:
                                st.image(
                                    viz.render_tree(nwk_str), use_container_width=True
                                )
                            else:
                                st.info(
                                    f"Tree has {n_leaves} leaves — too many for a static PNG. "
                                    "Switch to **Interactive D3** or **ETE4 smartview** above."
                                )

                        # -------------- Mode 2: D3 interactive, inline in Streamlit --------------
                        elif preview_mode == "Interactive D3 (inline)":
                            import interactive_tree_component as itc

                            html = itc.build_tree_html(newick_str=nwk_str, title=pfam)
                            st.components.v1.html(html, height=700, scrolling=True)
                            st.caption(
                                "Zoom/pan with mouse, click internal nodes to collapse, "
                                "hover leaves for details. Note the leaf names you want "
                                "to group and paste them into **Manual MRCA** above."
                            )

                        # -------------- Mode 3: ETE4 smartview, external server --------------
                        elif preview_mode == "ETE4 smartview (server)":

                            pfam_idx = list(results.keys()).index(pfam)
                            default_port = 5001 + pfam_idx

                            server_key = f"hrp_smartview_pid_{pfam}"
                            port_key = f"hrp_smartview_port_{pfam}"
                            running_pid = st.session_state.get(server_key)
                            running_port = st.session_state.get(port_key)

                            col_a, col_b, col_c = st.columns([1, 1, 3])

                            with col_a:
                                port = st.number_input(
                                    "Port",
                                    min_value=1024,
                                    max_value=65535,
                                    value=(
                                        running_port
                                        if running_pid is not None
                                        else default_port
                                    ),
                                    step=1,
                                    key=f"hrp_port_input_{pfam}",
                                    disabled=(running_pid is not None),
                                )

                            with col_b:
                                if running_pid is None:
                                    if st.button(
                                        "Launch",
                                        key=f"hrp_launch_{pfam}",
                                        use_container_width=True,
                                    ):
                                        # Probe the port first for a clean error
                                        sock = socket.socket(
                                            socket.AF_INET, socket.SOCK_STREAM
                                        )
                                        try:
                                            sock.bind(("localhost", int(port)))
                                            port_free = True
                                        except OSError:
                                            port_free = False
                                        finally:
                                            sock.close()

                                        if not port_free:
                                            st.error(
                                                f"Port {port} is already in use. Pick a different "
                                                f"port, or run `ss -tlnp | grep {port}` on the "
                                                "server to find/kill the process holding it."
                                            )
                                        else:
                                            os.system(
                                                f"fuser -k {int(port)}/tcp >/dev/null 2>&1"
                                            )
                                            run_env = os.environ.copy()
                                            run_env["QT_QPA_PLATFORM"] = "offscreen"
                                            # Open the RESOLVED tree so node paths in the viewer
                                            # match partition_by_node_path exactly.
                                            cmd = [
                                                sys.executable,
                                                "tree_from_db.py",
                                                "--pfam",
                                                pfam,
                                                "--version",
                                                st.session_state.get(
                                                    "hrp_build_version", "2026_01"
                                                ),
                                                "--prefix",
                                                r["prefix"],
                                                "--port",
                                                str(int(port)),
                                                # "--no_ncbi",
                                                "--use_resolved",
                                            ]
                                            proc = subprocess.Popen(cmd, env=run_env)

                                            with st.spinner(
                                                f"Starting ETE4 smartview on port {port}..."
                                            ):
                                                deadline = time.time() + 120
                                                connected = False
                                                while time.time() < deadline:
                                                    try:
                                                        with socket.create_connection(
                                                            ("127.0.0.1", int(port)),
                                                            timeout=2,
                                                        ):
                                                            connected = True
                                                            break
                                                    except OSError:
                                                        time.sleep(2)

                                            if not connected:
                                                st.error(
                                                    "ETE4 server did not start within 2 minutes. Check terminal."
                                                )
                                            else:
                                                st.session_state[server_key] = proc.pid
                                                st.session_state[port_key] = int(port)
                                                st.rerun()
                                else:
                                    if st.button(
                                        "Stop",
                                        key=f"hrp_stop_{pfam}",
                                        use_container_width=True,
                                    ):
                                        try:
                                            os.kill(running_pid, signal.SIGKILL)
                                        except ProcessLookupError:
                                            pass
                                        os.system(
                                            f"fuser -k {running_port}/tcp >/dev/null 2>&1"
                                        )
                                        del st.session_state[server_key]
                                        del st.session_state[port_key]
                                        st.rerun()

                            with col_c:
                                if running_pid is not None:
                                    st.success(
                                        f"Smartview running on port {running_port} (PID {running_pid})"
                                    )
                                    st.caption(
                                        f"Inline viewer below. SSH tunnel must include "
                                        f"`-L {running_port}:localhost:{running_port}`. "
                                        "Read a node path in the viewer, paste it into **Node path** above."
                                    )
                                else:
                                    st.caption(
                                        "Opens the resolved tree so viewer node paths match "
                                        "Node-path mode. Each Pfam needs its own port; stop when done."
                                    )

                            if running_pid is not None:
                                components.iframe(
                                    f"http://localhost:{running_port}/static/gui.html?tree=tree-1",
                                    width=1800,
                                    height=900,
                                )

                        # -------------- Newick download, available in all modes --------------
                        st.download_button(
                            "Download Newick",
                            data=nwk_str,
                            file_name=f"{pfam}.nwk",
                            mime="text/plain",
                            key=f"hrp_dl_nwk_{pfam}",
                        )
                    except Exception as e:
                        st.warning(f"Could not preview tree: {e}")

        # -------------------------------------------------------------
        # SECTION 4: Profile assembly + heatmap
        # -------------------------------------------------------------
        st.subheader("4. Assemble profile")

        n_locked = sum(1 for p in successful_pfams if partitions.get(p))
        if n_locked == 0:
            st.info(
                "Lock in at least one Pfam partition above to assemble the profile."
            )
        else:
            st.markdown(
                f"**{n_locked} / {len(successful_pfams)} Pfam(s) have a locked partition.**"
            )

            co1, co2, co3 = st.columns(3)
            with co1:
                binary = st.checkbox(
                    "Binary (0/1) instead of counts", value=False, key="hrp_binary"
                )
            with co2:
                log_scale = st.checkbox(
                    "Log-scale colors (counts only)",
                    value=False,
                    key="hrp_log_scale",
                    disabled=binary,
                )
            with co3:
                cluster_cols = st.checkbox(
                    "Cluster columns",
                    value=False,
                    key="hrp_cluster_cols",
                )

            if st.button(
                "Compute / refresh profile", type="primary", key="hrp_compute_btn"
            ):
                pfam_subclade_map = {
                    p: partitions[p] for p in successful_pfams if partitions.get(p)
                }
                with st.spinner("Querying database and assembling matrix..."):
                    out = uni.fetch_highres_profile(
                        version=st.session_state.get("hrp_build_version", "2026_01"),
                        pfam_subclade_map=pfam_subclade_map,
                        taxon_ids=st.session_state.get("hrp_build_taxids", None)
                        or None,
                        binary=binary,
                        db_config=config,
                    )
                st.session_state["hrp_profile_output"] = out
                st.session_state["hrp_profile_binary"] = binary
                st.session_state["hrp_profile_log"] = log_scale
                st.session_state["hrp_profile_cluster_cols"] = cluster_cols

            profile_output = st.session_state.get("hrp_profile_output", None)

            if profile_output is None:
                st.info(
                    "Click **Compute / refresh profile** to generate the matrix and heatmap."
                )

            elif profile_output["matrix"].empty:
                st.warning("Profile matrix is empty — no proteins found.")

            else:
                matrix = profile_output["matrix"]
                missing = profile_output["missing_accessions"]

                if missing:
                    st.warning(
                        f"{len(missing)} accession(s) from the input trees were not "
                        f"found in the database for this version and were excluded "
                        f"from the counts."
                    )
                    with st.expander("Show missing accessions"):
                        st.code("\n".join(sorted(missing)))

                # Heatmap
                st.markdown("**Heatmap**")
                buf = viz.draw_highres_profile_heatmap(
                    matrix,
                    column_origin=profile_output["column_origin"],
                    taxon_names=profile_output["taxon_names"],
                    missing_accessions=missing,
                    binary=st.session_state.get("hrp_profile_binary", False),
                    log_scale=st.session_state.get("hrp_profile_log", False),
                    cluster_cols=st.session_state.get(
                        "hrp_profile_cluster_cols", False
                    ),
                )
                st.image(buf, use_container_width=True)

                # Tabular matrix
                st.markdown("**Matrix**")
                display_matrix = matrix.copy()
                tn = profile_output["taxon_names"]
                display_matrix.index = [
                    f"{tx}  {tn.get(tx, '')}" for tx in display_matrix.index
                ]
                st.dataframe(display_matrix, use_container_width=True)

                # Downloads + reset
                dl1, dl2, dl3 = st.columns(3)
                with dl1:
                    st.download_button(
                        "Download Matrix (.csv)",
                        data=matrix.to_csv().encode(),
                        file_name="highres_profile.csv",
                        mime="text/csv",
                        key="hrp_dl_csv",
                    )
                with dl2:
                    st.download_button(
                        "Download Heatmap (.png)",
                        data=buf.getvalue(),
                        file_name="highres_profile.png",
                        mime="image/png",
                        key="hrp_dl_png",
                    )
                with dl3:
                    if st.button("Reset entire tab", key="hrp_reset_btn"):
                        for k in list(st.session_state.keys()):
                            if k.startswith("hrp_"):
                                del st.session_state[k]
                        st.rerun()
