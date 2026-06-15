# UniProt Lab Manager

A Streamlit-based tool for querying, visualizing, and analyzing a **local UniProt Reference Proteomes database** enriched with Pfam-A HMM search results.

Built at CGLab. The database is built once (see `setup/`) and the tool runs against it locally — no internet access needed for queries.

---

## What it does

The tool provides 10 functional tabs:

| Tab | Description |
|-----|-------------|
| **Standard Retrieval** | Fetch protein sequences with filters (taxon, GO term, Pfam, proteome ID) |
| **HMM Search** | Retrieve sequences by HMM profile name or accession |
| **Accession Lookup** | Batch fetch sequences by UniProt accession |
| **Domain Coordinate Lookup** | Fetch Pfam domain hits per protein + draw domain architecture diagram |
| **Database Info** | List available UniProt versions and stats |
| **Phylogenetic Tree** | Input HMM names/Pfam IDs → alignment → tree → interactive D3 viewer or ETE4 explorer |
| **GO → Domain Profiles** | For a GO term, return all HMM profiles found in annotated proteins |
| **Presence/Absence & Drill-down** | Cross-table of taxa × Pfam profiles with counts, clustered heatmap, and two-step drill-down into enzyme-level resolution |
| **Extract Downloaded Branch** | After selecting a branch of interest and downloading it as nwk , the user can load the branch file in this tab and extract the accessions and the corresponding protein sequences from the local DB|
| **High-Resolution Phylogenetic Profile** | Build a per-Pfam gene tree, cut it into subclades (paralog groups), and profile presence/counts of each subclade across taxa as a clustered heatmap |


The **Presence/Absence** tab implements a two-tier workflow:
- **Step 1**: Detect presence/absence of a broad Pfam profile across organisms (e.g. all decarboxylases in Euk/Bact/Arch)
- **Step 2**: Drill into a cell → sub-profile enrichment (TIGRFAM, specific Pfams) or domain architecture co-occurrence patterns → identifies specific enzymes (e.g. Histamine DC vs Tyramine DC vs Dopamine DC)

About the **Extract Downloaded Branch** tab: 
- **Step 1**: In the ETE4 Interactive Viewer, right-click the root of the clade you want.
- **Step 2**: Click 'Download branch as newick'.
- **Step 3**: Upload the file in the tab to extract all accessions in selected branch, and use for downstream analysis.

About the **High-Resolution Phylogenetic Profile** tab:

Where Presence/Absence treats a whole Pfam as one column, this tab resolves a Pfam *into its paralog groups* and profiles each one separately. The workflow:
- **Step 1**: For each Pfam, build (or reuse a cached) gene tree via `tree_builder.py`.
- **Step 2**: Partition each tree into subclades — either with a depth slider, by picking MRCAs, or by explicit node paths (`subclade_partition.py`). Subclades are labelled A, B, C, … in ladderized left-to-right order.
- **Step 3**: Assemble the profile matrix (taxa × `Pfam-subclade`) with `fetch_highres_profile()` and render it as a clustered heatmap, with a colored stripe grouping subclades by their parent Pfam. Export as CSV or PNG.


---

## Repository structure

```

uniprot-lab-manager
    ├── .env.example
    ├── .gitignore
    ├── API_REFERENCE.md
    ├── get_reference_uniprot_set_lib.py  # Backend retrieval library
    ├── homepage.py                       # Main entry point — run this with Streamlit
    ├── interactive_tree_component.py     # D3-based interactive tree (HTML/JS)
    ├── logo.png
    ├── pages
    │   └── retrieval_script_gui.py       # Main GUI page (all 10 tabs) — surfaced in the sidebar by Streamlit
    ├── README.md
    ├── setup                             # One-time setup scripts (run by DB admin only)
    │   ├── pyhmmer_hmmsearch.py          # Runs Pfam-A HMM search across the DB (PyHMMER)
    │   ├── README.md
    │   └── uniprot_sync_v7.py            # Builds the local UniProt reference proteomes DB
    ├── subclade_partition.py             # Cut a gene tree into subclades (depth / MRCA / node-path)
    ├── tree_builder.py                   # Per-Pfam tree orchestrator + caching for high-res profiling
    ├── tree_from_db.py                   # Alignment + tree pipeline (MAFFT/FastTree/IQ-TREE)
    ├── uniprot-lab-manager.yml       
    └── viz_utils.py

```

---

## Quickstart

### 1. Prerequisites — Database Setup

This tool requires a local UniProt Reference Proteomes MySQL database.

Build it first using the companion repository:
→ https://github.com/athenamarou/Ref_Proteomes_Local_DB

Run these two steps in order:
1. `uniprot_sync_v7.py` — builds the database
2. `pyhmmer_hmmsearch.py` — populates HMM search results

Once the database is running, proceed with the steps below to launch the GUI.

### 2. Install dependencies

We recommend using Conda to ensure all biological binaries (like MAFFT and FastTree) are installed correctly alongside Python.

```bash
conda env create -f uniprot-lab-manager.yml
conda activate bio_tools
```

### 3. Configure database connection

Create a `.env` file in the scripts directory:

```
DB_HOST=localhost
DB_USER=your_db_user
DB_PASSWORD=your_password
DB_NAME=uniprot_db_cglab
```

Alternatively, fill in the **DB Config** panel in the Streamlit sidebar at runtime.

### 4. Run the app

```bash
python -m streamlit run homepage.py
```

Then open `http://localhost:8501` (or the network URL shown in the terminal).

---

## Script overview

### `get_reference_uniprot_set_lib.py`
Core retrieval library. Can be used standalone or imported in other scripts.

```python
from get_reference_uniprot_set_lib import fetch_sequences, fetch_sequences_by_hmm_hit

# Fetch all human + mouse proteins
records = fetch_sequences("2026_01", taxon_ids=[9606, 10090])

# Fetch proteins with a Homeodomain hit
records = fetch_sequences_by_hmm_hit("2026_01", "Homeodomain", taxon_ids=[9606])
```

Key methods inside `UniProtRetriever`:
- `get_proteins()` — filtered sequence retrieval
- `get_proteins_by_hmm_hit()` — HMM-based retrieval
- `get_domains_by_accession()` — domain coordinates per protein
- `get_presence_absence_matrix()` — Step 1 of the presence/absence workflow
- `get_accessions_for_cell()` — bridge between Step 1 and Step 2
- `get_subprofile_hits()` — Step 2, Path A: deeper HMM resolution
- `get_domain_architectures()` — Step 2, Path B: domain co-occurrence patterns
- `get_accessions_with_taxids()` — bulk accession → taxon_id + scientific name lookup (chunked); maps tree leaves back to the taxa that own them
- `get_highres_profile()` — assemble the high-resolution profile matrix (taxa × `Pfam-subclade`) from a `{pfam: {subclade_label: accessions}}` map

### `tree_from_db.py`
Called by the GUI as a subprocess. Fetches sequences from the DB, aligns them, and builds a phylogenetic tree.

| Argument | Description | Default |
|---|---|---|
| `--pfam` | Pfam ID(s) or HMM name(s), comma-separated | required |
| `--version` | UniProt release (must match DB) | required |
| `--prefix` | Full path + run name for all output files | required |
| `--taxids` | NCBI taxon IDs to include, comma-separated or path to txt file | all taxa |
| `--exclude_taxids` | NCBI taxon IDs to exclude, comma-separated or txt file | none |
| `--aln` | Alignment tool: `mafft`, `einsi`, `clustalo` | `mafft` |
| `--ml` | Tree method: `fasttree`, `iqtree` | `fasttree` |
| `--gt` | TrimAl gap threshold | `0.01` |
| `--cpu` | Threads for alignment and tree building | `4` |
| `--evalue` | E-value cutoff for HMM hits | none |
| `--max_per_taxon` | Keep at most N sequences per taxon before alignment (keeps large gene families manageable) | none |
| `--colormap` | Path to a `taxid<TAB>color` file for the taxon colour strip (otherwise auto-assigned) | none |
| `--no_ncbi` | Skip NCBI taxonomy annotation (faster) | off |
| `--no_explore` | Build tree files only, no viewer | off |
| `--render_ete_static` | Render static PNG with domain shapes using ETE4 | off |
| `--static_layers` | Layers to include in static PNG: `names,domains,colors,gene,msa` | all |
| `--MSA` | Attach aligned sequences to leaves in the ETE viewer (SeqMotifFace) | off |
| `--use_resolved` | Load the pre-resolved `.resolved.nwk` instead of re-rooting/resolving on the fly, so ETE explorer node numbering matches the subclade partitioner | off |
| `--port` | Port for ETE4 interactive server | `5001` |
| `--output_dir` | Directory for all output files (alternative to full path in `--prefix`) | none |

One liner examples to call from CLI:
1. Build tree only (no viewer, fastest)
```bash
python tree_from_db.py \
  --pfam PF00041 \
  --version 2026_01 \
  --prefix /home/user/results/myrun \
  --taxids 9606,10090,7227,7955,6239 \
  --aln mafft \
  --ml fasttree \
  --gt 0.01 \
  --cpu 8 \
  --evalue 1e-5 \
  --no_ncbi \
  --no_explore
```

2. Static image with domain architectures
```bash
python tree_from_db.py \
  --pfam PF00041 \
  --version 2026_01 \
  --prefix /home/user/results/myrun \
  --taxids 9606,10090,7227,7955,6239 \
  --aln mafft \
  --ml fasttree \
  --gt 0.01 \
  --cpu 8 \
  --evalue 1e-5 \
  --no_ncbi \
  --render_ete_static \
  --static_layers names,domains,colors,gene
```

3. Live interactive browser (ETE4) (if run on the server needs ssh -L 5001:localhost:5001 user@server)

```bash
python tree_from_db.py \
  --pfam PF00041 \
  --version 2026_01 \
  --prefix /home/user/results/myrun \
  --taxids 9606,10090,7227,7955,6239 \
  --aln mafft \
  --ml fasttree \
  --gt 0.01 \
  --cpu 8 \
  --evalue 1e-5 \
  --no_ncbi \
  --port 5001
```

4. Reuse cached sequences (skip DB query, change aligner or tree method)

```bash
python tree_from_db.py \
  --pfam PF00041 \
  --version 2026_01 \
  --prefix /home/user/results/myrun \
  --aln einsi \
  --ml iqtree \
  --gt 0.01 \
  --cpu 8 \
  --no_ncbi \
  --no_explore
```

5. Exclude taxa

```bash
python tree_from_db.py \
  --pfam PF00041 \
  --version 2026_01 \
  --prefix /home/user/results/myrun \
  --exclude_taxids 9615,9913 \
  --aln mafft \
  --ml fasttree \
  --gt 0.01 \
  --cpu 8 \
  --evalue 1e-5 \
  --no_ncbi \
  --no_explore
```

6. Cap sequences per taxon (large gene families)

```bash
python tree_from_db.py \
  --pfam PF00046 \
  --version 2026_01 \
  --prefix /home/user/results/homeodomain \
  --taxids 9606,10090,7227,7955,6239 \
  --max_per_taxon 5 \
  --aln mafft \
  --ml fasttree \
  --gt 0.01 \
  --cpu 8 \
  --evalue 1e-5 \
  --no_ncbi \
  --no_explore
```

Outputs: `.fa`, `.mft`, `.mft.gt01`, `.mft.gt01.lg.fasttree`, `.itol_colors.txt`, `.itol_domains.txt`

Use `--no_explore` to skip the ETE4 server (D3 viewer mode in GUI).

### `subclade_partition.py`
Partitions a gene tree (ETE4 `Tree`) into subclades for high-resolution phylogenetic profiling. All modes return the same structure: `{subclade_label: set(leaf_names)}`, where labels run A, B, C, … in ladderized left-to-right order so the column order in the final matrix is reproducible.

Three partitioning modes:
- `partition_by_depth(tree, threshold)` — **Mode 1, depth slider**: cut the tree at a chosen root-to-node distance. Every node that crosses the threshold becomes a subclade root.
- `partition_by_mrca(tree, mrca_specs, include_unassigned=True)` — **Mode 2, manual MRCA picking**: each group of leaf names defines one subclade via its MRCA. Overlaps are first-come-first-served; leftover leaves go to `unassigned`.
- `partition_by_node_path(tree, paths)` — **explicit node paths**: pick only the nodes you name, as child-index paths from the root (`[]` = root, `[1, 1]` = root's second child's second child).

Helpers: `get_max_root_distance(tree)` (upper bound for the depth slider) and `list_internal_nodes(tree)` (internal nodes with their ETE4 path, leaf count, and a sample of leaves — useful for finding paths to pick).

```python
from ete4 import Tree
import subclade_partition as sp

tree = Tree(open("PF00041.nwk").read())
# Mode 1
parts = sp.partition_by_depth(tree, threshold=0.8)
# Mode 2
parts = sp.partition_by_mrca(tree, [
    ["P12345_HUMAN", "Q67890_MOUSE"],
    ["O11111_DROME"],
])
```

Run `python subclade_partition.py` to see all modes demonstrated on a toy tree.

### `tree_builder.py`
Orchestrates per-Pfam gene tree building for the high-resolution profiling pipeline. Calls `tree_from_db.py` as a subprocess (with `--no_ncbi --no_explore`), caches results, and parses the resulting Newick into ETE4 `Tree` objects.

The pipeline:
```
build_trees(pfams, ...)
        │  {pfam: {"tree": ete4.Tree, "leaves": [...], ...}}
        ▼
subclade_partition.partition_by_depth/by_mrca   (per tree)
        │  {pfam: {"A": {accs}, "B": {accs}, ...}}
        ▼
library.fetch_highres_profile(...)
        │
        ▼  pandas matrix
```

**Caching**: each `(pfam, taxids, exclude_taxids, version, evalue, aln, ml, gt)` combination gets its own subdirectory under `output_root`, named after a short MD5 hash of the parameters. Re-running with identical parameters reuses the existing tree; different parameters build in a new subdirectory (old runs are never deleted).

Key functions:
- `build_one_tree(...)` / `build_trees(pfams, ...)` — build (or load from cache) one tree or many.
- `parse_leaf_to_accession()` / `parse_leaf_to_taxid()` / `strip_leaf_prefix_in_subclades()` — convert tree leaves (written as `"{taxid}.{accession}"`) back to bare accessions for the profile-assembly step.
- `cache_key(...)` — deterministic short ID for one tree-build configuration.

FastTree's unrooted, trifurcating root is resolved (`resolve_polytomy` + `ladderize`) and written out as `*.resolved.nwk` so the node numbering shown in the ETE explorer matches what `partition_by_node_path` resolves against.

### `viz_utils.py`
Standalone visualization helpers. Returns `BytesIO` PNG buffers, compatible with `st.image()`.

- `draw_domain_architecture(domain_records)` — horizontal bar diagram, one row per protein, colored domain blocks at alignment positions. Color is deterministic by domain name (MD5 hash → HSV).
- `render_tree(newick_string)` — static Bio.Phylo + matplotlib tree rendering, auto-scales to leaf count.
- `draw_presence_absence_heatmap(matrix_df)` — seaborn clustermap with hierarchical clustering of both rows and columns.
- `draw_highres_profile_heatmap(matrix_df, ...)` — clustered heatmap of the high-resolution profile (rows = taxa, columns = `Pfam-subclade` pairs), with a colored stripe above the columns grouping subclades by parent Pfam. Supports `binary`, `log_scale`, and row/column clustering toggles.

### `interactive_tree_component.py`
Builds a self-contained HTML string with a D3 v7 phylogenetic tree. Pass the result to `st.components.v1.html()`.

Features: zoom/pan, clade collapse/expand, hover tooltips, taxonomy color strip, optional presence/absence heatmap columns, SVG export, ladderize, align-leaves toggle.

```python
import interactive_tree_component as itc

html = itc.build_tree_html(
    newick_str=open("myrun.nwk").read(),
    leaf_colors=itc.parse_itol_colors("myrun.itol_colors.txt"),
    title="Homeodomain family",
)
st.components.v1.html(html, height=800, scrolling=True)
```

---

## Database schema

The local MySQL database (`uniprot_db_cglab`) contains 9 tables:

| Table | Description |
|-------|-------------|
| `proteins` | Core protein table: accession, taxon, proteome, gene/protein name, length |
| `sequences` | Deduplicated sequences (MD5 hash deduplication) |
| `proteomes` | Proteome metadata: proteome_id, taxon_id, description |
| `hmm_search_results` | PyHMMER Pfam-A search results: all domain hits with coordinates and scores |
| `pfam_domains` | Pfam domain metadata: name, description, clan |
| `protein_pfam` | Protein ↔ Pfam domain assignments |
| `go_terms` | GO term definitions |
| `protein_go` | Protein ↔ GO term annotations |
| `sequence_changes` | Tracks accessions that changed sequence between versions |

`hmm_search_results` is the central table for all HMM-based queries. It stores full and domain-level e-values, scores, alignment coordinates (`ali_from`, `ali_to`), HMM model coordinates, and envelope coordinates per domain occurrence.

---

## Requirements

```
streamlit
pandas
matplotlib
seaborn
biopython
mysql-connector-python
python-dotenv
ete4>=4.4.0
numpy
```

---

## For Developers: Python API & CLI

Want to use our retrieval backend in your own Python scripts or Jupyter notebooks? You can import `get_reference_uniprot_set_lib.py` as a standalone library to fetch BioPython `SeqRecord` objects or generate FASTA strings programmatically.

**[Read the API Reference & CLI Guide here](API_REFERENCE.md)**

## Notes

- Output files (`.fa`, `.mft`, `.nwk`, `_tree_domains.png`, `.itol_colors.txt`, `.itol_domains.txt`) are generated at the path specified by `--prefix` and excluded from version control via `.gitignore`.
- The `.env` file contains credentials and must never be committed.
- For better visualizing trees, download the `.nwk` file and upload to iTOL (https://itol.embl.de) with the `.itol_colors.txt` taxon colour strip and `.itol_domains.txt` domain annotation files.
- ETE4 interactive explorer requires port access (default 5001). The D3 viewer (default mode) has no server dependency.
- High-resolution profiling reuses cached trees keyed by build parameters (see `tree_builder.py`); changing taxa, e-value, aligner, or tree method triggers a fresh build in a new subdirectory rather than overwriting.
