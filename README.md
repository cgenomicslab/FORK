# UniProt Lab Manager

A Streamlit-based tool for querying, visualizing, and analyzing a **local UniProt Reference Proteomes database** enriched with Pfam-A HMM search results.

Built at CGLab. The database is built once (see `setup/`) and the tool runs against it locally — no internet access needed for queries.

---

## What it does

The tool provides 8 functional tabs:

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

The **Presence/Absence** tab implements a two-tier workflow:
- **Step 1**: Detect presence/absence of a broad Pfam profile across organisms (e.g. all decarboxylases in Euk/Bact/Arch)
- **Step 2**: Drill into a cell → sub-profile enrichment (TIGRFAM, specific Pfams) or domain architecture co-occurrence patterns → identifies specific enzymes (e.g. Histamine DC vs Tyramine DC vs Dopamine DC)

---

## Repository structure

```
uniprot-lab-manager/
│
├── README.md
├── requirements.txt
├── .gitignore
│
├── retrieval_script_gui.py        # Main entry point — run this with Streamlit
├── get_reference_uniprot_set_lib.py  # Backend retrieval library
├── tree_from_db.py                # Alignment + tree pipeline (MAFFT/FastTree/IQ-TREE)
├── viz_utils.py                   # Visualization: domain architecture, static tree, heatmap
├── interactive_tree_component.py  # D3-based interactive tree (HTML/JS)
│
└── setup/                         # One-time setup scripts (run by DB admin only)
    ├── README.md
    ├── uniprot_sync_v7.py         # Builds the local UniProt reference proteomes DB
    └── pyhmmer_hmmsearch.py      # Runs Pfam-A HMM search across the DB (PyHMMER)
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
python -m streamlit run retrieval_script_gui.py
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

### `tree_from_db.py`
Called by the GUI as a subprocess. Fetches sequences from the DB, aligns them, and builds a phylogenetic tree.

```bash
python tree_from_db.py \
  --pfam PF00046,Homeodomain \
  --version 2026_01 \
  --taxids 9606,10090,7227 \
  --prefix myrun \
  --aln mafft \
  --ml fasttree \
  --cpu 32
```

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

Outputs: `.fa`, `.mft`, `.mft.gt01`, `.mft.gt01.lg.fasttree`, `.itol_colors.txt`

Use `--no_explore` to skip the ETE4 server (D3 viewer mode in GUI).

### `viz_utils.py`
Standalone visualization helpers. Returns `BytesIO` PNG buffers, compatible with `st.image()`.

- `draw_domain_architecture(domain_records)` — horizontal bar diagram, one row per protein, colored domain blocks at alignment positions. Color is deterministic by domain name (MD5 hash → HSV).
- `render_tree(newick_string)` — static Bio.Phylo + matplotlib tree rendering, auto-scales to leaf count.
- `draw_presence_absence_heatmap(matrix_df)` — seaborn clustermap with hierarchical clustering of both rows and columns.

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

- Output files (`.fa`, `.mft`, `.nwk`, etc.) are generated in the working directory and excluded from version control via `.gitignore`.
- The `.env` file contains credentials and must never be committed.
- For publication-quality trees, download the `.nwk` file and upload to [iTOL](https://itol.embl.de) with the `.itol_colors.txt` color strip file.
- ETE4 interactive explorer requires port access (default 5001). The D3 viewer (default mode) has no server dependency.
