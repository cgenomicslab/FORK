# UniProt Lab Manager

A Flask web app for querying, visualizing, and comparing proteins across a **local UniProt Reference Proteomes database** enriched with Pfam-A HMM search results. Built at CGLab — no internet needed for queries once the database is set up.

![home tab](figures/home_tab.png)

---

## What it does

Three main analysis modules, plus sequence/domain utilities:

| Module | What you get |
|--------|-------------|
| **Phylogenetic Tree** | Fetch sequences by Pfam/taxon → align (MAFFT) → tree (FastTree/IQ-TREE) → interactive D3 viewer or ETE4 explorer with domain shapes |
| **Presence / Absence** | Taxa × Pfam profile heatmap; drill into any cell for sub-profiles or domain architecture breakdown |
| **High-Res Profile** | Partition gene trees into subclade (paralog) groups and profile each subclade separately across taxa |
| **Utilities** | Standard retrieval, HMM search, accession lookup, domain coordinates, GO→domain profiles, branch extraction |

---

## Quickstart

### 1. Build the database

This tool requires a local MySQL database. Build it first with the companion repo:
→ https://github.com/athenamarou/Ref_Proteomes_Local_DB

Run in order:
1. `uniprot_sync_v7.py` — builds the database
2. `pyhmmer_hmmsearch.py` — populates HMM search results

### 2. Install dependencies

```bash
conda env create -f uniprot-lab-manager.yml
conda activate bio_tools
```

### 3. Configure database connection

Create a `.env` file:

```
DB_HOST=localhost
DB_USER=your_db_user
DB_PASSWORD=your_password
DB_NAME=uniprot_db_cglab
```

Or fill in the **DB Config** panel in the Streamlit sidebar at runtime.

### 4. Run

```bash
bash run_webapp.sh
# or directly:
python app.py
```

Open `http://localhost:5000`.

---

## Toy example

### GUI — Presence/Absence workflow

Query a few metabolic Pfam profiles across 5 model organisms, get a clustered heatmap, then drill into any cell to see domain architectures:

**Step 1** — enter taxon IDs (`9606, 10090, 9598, 7227, 6087`) and Pfam names (`Glycolytic, His_Phos_1, COX2, PFK`) in the Presence/Absence tab → clustered heatmap:

![presence absence heatmap](figures/presence_absence.png)

**Step 2** — click any cell (e.g. *Glycolytic* × Mouse) → domain architecture breakdown:

![domain architecture drill-down](figures/domain_arch_presence_abs.png)

---

### CLI — build a phylogenetic tree

```bash
python tree_from_db.py \
  --pfam Glycolytic \
  --version 2026_01 \
  --prefix /tmp/glycolytic_tree \
  --taxids 9606,10090,7227 \
  --aln mafft \
  --ml fasttree \
  --no_ncbi \
  --no_explore
```

Outputs: `.fa`, `.mft`, `.mft.gt01`, `.nwk`, `.itol_colors.txt`, `.itol_domains.txt`

To open the interactive ETE4 viewer (with domain shapes and node popups):

```bash
python tree_from_db.py \
  --pfam Glycolytic \
  --version 2026_01 \
  --prefix /tmp/glycolytic_tree \
  --taxids 9606,10090,7227 \
  --aln mafft \
  --ml fasttree
```

![phylogenetic tree ETE viewer](figures/phylogenetic_tree_ete.png)

---

### Python API — fetch sequences programmatically

```python
from get_reference_uniprot_set_lib import fetch_sequences, fetch_sequences_by_hmm_hit

# Fetch all human + mouse proteins
records = fetch_sequences("2026_01", taxon_ids=[9606, 10090])

# Fetch proteins with a Glycolytic domain hit
records = fetch_sequences_by_hmm_hit("2026_01", "Glycolytic", taxon_ids=[9606, 10090, 7227])

for r in records:
    print(r.id, len(r.seq))
```

Returns BioPython `SeqRecord` objects, ready for downstream analysis or writing to FASTA.

---

### High-resolution profile

Split the Glycolytic gene tree into subclade groups (paralogs A/B/C…) and profile each separately:

![high-res profile heatmap](figures/high_res_heatmap.png)

Columns are `Pfam·Subclade` pairs; the color stripe groups subclades by parent Pfam. Export as CSV or PNG.

Species tree with Phylogenetic Profile each taxon for the subclades of interest.

![high-re profile ete tree](figures/phylogenetic_tree_ete.png)


---

## Repository structure

```
uniprot-lab-manager-copy/
├── app.py                            # Entry point — Flask app
├── run_webapp.sh                     # Convenience launcher
├── templates/                        # HTML templates (Jinja2)
│   ├── base.html
│   ├── index.html
│   ├── tree.html
│   ├── presence.html
│   ├── highres.html
│   └── utilities.html
├── static/                           # CSS / JS assets
├── get_reference_uniprot_set_lib.py  # Backend retrieval library (importable)
├── tree_from_db.py                   # CLI: fetch → align → tree → viewer
├── subclade_partition.py             # Partition gene trees into subclades
├── tree_builder.py                   # Per-Pfam tree orchestration + caching
├── interactive_tree_component.py     # D3-based tree viewer (HTML/JS)
├── viz_utils.py                      # Heatmap, domain diagram, tree rendering
├── utils.py                          # Shared helpers
├── figures/                          # Screenshots for README
└── setup/                            # One-time DB build scripts (admin only)
    ├── uniprot_sync_v7.py
    └── pyhmmer_hmmsearch.py
```

---

## Requirements

```
flask, pandas, matplotlib, seaborn, biopython,
mysql-connector-python, python-dotenv, ete4>=4.4.0, numpy
```

All managed via `uniprot-lab-manager.yml`.

---

## Notes

- Output files (`.fa`, `.mft`, `.nwk`, `.itol_*.txt`) are written to the path given by `--prefix` and excluded from version control via `.gitignore`.
- The `.env` file contains credentials — never commit it.
- As an alternative for visualizing trees, load the `.nwk` + `.itol_colors.txt` + `.itol_domains.txt` files into [iTOL](https://itol.embl.de).
- ETE4 interactive viewer requires port 5001 (default). The D3 viewer has no server dependency.
- High-res profiling caches trees by build parameters; rerunning with the same settings reuses the cache.

**[Full API reference & CLI guide](API_REFERENCE.md)**
