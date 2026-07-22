# FORK

A Flask web app for querying, visualizing, and comparing proteins across a **local UniProt Reference Proteomes database** enriched with Pfam-A HMM search results. Built at CGLab — no internet needed for queries once the database is set up.

![home tab](figures/home_tab.png)

---

## What it does

Four main analysis modules, plus sequence/domain utilities:

| Module | What you get |
|--------|-------------|
| **Phylogenetic Tree** | Fetch sequences by Pfam/taxon → align (MAFFT) → tree (FastTree/IQ-TREE) → interactive ETE4 explorer with domain shapes; download the tree as plain Newick |
| **Presence / Absence** | Taxa × Pfam profile heatmap; drill into any cell for sub-profiles or domain architecture breakdown |
| **High-Res Profile** | Partition gene trees into subclade (paralog) groups — by depth, manual MRCA, node path, or automatic duplication — and profile each subclade separately across taxa; optionally combine the DB Pfams with an uploaded FASTA |
| **Comparative Analysis** | "Concept check": find protein families present in one taxon group but absent in another (e.g. Mucorales vs Human) — tables (present-in-A-not-B, the reverse, shared), a heatmap, a group-coloured species tree, and protein-accession drill-down |
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
conda env create -f FORK.yml
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

Or fill in the **Database** panel (top bar) at runtime.

### 4. Run

```bash
bash run_webapp.sh
# or directly:
python app.py
```

Open `http://localhost:8080` (if 8080 is taken the app picks the next free port and prints the URL).

---

## Deployment (production)

FORK runs behind **Gunicorn + Nginx + systemd** with Let's Encrypt TLS. The full step-by-step
runbook is in **[HOSTING_in_laptop.md](HOSTING_in_laptop.md)**. Two things that are easy to get
wrong:

- **Exactly one Gunicorn worker, threaded.** Job state lives in memory, so multiple worker
  processes would break job polling and the ETE4 viewers. Concurrency comes from threads:
  ```bash
  gunicorn --worker-class gthread --workers 1 --threads 8 --timeout 300 --bind 127.0.0.1:8080 app:app
  ```
- **`taxa.sqlite` must exist before the server starts** — the app instantiates `NCBITaxa()` at
  import time (see Notes).

### Users & privacy

The app supports **optional accounts**. Without logging in, each visitor's **Recent runs** are
private to their browser session (cleared when the browser closes); registering an account keeps
them saved and private across sessions and devices. Runs and their results are isolated per user.
Passwords are hashed with werkzeug; accounts live in a local `.users.json` (never committed).

---

## Toy example

The figures below walk one worked example: **neuronal & cell-signalling gene families across animal
evolution** — a panel of ~18 organisms (chordates, arthropods, other invertebrates, a
choanoflagellate, and a plant outgroup) profiled for the Pfam domains of ligand-gated ion channels,
GPCRs, gap junctions, synaptic proteins, and developmental transcription factors.

The organism panel, as a taxon-coloured NCBI species tree:

![species tree of the toy-example organisms](figures/ncbi_tree.png)

### GUI — Presence/Absence workflow

Query several neuronal/signalling Pfam profiles across the panel, get a clustered heatmap, then drill into any cell to see domain architectures:

**Step 1** — enter the taxon IDs
(`9606, 7739, 7668, 7227, 7159, 7165, 7091, 7070, 7460, 6669, 126957, 6945, 6239, 45351, 10228, 400682, 946362, 3702`)
and Pfam names
(`Homeodomain, HLH, Neur_chan_LBD, Neur_chan_memb, Lig_chan, Lig_chan-Glu_bd, 7tm_6, 7tm_7, Innexin, Synaptobrevin, C2`)
in the Presence/Absence tab → clustered heatmap of protein counts per taxon × profile:

![presence absence heatmap](figures/presence_absence.png)

**Step 2** — click any cell (e.g. *Neur_chan_memb* × *Anopheles gambiae*) → domain-architecture breakdown of the matching proteins (here the two-domain **Neur_chan_LBD + Neur_chan_memb** ligand-gated ion-channel architecture):

![domain architecture drill-down](figures/domain_arch_presence_abs.png)

---

### CLI — build a phylogenetic tree

Build a gene tree for the ligand-gated ion-channel family (`Neur_chan_memb`) across the panel:

```bash
python tree_from_db.py \
  --pfam Neur_chan_memb \
  --version 2026_01 \
  --prefix /tmp/neur_chan_tree \
  --taxids 9606,7227,6239,45351,7668,7739,7165,7091,7070 \
  --aln mafft \
  --ml fasttree \
  --no_ncbi \
  --no_explore
```

Outputs: `.fa`, `.mft`, `.mft.gt01`, `.nwk`, `.itol_colors.txt`, `.itol_domains.txt`

To open the interactive ETE4 viewer (branches coloured by taxon, bold species names, and domain shapes aligned on the right):

```bash
python tree_from_db.py \
  --pfam Neur_chan_memb \
  --version 2026_01 \
  --prefix /tmp/neur_chan_tree \
  --taxids 9606,7227,6239,45351,7668,7739,7165,7091,7070 \
  --aln mafft \
  --ml fasttree
```

![phylogenetic tree ETE viewer](figures/phylogenetic_tree_ete.png)

---

### Python API — fetch sequences programmatically

```python
from get_reference_uniprot_set_lib import fetch_sequences, fetch_sequences_by_hmm_hit

# Fetch all human + fruit-fly proteins
records = fetch_sequences("2026_01", taxon_ids=[9606, 7227])

# Fetch proteins with a Neur_chan_memb (ligand-gated ion channel) domain hit
records = fetch_sequences_by_hmm_hit("2026_01", "Neur_chan_memb", taxon_ids=[9606, 7227, 6239])

for r in records:
    print(r.id, len(r.seq))
```

Returns BioPython `SeqRecord` objects, ready for downstream analysis or writing to FASTA.

---

### High-resolution profile

Build gene trees for several neuronal/signalling Pfams at once, split each into subclade (paralog) groups, and profile every subclade separately across the panel:

![high-res profile heatmap](figures/high_res_heatmap.png)

Columns are `Pfam·Subclade` pairs (e.g. `PF02931-A` = Neur_chan_LBD subclade A, `PF00957-A/B` = Synaptobrevin subclades); the colour stripe groups subclades by parent Pfam. Export as CSV or PNG.

You can also **combine the DB Pfams with an uploaded FASTA** (headers `{taxid}.{accession}`) — it is built into its own gene tree and joins the same profile. In the node list, **click a row to add its node path** instead of typing it; in the ETE4 tree preview, **right-click a branch → "Use branch for profiling"** to send it straight to the Node-path list. Any built gene tree can be **downloaded as a plain Newick file**.

Four ways to define the subclades:

| Mode | How subclades are chosen |
|------|--------------------------|
| **Depth slider** | Cut the tree at a chosen root-to-node distance; every branch crossing that depth starts a subclade |
| **Manual MRCA** | Pick groups of leaves; each group's most recent common ancestor becomes one subclade |
| **Node path** | Name internal nodes explicitly by their child-index path from the root |
| **Auto duplication** | Give a taxonomic group (NCBI taxid, e.g. 33213 = Bilateria); nodes where the same species appears on both sides of a split — a duplication signature — become the split points |

Auto duplication works best for families with a handful of ancient paralog groups. For large superfamilies (thousands of leaves) the outermost duplication sits near the root and yields only a coarse 2-way split — use the depth slider there instead.

Species tree with Phylogenetic Profile each taxon for the subclades of interest.

![high-re profile ete tree](figures/profile_ete.png)


---

### Comparative analysis ("concept check")

Ask which protein families are present in one taxon group but absent in another — e.g. *families in **Mucorales** (`4827`) but not in **Human** (`9606`)*. Give each group a taxon or clade ID (a clade is expanded to its member species via NCBI taxonomy), and/or derive a group's taxa from an uploaded tree.

The report gives you:

- **Present in A, absent in B**, the **reverse**, and **shared** — family tables (name, accession, type, taxa count, protein count), each downloadable as CSV.
  ![comparative analysis results](figures/comparative_analysis_results.png)

- **Protein drill-down** — click a family's protein count to list the actual accessions (accession · taxon · organism · best e-value).
- **Heatmap** of the top differential families across the group-A taxa (capped at 60 taxa for readability).
  ![comparative analysis protein drill-down and heatmap](figures/comparative_analysis_heatmap.png)
  
- **Species tree coloured by group** — static PNG *and* the interactive ETE4 explorer (green = A only, amber = B only, blue = both; capped at 150 taxa).
  ![comparative analysis species-tree](figures/comparative_analysis_ncbi_tree.png)
  ![comparative analysis species-tree](figures/comparative_analysis_ete4_tree.png)


---

## Repository structure

```
FORK/
├── app.py                            # Entry point — Flask app
├── run_webapp.sh                     # Convenience launcher
├── templates/                        # HTML templates (Jinja2)
│   ├── base.html
│   ├── index.html
│   ├── about.html
│   ├── tree.html
│   ├── presence.html
│   ├── highres.html
│   ├── profiling.html                # Combined Presence/Absence + High-Res page
│   ├── compare.html                  # Comparative "concept check" page
│   ├── utilities.html
│   └── login.html                    # Log in / register (optional accounts)
├── static/                           # CSS / JS assets
│   └── ete4_overrides/contextmenu.js # Repo-served ETE4 right-click menu override
├── get_reference_uniprot_set_lib.py  # Backend retrieval library (importable)
├── tree_from_db.py                   # CLI: fetch → align → tree → viewer
├── subclade_partition.py             # Partition gene trees into subclades
├── tree_builder.py                   # Per-Pfam tree orchestration + caching
├── ete_profile.py                    # ETE4 viewer — presence/absence on NCBI tree
├── ete_highres_profile.py            # ETE4 viewer — high-res profile on NCBI tree
├── ete_species_tree.py               # ETE4 viewer — comparison species tree (by group)
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
flask, gunicorn, pandas, matplotlib, seaborn, biopython,
mysql-connector-python, python-dotenv, ete4>=4.4.0, PyQt6, numpy
```

CLI tools (**MAFFT**, **trimAl**, **FastTree**, **IQ-TREE**) and the NCBI taxonomy database are
provided by the same environment. All managed via `FORK.yml`.

---

## Notes

- Output files (`.fa`, `.mft`, `.nwk`, `.itol_*.txt`) are written to the path given by `--prefix` and excluded from version control via `.gitignore`.
- The `.env` file contains credentials — never commit it.
- As an alternative for visualizing trees, load the `.nwk` + `.itol_colors.txt` + `.itol_domains.txt` files into [iTOL](https://itol.embl.de).
- ETE4 interactive viewer starts a local server on the first free port in the range 5001–5050 (auto-selected, so multiple trees can be opened at once). The D3 viewer has no server dependency.
- The **Comparative** tab, the species-tree renders, and **auto-duplication** partitioning use ETE4's NCBI taxonomy database (`taxa.sqlite`). It downloads on first use, or build it once with `python -c "from ete4 import NCBITaxa; NCBITaxa().update_taxonomy_database()"`.
- High-res profiling caches trees by build parameters; rerunning with the same settings reuses the cache.
- In the hosted deployment, **Recent runs and results are per-user** (an anonymous browser session, or an optional account) — see **Deployment** above.
- Heatmaps render at print quality with a perceptually-uniform, colorblind-safe colormap.

**[Full API reference & CLI guide](API_REFERENCE.md)**
