# Full Evolution Log: from the original Streamlit tool to the current Flask web app

This document compares the **original tool** at
`/local_db/scripts/uniprot-lab-manager/`
against the **current project** at
`/local_db/scripts/uniprot-lab-manager_copy/`
and records every change, addition, and architectural decision made along the way.

---

## 1. Original tool — what existed before

### Entry point and structure

The original is a **Streamlit multi-page application**:

```
uniprot-lab-manager/
├── homepage.py                   ← Streamlit landing page (logo + welcome text)
├── pages/
│   └── 1_UniProt_Lab_Manager.py  ← one-liner: exec(open("retrieval_script_gui.py").read())
├── retrieval_script_gui.py       ← all UI logic (2 229 lines, Streamlit)
├── get_reference_uniprot_set_lib.py
├── viz_utils.py
├── tree_from_db.py
├── subclade_partition.py
├── tree_builder.py
├── interactive_tree_component.py
├── setup/
│   ├── uniprot_sync_v7.py
│   └── pyhmmer_hmmsearch.py
└── uniprot-lab-manager.yml
```

The app is launched with `streamlit run homepage.py`.
The sidebar lets users navigate between pages; the single content page
executes `retrieval_script_gui.py` through the `pages/` router.

### Tabs / sections in the original Streamlit GUI

All ten functional sections live in `retrieval_script_gui.py` and are
selected via a `st.selectbox`:

| Tab | What it does |
|---|---|
| **Standard Retrieval** | Filter sequences by version, taxonomy IDs (typed or file), proteome ID, GO ID, Pfam ID. Returns a dataframe + FASTA download. |
| **HMM Search** | Query by HMM name or accession, optional e-value cutoff, optional taxonomy filter. Returns hits + FASTA download. |
| **Accession Lookup** | Batch fetch by accessions or protein names. Returns records + FASTA download. |
| **Domain Coordinate Lookup** | Fetch domain architecture for a list of accessions. Displays table + draw-on-demand domain architecture diagram (PNG download). |
| **Database Info** | List available UniProt versions and table statistics. |
| **Phylogenetic Tree** | Full tree-building pipeline via `tree_from_db.py`. Three viewer modes: interactive D3 (default), ETE4 static PNG (domain shapes), ETE4 live server (SSH-tunnel iframe). Alignment tool, tree method, trimAl threshold, CPU threads, taxon/Pfam colouring, custom colormap, local FASTA, MSA overlay, and output directory are all configurable. Downloads: Newick, iTOL colours file, iTOL domains file. |
| **GO → Domain Profiles** | Given a GO term, returns all HMM domain profiles associated with proteins annotated with that GO term. |
| **Presence/Absence & Drill-down** | Builds a presence/absence matrix (organisms × Pfam profiles, cell value = protein count). On-demand clustered heatmap. Cell drill-down: sub-profile enrichment table and domain architecture table for the selected (organism, profile) combination. Architecture diagram can be drawn for any selected pattern. |
| **Extract Downloaded Branch** | Upload a Newick branch exported from the ETE4 viewer. Extracts leaf accessions, fetches their sequences from the DB, offers `.txt` (accession list) and `.fasta` downloads. |
| **High-Resolution Phylogenetic Profile** | Multi-step workflow: (1) build one gene tree per Pfam, (2) partition each tree into subclades (depth slider / manual MRCA / node path), (3) assemble a matrix where columns = subclades (not whole Pfam families). Includes: NCBI species tree reference, ETE4 smartview per-Pfam, profile heatmap with binary / log-scale / column-clustering options, subclade FASTA download, per-run cache keying. |

### Supporting modules in the original

| File | Role |
|---|---|
| `get_reference_uniprot_set_lib.py` | All database interaction: `UniProtRetriever` class, `fetch_sequences`, `fetch_sequences_by_hmm_hit`, `fetch_fasta_string`, `fetch_domains_by_accession`, `fetch_domains_by_go`, `fetch_sequences_by_accession`, `fetch_presence_absence_matrix`, `fetch_accessions_for_cell`, `fetch_subprofile_hits`, `fetch_domain_architectures`, `fetch_accessions_with_taxids`, `fetch_highres_profile`. Also a CLI entry point (`_build_parser` / `main`). |
| `viz_utils.py` | `draw_domain_architecture`, `render_tree` (Bio.Phylo → PNG), `draw_presence_absence_heatmap` (seaborn clustermap), `draw_highres_profile_heatmap` (custom matplotlib heatmap with optional log scale and column clustering). |
| `tree_from_db.py` | CLI script: fetches sequences from DB, runs alignment (MAFFT/EINSI/ClustalO), trimAl, FastTree/IQ-TREE, annotates with NCBI lineage, writes Newick + iTOL colour/domain files, optionally launches ETE4 smartview server or renders static PNG. |
| `subclade_partition.py` | `partition_by_depth`, `partition_by_mrca`, `partition_by_node_path`, `list_internal_nodes`, `get_max_root_distance`. |
| `tree_builder.py` | `build_one_tree`, `build_trees` (parallel/sequential tree building with caching), `cache_key`, `strip_leaf_prefix_in_subclades`, `parse_leaf_to_accession`, `parse_leaf_to_taxid`. |
| `interactive_tree_component.py` | D3-based inline Newick viewer: `build_tree_html`, `parse_itol_colors`. |
| `setup/uniprot_sync_v7.py` | Syncs a local UniProt reference proteome database from a downloaded `.tar.gz`. Logs to `update_history.log`. |
| `setup/pyhmmer_hmmsearch.py` | Runs pyHMMER against the local DB to populate HMM hit tables. |

---

## 2. What changed in the copy — Streamlit-level additions

The copy started from an earlier (878-line) snapshot of `retrieval_script_gui.py`
and all ten final tabs were built up during development. The backend modules
(`get_reference_uniprot_set_lib.py`, `viz_utils.py`, `tree_from_db.py`,
`subclade_partition.py`, `tree_builder.py`, `interactive_tree_component.py`)
are **byte-for-byte identical** between the original and the copy, with one
exception (see below).

### Changes to `retrieval_script_gui.py` during development

| Commit | Change |
|---|---|
| `cc78328` | Query extended to include protein names in `fetch_sequences` conditions |
| `1cdb3f1` | DB connection configurable via CLI args (`--host`, `--user`, `--database`) |
| `3acb2e7` | Taxonomy ID input: added file-upload option; output directory for tree made configurable |
| `9ba98ee`, `5089c20`, `ff2fd31` | General UI and retrieval script updates; merge conflicts resolved |
| `82414fc` | ETE4 features added to Phylogenetic Tree tab and `tree_from_db.py` |
| `dab7f4b` | New tab: **Extract Downloaded Branch** (accession extraction from uploaded Newick) |
| `3187e2d` | Extract Downloaded Branch tab extended to also fetch and download sequences |
| `14764e3` | Extract Downloaded Branch: further enhancement |
| `adc1273` | High-res profile tab: **Node path** partition mode — paste ETE4 node path directly |
| `7ba151a` | Presence/Absence tab: checkbox label for clustering columns renamed |
| `641687f` | High-res profile: `ete_highres_profile.py` and `ete_profile.py` integrated |
| `398b4c0` | UI and `viz_utils.py` optimisation |
| `870422c` | ETE improvements in Phylogenetic Tree tab (rendering quality, popup properties) |
| `30afb88` | High-res profile tab: **NCBI Taxonomy species tree** added (section 2b) |
| `973a274` | README updated; NCBI taxa annotation added (`sp_tree.annotate_ncbi_taxa`) |

### The one code difference between original and copy

`retrieval_script_gui.py` in the copy has one additional line vs. the original
(at line 1623 in the species tree rendering block):

```python
sp_tree.annotate_ncbi_taxa(taxid_attr='name')
```

This call annotates species tree nodes with scientific names so they are
displayed in the rendered PNG instead of bare taxon IDs.

### Changes to `tree_from_db.py` CLI flags (vs. early version)

All flags added during development (all present in both original and copy at their
current state):

| Flag | Purpose |
|---|---|
| `--evalue` | E-value cutoff for HMM hits |
| `--no_ncbi` | Skip NCBI lineage annotation (faster) |
| `--no_explore` | Do not launch ETE4 server (run headless) |
| `--exclude_taxids` | Comma-separated taxon IDs to exclude |
| `--output_dir` | Write all output files to this directory |
| `--pfam_source` | `hmmsearch` (local) or `uniprot` (DB annotation) |
| `--pfam_logic` | `or` (any Pfam) / `and` (all Pfams) for multi-Pfam queries |
| `--color_by` | `taxon` or `pfam` branch colouring |
| `--MSA` | Attach MSA to leaves in ETE4 viewer |
| `--positions` | MSA display column range (e.g. `30:60`) |
| `--port` | Port for ETE4 smartview server |
| `--colormap` | Path to custom taxon→colour mapping file |
| `--local_fasta` | Use local FASTA instead of DB fetch |
| `--render_ete_static` | Render a static ETE4 PNG (no server) |
| `--static_layers` | Comma-separated layers: `names,domains,colors,gene,msa` |
| `--use_resolved` | Open the resolved (ladderized) tree in the viewer |

---

## 3. The Flask web application — new in the copy only

The major addition in the copy is a complete **Flask-based web application**
that replaces Streamlit as the user-facing interface. The Streamlit file
(`retrieval_script_gui.py`) is kept for reference but is no longer the
entry point.

### New files added

| File | Lines | Role |
|---|---|---|
| `app.py` | 1 263 | Flask server — all routes and background job management |
| `utils.py` | 32 | Shared colour utility (`color_gradient`) used by ETE rendering scripts |
| `ete_profile.py` | 220 | Standalone ETE4 profile viewer — launched as a subprocess by Flask; takes an HMM tblout file + taxid file + colormap and opens the ETE4 smartview server |
| `ete_highres_profile.py` | 158 | High-res variant of the ETE4 profile viewer for use in the High-res Profiling page |
| `run_webapp.sh` | 6 | Shell launcher: `cd` to project root → `python app.py` |
| `templates/base.html` | — | Shared HTML layout with navigation bar |
| `templates/index.html` | — | Standard Retrieval / home page |
| `templates/tree.html` | — | Phylogenetic Tree page |
| `templates/presence.html` | — | Presence/Absence & Drill-down page |
| `templates/highres.html` | — | High-Resolution Phylogenetic Profile page |
| `templates/utilities.html` | — | Utilities: HMM search, accession lookup, domain lookup, GO profiles, Extract Branch |
| `static/css/main.css` | — | Application stylesheet |
| `static/js/app.js` | — | Frontend logic — all AJAX calls to the Flask API |
| `static/logo.png` | — | App logo — the CGLab green tree mark (embedded by the `#fork-logo` SVG symbol in `base.html`) |

### Why Flask instead of Streamlit

Streamlit re-runs the entire script on every widget interaction. This is
incompatible with long-running operations (tree building, HMM search,
high-res profiling) because the app would block the browser until the job
finished. Flask solves this by dispatching each job to a background thread
and exposing a `/api/job/<job_id>` polling endpoint — the browser submits a
job, receives a job ID, then polls for status/result without blocking.

### Flask pages and their API endpoints

#### Infrastructure (`app.py` lines 40–165)

| Route | Method | Purpose |
|---|---|---|
| `/` | GET | Serve `index.html` |
| `/tree` | GET | Serve `tree.html` |
| `/presence` | GET | Serve `presence.html` |
| `/highres` | GET | Serve `highres.html` |
| `/utilities` | GET | Serve `utilities.html` |
| `/api/db-config` | POST | Set DB connection for the session |
| `/api/db-defaults` | GET | Read `.env` defaults |
| `/api/db-info` | GET | Table statistics |
| `/api/job/<job_id>` | GET | Poll background job status / collect result |

#### Phylogenetic Tree (`/tree`)

| Route | Method | Purpose |
|---|---|---|
| `/api/run-tree` | POST | Launch tree-build pipeline as a background job |
| `/api/tree-html/<job_id>` | GET | Return rendered ETE4 interactive HTML |
| `/api/tree-static-png/<job_id>` | GET | Return static PNG (base64) |
| `/api/download/tree/<job_id>/<filetype>` | GET | Download Newick / PNG / SVG |

#### Presence/Absence & Drill-down (`/presence`)

| Route | Method | Purpose |
|---|---|---|
| `/api/presence-matrix` | POST | Compute presence/absence matrix |
| `/api/drill-down` | POST | Drill into a selected (organism × profile) cell |
| `/api/draw-architecture` | POST | Render domain architecture figure |

#### High-Resolution Phylogenetic Profile (`/highres`)

| Route | Method | Purpose |
|---|---|---|
| `/api/highres/build-trees` | POST | Build gene trees from DB Pfams and/or an uploaded external FASTA (`combine_fasta`, headers `{taxid}.{accession}`) — background job |
| `/api/highres/job/<job_id>` | GET | Poll build job |
| `/api/highres/partition` | POST | Partition a tree into subclades |
| `/api/highres/list-nodes` | POST | List internal nodes (for node-path mode) |
| `/api/highres/tree-html` | POST | Render a subclade as ETE4 interactive HTML |
| `/api/highres/launch-ete4` | POST | Launch live ETE4 smartview server for a subclade |
| `/api/highres/stop-ete4` | POST | Kill the live ETE4 server |
| `/api/highres/compute-profile` | POST | Assemble and return the profile matrix |
| `/api/highres/download-fasta` | POST | Download FASTA for all accessions in a subclade |
| `/api/highres/species-tree` | POST | Build and render NCBI taxonomy species tree |
| `/api/highres/ete-matrix-viz` | POST | Render the profile matrix as an ETE4 figure |
| `/api/ete-profile/launch` | POST | Launch `ete_profile.py` as a standalone subprocess viewer |

#### Utilities (`/utilities`)

| Route | Method | Purpose |
|---|---|---|
| `/api/fetch-sequences` | POST | Standard sequence retrieval |
| `/api/hmm-search` | POST | HMM hit retrieval |
| `/api/accession-lookup` | POST | Batch accession / protein name lookup |
| `/api/domain-lookup` | POST | Domain coordinate lookup |
| `/api/go-domains` | POST | GO term → HMM domain profiles |
| `/api/extract-branch` | POST | Extract accessions + sequences from an uploaded Newick branch |

---

## 4. Summary: what the copy has that the original does not

| What | Where |
|---|---|
| Flask web app with async job management | `app.py` |
| HTML templates for all pages | `templates/` |
| Frontend JS/CSS | `static/js/app.js`, `static/css/main.css` |
| Standalone ETE4 profile viewers | `ete_profile.py`, `ete_highres_profile.py` |
| Colour gradient utility | `utils.py` |
| Shell launcher | `run_webapp.sh` |
| NCBI taxa annotation on species tree nodes | `retrieval_script_gui.py` line 1623 |

All backend Python modules (`get_reference_uniprot_set_lib.py`, `viz_utils.py`,
`tree_from_db.py`, `subclade_partition.py`, `tree_builder.py`,
`interactive_tree_component.py`) are identical in both locations.

---

## 5. How to run

**Streamlit (original interface, still works in the copy):**
```bash
streamlit run homepage.py   # not present in copy — use retrieval_script_gui.py directly
streamlit run retrieval_script_gui.py
```

**Flask web app (new interface, copy only):**
```bash
bash run_webapp.sh
# or directly:
python app.py
```
