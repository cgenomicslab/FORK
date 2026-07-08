# API Reference Guide

This document provides complete documentation for the backend retrieval tools, visualization scripts, and command-line pipelines of the **FORK**.

---

## `get_reference_uniprot_set_lib.py`

### Overview

A retrieval script-interface for the Local UniProt Reference Proteome database. It connects to the local MySQL database and allows querying of protein sequences from any loaded UniProt reference proteome release.

The script serves two purposes simultaneously:

- **Command-Line Interface (CLI):** Run directly from the command line to export FASTA files.
- **Python Library API:** Import as a library into other scripts for programmatic access to sequences, enabling downstream analyses such as multiple sequence alignments (MSA) and phylogenetic tree construction.

> As of release `2026_01`, the database contains **133,566,280 canonical proteins** from **34,230 reference proteomes** across Bacteria, Eukaryota, Viruses, and Archaea. Domain annotations from a complete Pfam-A hmmsearch using pyHMMER are available for retrieval, enabling rapid, HMM-based protein set construction.

---

### Requirements

**Option 1**: Setting up the Conda environment from the uniprot-lab-manager.yml -- Needed if planning to run Streamlit UI tool.

```bash
conda env create -f uniprot-lab-manager.yml
conda activate bio_tools
```

**Option 2**: If you only want to use the get_reference_uniprot_set_lib.py script as a Python library in your own project or as a command-line tool, you can skip the heavy Conda environment and install only the minimal requirements via pip:


#### Python Libraries

- `mysql-connector-python`
- `python-dotenv`
- `biopython`*

```bash
pip install mysql-connector-python python-dotenv biopython
```

#### Database

- MySQL 8.0+
- Tables expected: `proteins`, `sequences` (joined on `seq_id`), `protein_go`, `protein_pfam`, `hmm_search_results`
- User Permissions: Read-only (`SELECT`) privileges.

#### Environment Variables

A `.env` file must be present in the same directory as your execution script:

```env
DB_HOST=localhost
DB_USER=cglab_user
DB_PASSWORD=your_password
DB_NAME=uniprot_db_cglab
```

---

### Command-Line Usage

```bash
python get_reference_uniprot_set_lib.py -version VERSION [filters]
```

#### Required Argument

| Argument | Description |
|---|---|
| `-version VERSION` | UniProt release string (e.g., `2026_01`). |

#### Optional Filters

| Argument | Description |
|---|---|
| `-taxonomy ID [ID ...]` | One or more NCBI Taxonomy IDs. |
| `--proteome-id ID` | UniProt Proteome ID (e.g., `UP000005640`). |
| `--go-id ID` | Gene Ontology term (e.g., `GO:0005634`). |
| `--pfam-id ID` | Pfam domain accession via the `protein_pfam` table (e.g., `PF00870`). |
| `--hmm-name NAME` | Pfam HMM name or accession via local `hmm_search_results` (e.g., `Homeodomain` or `PF00046`). |
| `--evalue FLOAT` | Maximum full-sequence E-value for HMM hits (default: `1e-5`). Applies to `--hmm-name`, `--domains`, and `--go-domains`. |
| `--accession ID [ID ...]` | Fetch full sequences for one or more protein accessions (e.g., `P04637 P10275`). |
| `--domains ID [ID ...]` | Get all Pfam domain hits and coordinates for one or more protein accessions. |
| `--go-domains GO_ID` | Get all HMM profiles found in proteins annotated with a specific GO term. |
| `--output-dir DIR` | Directory where the FASTA file will be written. |
| `--list-versions` | Print all versions present in the database and exit. |
| `--list-proteomes` | List all proteome IDs for the given version and exit. |

> **Note:** All filters are optional and combinable using AND logic.

---

### CLI Examples

#### 1. Database Exploration

List all loaded database versions:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --list-versions
```

List all proteome IDs for a specific version:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --list-proteomes
```

#### 2. Basic Sequence Retrieval

Retrieve the full human reference proteome:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606
```

Retrieve proteins across multiple taxonomy IDs:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606 10090 10116
```

Retrieve the human proteome filtered by a specific Proteome ID:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --proteome-id UP000005640
```

Retrieve human proteins annotated with a specific GO term (e.g. nucleus, GO:0005634):

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606 --go-id GO:0005634
```

Retrieve human proteins carrying a specific Pfam domain (via `protein_pfam` table):

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606 --pfam-id PF00870
```

Retrieve sequences for specific accessions:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --accession P04637 P10275
```

Retrieve a single accession and write the FASTA to a specific directory:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --accession P04637 --output-dir ./fastas
```

#### 3. HMM & Domain-Based Retrieval

Retrieve all proteins with a Homeodomain hit across all taxa:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --hmm-name Homeodomain
```

Retrieve human Homeodomain proteins with a strict E-value cutoff:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --hmm-name Homeodomain -taxonomy 9606 --evalue 1e-10
```

Retrieve by Pfam accession (prefix match — `PF00046` automatically matches `PF00046.36`):

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --hmm-name PF00046
```

Retrieve kinase-domain proteins from Human, Mouse, and Rat:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --hmm-name PF00069 -taxonomy 9606 10090 10116
```

#### 4. Fetching Domain Coordinates

Retrieve all Pfam domain hits and coordinates for a single protein:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --domains P04637
```

Retrieve domain profiles for multiple accessions:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --domains P04637 P10275 Q8BFR5
```

Retrieve domain profiles with an E-value filter:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --domains P04637 --evalue 1e-10
```

#### 5. GO → Domain Profile Lookup

List all HMM profiles found in proteins annotated with a GO term:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --go-domains GO:0005634
```

Same, filtered by E-value:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --go-domains GO:0005634 --evalue 1e-10
```

#### 6. Export Configurations

Combine filters and write output to a target directory:

```bash
python get_reference_uniprot_set_lib.py -version 2026_01 --hmm-name Homeodomain -taxonomy 9606 --output-dir ./fastas
```

#### Output File Format

- **Filename:** `uniprot_<identifier>_<version>.fasta`
- **FASTA Header Format:** `>{taxon_id}.{accession}`

---

### Library Usage

#### Import Statement

```python
from get_reference_uniprot_set_lib import (
    UniProtRetriever,
    get_db_config,
    fetch_sequences,
    fetch_sequences_by_hmm_hit,
    fetch_fasta_string,
    fetch_fasta_string_by_hmm_hit,
)
```

#### Standard Retrieval (One-Liner)

```python
from get_reference_uniprot_set_lib import fetch_sequences

records = fetch_sequences(version="2026_01", taxon_ids=[9606, 10090])
for r in records:
    print(r["accession"], r["organism"])
```

#### HMM-Based Retrieval (One-Liner)

Get sequences for a domain family directly from local `hmm_search_results`:

```python
from get_reference_uniprot_set_lib import fetch_sequences_by_hmm_hit

# Get all Homeodomain proteins across all taxa
records = fetch_sequences_by_hmm_hit(version="2026_01", hmm_query="Homeodomain")

# Human + mouse kinases with strict E-value
records = fetch_sequences_by_hmm_hit(
    version="2026_01",
    hmm_query="PF00069",
    evalue_cutoff=1e-10,
    taxon_ids=[9606, 10090],
)

# All bacterial kinases, excluding a specific clade
records = fetch_sequences_by_hmm_hit(
    version="2026_01",
    hmm_query="PF00069",
    taxon_ids=[2],             # Bacteria
    exclude_taxon_ids=[1224],  # Exclude Proteobacteria
)
```

#### HMM-Based Retrieval → FASTA String (One-Liner)

Pipe sequences directly into BioPython without writing a physical file to disk:

```python
import io
from Bio import SeqIO
from get_reference_uniprot_set_lib import fetch_fasta_string_by_hmm_hit

fasta = fetch_fasta_string_by_hmm_hit("2026_01", "Homeodomain", taxon_ids=[9606])
seqs = list(SeqIO.parse(io.StringIO(fasta), "fasta"))
```

#### Fetch Sequences by Accession

```python
from get_reference_uniprot_set_lib import fetch_sequences_by_accession

# Single accession
records = fetch_sequences_by_accession("2026_01", "P04637")

# Multiple accessions — auto-chunked for large lists
records = fetch_sequences_by_accession("2026_01", ["P04637", "P10275", "Q8BFR5"])
```

#### Fetch Domain Coordinates by Accession

```python
from get_reference_uniprot_set_lib import fetch_domains_by_accession

# All Pfam domains for p53 and androgen receptor
domains = fetch_domains_by_accession("2026_01", ["P04637", "P10275"])
for d in domains:
    print(d["accession"], d["hmm_name"], d["ali_from"], d["ali_to"], d["full_evalue"])

# With strict E-value filter
domains = fetch_domains_by_accession("2026_01", ["P04637"], evalue_cutoff=1e-10)
```

Each result dict contains: `accession`, `protein_name`, `organism`, `taxon_id`, `proteome_id`, `hmm_name`, `hmm_accession`, `domain_number`, `domain_count`, `full_evalue`, `full_score`, `domain_evalue`, `domain_score`, `ali_from`, `ali_to`, `hmm_from`, `hmm_to`, `env_from`, `env_to`.

#### Fetch HMM Profiles by GO Term

```python
from get_reference_uniprot_set_lib import fetch_domains_by_go

# All Pfam profiles found in nuclear proteins (GO:0005634)
results = fetch_domains_by_go("2026_01", "GO:0005634")
for r in results:
    print(r["hmm_name"], r["hmm_accession"], r["protein_count"])

# With E-value filter
results = fetch_domains_by_go("2026_01", "GO:0005634", evalue_cutoff=1e-5)
```

#### Advanced Retrieval: Presence/Absence & Drill-Down (One-Liners)

The library exposes the backend logic used in the GUI's Presence/Absence tab, allowing you to programmatically build complex analytical matrices and domain architectures.

| Function | Description |
|---|---|
| `fetch_presence_absence_matrix(version, pfam_queries, ...)` | Step 1: Returns a list of dicts mapping taxa to domain counts. Ready to be converted into a pandas DataFrame and pivoted. |
| `fetch_accessions_for_cell(version, pfam_query, taxon_id, ...)` | Bridge: Returns the specific accessions that make up a single cell in your matrix. |
| `fetch_subprofile_hits(version, accessions, ...)` | Step 2A: Given a list of accessions, returns all deeper HMM profiles (TIGRFAMs, sub-Pfams) enriched in that set. |
| `fetch_domain_architectures(version, accessions, ...)` | Step 2B: Characterizes the full left-to-right domain architecture patterns for a set of proteins. |

**Example: Presence/Absence Matrix**

```python
from get_reference_uniprot_set_lib import fetch_presence_absence_matrix
import pandas as pd

rows = fetch_presence_absence_matrix(
    version="2026_01",
    pfam_queries=["Homeodomain", "PF00069"],
    taxon_ids=[9606, 10090, 7227],
    evalue_cutoff=1e-5,
)
df = pd.DataFrame(rows)
# pivot to get taxa as rows, profiles as columns
matrix = df.pivot_table(index="taxon_id", columns="hmm_name", values="protein_count", fill_value=0)
```

**Example: Two-Tier Drill-Down Workflow**

```python
from get_reference_uniprot_set_lib import (
    fetch_accessions_for_cell,
    fetch_domain_architectures,
    fetch_subprofile_hits,
)

# 1. Get accessions for Human (9606) proteins with a Homeodomain
cell_data = fetch_accessions_for_cell("2026_01", "Homeodomain", 9606)
acc_list = [r["accession"] for r in cell_data]

# 2A. Find co-occurring sub-profiles in this set (TIGRFAMs, sub-Pfams, etc.)
subprofiles = fetch_subprofile_hits(
    "2026_01", acc_list,
    exclude_queries=["Homeodomain", "PF00046"],  # hide the original broad query
)
for s in subprofiles[:5]:
    print(s["hmm_name"], s["protein_count"], s["coverage"])

# 2B. Get the full domain architectures for these proteins
archs = fetch_domain_architectures("2026_01", acc_list)
for a in archs[:5]:
    print(f"{a['architecture']} -> {a['protein_count']} proteins ({a['coverage']:.0%})")
# Homeodomain+PBC → 38 proteins (51%)
# Homeodomain     → 12 proteins (16%)
# TALE+Homeodomain → 4 proteins (5%)

# Use collapse_repeats=False when repeat count is biologically meaningful
archs_full = fetch_domain_architectures("2026_01", acc_list, collapse_repeats=False)
```

#### Context Manager (Recommended for Multiple Queries)

```python
from get_reference_uniprot_set_lib import UniProtRetriever, get_db_config

with UniProtRetriever(get_db_config()) as db:
    # Retrieve by HMM hit
    records = db.get_proteins_by_hmm_hit(
        version="2026_01",
        hmm_query="Homeodomain",
        evalue_cutoff=1e-5,
        taxon_ids=[9606, 10090],
    )

    # Convert to BioPython SeqRecord objects for downstream pipelines
    seqrecords = db.to_biopython(records)

    # Export directly to FASTA file
    db.export_fasta(records, "2026_01", "Homeodomain_human_mouse", output_dir="./fastas")
```

#### Full Tree-Building Pipeline Example

```python
from get_reference_uniprot_set_lib import UniProtRetriever, get_db_config

with UniProtRetriever(get_db_config()) as db:
    # 1. Retrieve all bacterial proteins with a specific kinase domain
    records = db.get_proteins_by_hmm_hit(
        version="2026_01",
        hmm_query="PF00069",   # kinase domain
        evalue_cutoff=1e-5,
        taxon_ids=[2],         # NCBI taxid for Bacteria
    )
    print(f"Retrieved {len(records):,} sequences")

    # 2. Convert to BioPython SeqRecord
    seqrecords = db.to_biopython(records)

    # 3. Write to FASTA for alignment
    db.export_fasta(records, "2026_01", "kinase_bacteria", output_dir="./fastas")

    # 4. Pass seqrecords directly into alignment/tree pipeline
    # e.g. with MAFFT via subprocess, or ETE3, or DendroPy
```

#### Overriding Credentials at Runtime

```python
config = get_db_config(host="10.0.0.5", user="myuser", password="mypass")
with UniProtRetriever(config) as db:
    records = db.get_proteins(version="2026_01", taxon_ids=[6087])
```

---

### Architecture

```
INPUT ARGUMENTS (-version, -taxonomy, --hmm-name, etc.)
                           |
                           V
CONNECT to MySQL (via UniProtRetriever or convenience function)
                           |
                           V
BUILD QUERY
  ├── Standard: proteins JOIN sequences
  │   Optional joins: protein_go, protein_pfam
  │   WHERE clauses from provided filters (AND logic)
  │
  └── HMM-based: hmm_search_results JOIN proteins JOIN sequences
      WHERE hmm_name or hmm_accession LIKE query
      AND full_evalue <= cutoff
      SELECT DISTINCT (deduplicates multi-domain proteins)
                           |
                           V
EXECUTE QUERY → fetchall() → list of dicts
                           |
                           V
OUTPUT (choose one or combine)
  ├── to_fasta_string()         → in-memory FASTA string
  ├── to_biopython()            → list of BioPython SeqRecord objects
  └── export_fasta()            → write .fasta file to disk
```

---

### API Reference

#### Convenience Functions (Class-less)

| Function | Description |
|---|---|
| `fetch_sequences(...)` | Standard filter retrieval → list of dicts |
| `fetch_sequences_by_hmm_hit(...)` | HMM hit retrieval → list of dicts |
| `fetch_fasta_string(...)` | Standard filter retrieval → FASTA string |
| `fetch_fasta_string_by_hmm_hit(...)` | HMM hit retrieval → FASTA string |
| `fetch_domains_by_accession(...)` | Fetch Pfam domain coordinates/scores for specific protein accessions |
| `fetch_domains_by_go(...)` | Fetch all HMM profiles found in proteins annotated with a specific GO term |
| `fetch_sequences_by_accession(...)` | Fetch full sequence records by UniProt accession (auto-chunked for large lists) |
| `fetch_presence_absence_matrix(...)` | Build taxa vs. profile matrix (Step 1 of drill-down) |
| `fetch_accessions_for_cell(...)` | Get accessions for one matrix cell (Bridge to Step 2) |
| `fetch_subprofile_hits(...)` | Resolve sub-profiles / TIGRFAMs for a set of accessions (Step 2A) |
| `fetch_domain_architectures(...)` | Get left-to-right domain patterns for a set of accessions (Step 2B) |

---

#### `get_db_config(host, user, password, database)`

Builds the database configuration dictionary from environment variables or explicit arguments. Explicit arguments take priority over `.env` values, which take priority over lab defaults.

---

#### `class UniProtRetriever`

##### Connection & Meta

- `connect()` / `close()`: Manual connection management. Not needed when using the `with` context manager.
- `list_available_versions()`: Prints and returns a summary of all UniProt versions in the database. Returns `list[dict]` with keys `version`, `protein_count`, `taxon_count`, `proteome_count`.
- `get_proteome_ids(version)`: Returns a list of all UniProt Proteome IDs for the given version.

##### Sequence Retrieval

- `get_proteins(version, taxon_ids=None, proteome_id=None, go_id=None, pfam_id=None)`: Standard retrieval with optional filters. All filters combinable with AND logic. Returns `list[dict]` with keys `accession`, `name`, `organism`, `taxon_id`, `proteome_id`, `sequence`.

- `get_proteins_by_accession(version, accessions, chunk_size=5000)`: Fetch full protein records for specific accessions. Large lists are chunked automatically to prevent SQL bottlenecks. Accepts a single string or a list.

- `get_proteins_by_hmm_hit(version, hmm_query, evalue_cutoff=None, taxon_ids=None, exclude_taxon_ids=None)`: Retrieves sequences for proteins with a hit to a given Pfam HMM profile. `hmm_query` matches by exact name or prefix accession (e.g., `PF00046` matches `PF00046.36`). Results are deduplicated — a protein with five Homeodomain repeats appears once. `exclude_taxon_ids` accepts a list of NCBI TaxIDs to remove from results.

##### Domain & GO Analysis

- `get_domains_by_accession(version, accessions, evalue_cutoff=None, chunk_size=5000)`: Returns one row per domain occurrence (including coordinates and scores) for the given accessions. A protein with five Ankyrin repeats produces five rows. Keys: `accession`, `protein_name`, `organism`, `taxon_id`, `proteome_id`, `hmm_name`, `hmm_accession`, `domain_number`, `domain_count`, `full_evalue`, `full_score`, `domain_evalue`, `domain_score`, `ali_from`, `ali_to`, `hmm_from`, `hmm_to`, `env_from`, `env_to`.

- `get_domains_by_go(version, go_id, evalue_cutoff=None)`: Returns all HMM profiles found in proteins annotated with a specific GO term, with counts. Keys: `hmm_name`, `hmm_accession`, `protein_count`.

##### Presence/Absence & Drill-Down

- `get_presence_absence_matrix(version, pfam_queries, taxon_ids=None, evalue_cutoff=None)`: Builds a matrix mapping taxa to domain counts. `pfam_queries` is a list of Pfam names or accessions. Returns `list[dict]` with keys `taxon_id`, `scientific_name`, `hmm_name`, `hmm_accession`, `hmm_type`, `protein_count`.

- `get_accessions_for_cell(version, pfam_query, taxon_id, evalue_cutoff=None)`: Returns the accessions that make up a single matrix cell (one taxon + one profile). Returns `list[dict]` with keys `accession`, `protein_name`, `organism`, `taxon_id`, `proteome_id`. No sequences — use `get_proteins_by_accession()` if sequences are needed.

- `get_subprofile_hits(version, accessions, evalue_cutoff=None, exclude_queries=None, chunk_size=5000)`: Returns deeper HMM profiles (TIGRFAMs, sub-Pfams) enriched in the input set, sorted by `protein_count` descending. `exclude_queries` accepts a list of profile names/accessions to suppress (e.g., hide the original broad Pfam so sub-profiles stand out). Keys: `hmm_name`, `hmm_accession`, `hmm_type`, `protein_count`, `coverage`, `best_evalue`, `best_score`.

- `get_domain_architectures(version, accessions, evalue_cutoff=None, collapse_repeats=True, chunk_size=5000)`: Characterizes the ordered, left-to-right domain architecture patterns. With `collapse_repeats=True` (default), `ANK+ANK+ANK+KH` becomes `ANK+KH` — use `False` when exact repeat count is biologically meaningful. Keys: `architecture`, `arch_accessions`, `protein_count`, `coverage`, `example_accessions` (up to 5).

##### Output Formatting

- `to_fasta_string(records)`: Converts records to a FASTA-formatted string in memory.
- `to_biopython(records)`: Converts records to a list of BioPython `SeqRecord` objects. Requires BioPython.
- `export_fasta(records, version, identifier, output_dir, filename)`: Writes records to a FASTA file.

---

### Notes

- All records returned are plain Python dicts — no special dependencies to access or iterate.
- `to_biopython()` is the recommended bridge to alignment tools (MUSCLE, MAFFT) and tree-building libraries (ETE3/ETE4).
- When issuing multiple queries in the same script, always prefer the context manager to avoid connection overhead.
- The database stores only **canonical sequences**. Isoforms are not present.
- `get_proteins_by_hmm_hit()` queries `hmm_search_results` which was populated by a full Pfam-A hmmsearch using **gathering thresholds**. Results are therefore already filtered at the profile-specific trusted cutoff — `evalue_cutoff` provides an additional filter on top of that.
- For large taxon sets, consider filtering by HMM hit or GO term first to reduce the result set before loading sequences into memory.

---

## `viz_utils.py`

### Overview

A standalone visualization library that generates publication-ready diagrams from database queries. All functions return a `BytesIO` PNG image buffer. This means no files are written to disk, making it perfect for Streamlit dashboards or rendering directly inside Jupyter Notebooks.

### Import

```python
import viz_utils as viz
```

### API Reference

| Function | Description |
|---|---|
| `draw_domain_architecture(domain_records, ...)` | Draws horizontal protein bars with colored domain blocks. Colors are deterministic based on domain name (MD5 hash). |
| `render_tree(newick_string, ...)` | Renders a static phylogenetic tree using `Bio.Phylo` and matplotlib. Auto-scales figure height based on leaf count. |
| `draw_presence_absence_heatmap(matrix_df, ...)` | Draws a seaborn clustered heatmap from a pivoted presence/absence pandas DataFrame. |

### Example: Rendering in Jupyter Notebook

```python
from IPython.display import Image
import viz_utils as viz
from get_reference_uniprot_set_lib import fetch_domains_by_accession

# Fetch domain data
domains = fetch_domains_by_accession("2026_01", ["P04637", "P10275"])

# Generate diagram buffer
img_buffer = viz.draw_domain_architecture(domains, title="P53 Domain Architecture")

# Display in notebook
Image(data=img_buffer.getvalue())
```

---

## `tree_from_db.py`

### Overview

A robust command-line wrapper for generating phylogenetic trees directly from local database queries. It automatically fetches sequences, aligns them (MAFFT, E-INS-i, ClustalO), infers a maximum likelihood tree (FastTree, IQ-TREE), and generates color strips for iTOL.

Optionally launches an interactive ETE4 local web server to explore the tree.

### Command-Line Usage

```bash
python tree_from_db.py --pfam PFAM_NAME --version VERSION --prefix PREFIX [options]
```

### Required Arguments

| Argument | Description |
|---|---|
| `--pfam` | Pfam domain name or accession (comma-separated for multiple). |
| `--version` | UniProt database version (e.g., `2026_01`). |
| `--prefix` | Prefix for all generated output files. |

### Key Optional Arguments

| Argument | Description |
|---|---|
| `--taxids` | Comma-separated list of TaxIDs (or path to a file). |
| `--exclude_taxids` | Comma-separated list of TaxIDs to exclude. |
| `--evalue` | Maximum E-value cutoff. |
| `--aln` | Alignment software to use (`mafft` [default], `einsi`, `clustalo`). |
| `--ml` | Tree building software to use (`fasttree` [default], `iqtree`). |
| `--cpu` | Number of threads to use (default: `32`). |
| `--no_explore` | Skips the ETE4 interactive server launch. |

### Example: Full Tree Pipeline

Fetch all Homeodomain proteins in Human and Mouse, align with MAFFT, build with FastTree, and skip the ETE4 viewer:

```bash
python tree_from_db.py \
  --pfam Homeodomain \
  --version 2026_01 \
  --taxids 9606,10090 \
  --prefix human_mouse_homeo \
  --aln mafft \
  --ml fasttree \
  --no_explore
```

**Outputs generated:**

| File | Description |
|---|---|
| `human_mouse_homeo.fa` | Raw sequences |
| `human_mouse_homeo.mft` | Alignment |
| `human_mouse_homeo.mft.gt01.lg.fasttree` | Newick Tree |
| `human_mouse_homeo.mft.gt01.lg.fasttree.itol_colors.txt` | Drag-and-drop iTOL color strip |
