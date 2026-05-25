# Setup — One-time scripts

These scripts are run **once by the database administrator** to build and populate the local database. Regular users of the Streamlit tool do not need to run these.

---

## 1. `uniprot_sync_v7.py` — Build the local database

Downloads UniProt Reference Proteomes and populates the MySQL database.

**Run once per UniProt release** (or when you want to add a new version):

```bash
python uniprot_sync_v7.py --version 2026_01
```

This script:
- Downloads reference proteome FASTA files from UniProt FTP
- Parses protein metadata (accession, gene name, organism, taxon ID, proteome ID)
- Deduplicates sequences by MD5 hash
- Populates: `proteins`, `sequences`, `proteomes`, `go_terms`, `protein_go`, `protein_pfam`, `pfam_domains`

Requires: MySQL connection configured in `.env` (see main README).

Expected runtime: several hours for a full release (hundreds of proteomes).

---

## 2. `py_hmmer_hmmsearch.py` — Run HMM search across the database

Scans all protein sequences in the DB against Pfam-A HMMs using PyHMMER.

**Run once after `uniprot_sync_v7.py` completes** for a new version:

```bash
python py_hmmer_hmmsearch.py --version 2026_01 --pfam_db /path/to/Pfam-A.hmm --cpu 32
```

This script:
- Fetches all sequences for the given version from the DB
- Runs PyHMMER `hmmsearch` against Pfam-A with gathering thresholds
- Stores all hits (including multi-domain proteins) in `hmm_search_results`

The `hmm_search_results` table is what powers all HMM-based queries, presence/absence analysis, and phylogenetic tree building in the main tool.

Expected runtime: many hours for ~130M proteins. Designed to run on the lab server with multiple cores.

---

## Database connection

Both scripts read connection parameters from a `.env` file in the same directory:

```
DB_HOST=localhost
DB_USER=your_db_user
DB_PASSWORD=your_password
DB_NAME=uniprot_db_cglab
```
