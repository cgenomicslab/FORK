# FORK — Session Handoff

A portable summary of the Claude Code session so it can be continued in another
Claude chat (claude.ai) or a fresh Claude Code session. Paste or upload this
file to give the next assistant full context.

- **Project:** FORK — a Flask web app for CGLab reference-proteome analysis
  (phylogenetic trees, high-resolution phylogenetic profiling, presence/absence).
- **Local dev path (Mac):** `/Users/athenamarounka/Documents/master/CG_LAB/local_db/scripts/FORK`
- **Server (bare metal):** `nada-P8`, app checked out at `~/uniprot_scripts/FORK`,
  conda env `bio_tools`.
- **Goal in progress:** host it at **fork.cgenomicslab.org** for testing, using
  **Gunicorn + Nginx + systemd**.

## How to resume the ORIGINAL Claude Code session (Mac)
The full transcript is saved locally as JSONL:
```
~/.claude/projects/-Users-athenamarounka-Documents-master-CG-LAB-local-db-scripts-FORK/08296c9f-5bff-438d-aa07-71f417ea016a.jsonl
```
From the project directory run `claude --resume` (or `/resume` inside Claude Code)
and pick the session — the entire history comes back with context intact. This
file is just a human-readable digest for continuing *elsewhere*.

---

## 1) Code changes made this session (all in the repo)

1. **New logo.** Replaced the wooden-fork SVG with the CGLab green tree image
   (`static/logo.png`). The `#fork-logo` `<symbol>` in `templates/base.html`
   embeds it via `<image href="/static/logo.png">` on a square `viewBox`; all
   brand spots use it via `<use>`. Marks squared to `0 0 48 48`.

2. **Click-to-select node paths** (Node-path partition mode). In
   `templates/highres.html` and `templates/profiling.html`, the "Show node list"
   table rows are clickable — click adds/removes that node's path in the
   textarea (✓ highlight), real `addEventListener` handlers. Styles in
   `static/css/main.css` (`.nodelist-row`, `.nl-check`, `.selected`).

3. **Combine DB Pfams with an uploaded FASTA** (high-res profiling). New
   "External FASTA (optional)" upload; you can supply Pfams, a FASTA, or both.
   FASTA is built into its own gene tree via `tree_from_db.py --local_fasta`.
   - Headers **must** be `{taxid}.{accession}` (e.g. `>9606.P04637`); rejected
     up front otherwise.
   - `tree_builder.py`: `build_one_tree`/`build_trees` take `local_fasta(s)`;
     FASTA content hash folded into the cache key.
   - `get_reference_uniprot_set_lib.py`: `get_highres_profile` /
     `fetch_highres_profile` take optional `acc_to_taxon_override`; external
     sequences attributed via the taxon in the leaf name; **DB lookup always
     wins** so pure-DB results are unchanged. Added `_taxon_names_by_ids`.
   - `app.py`: `build-trees` accepts `combine_fasta`; `compute-profile` builds
     the taxid override; `download-fasta` falls back to uploaded sequences.

4. **"Cannot draw tree with width 0" fix.** ETE4's native dialog for
   zero-width trees (identical sequences / topology-only NCBI trees). In
   `tree_from_db.py` (MODE B explorer), `ete_highres_profile.py`, and
   `ete_profile.py`: after bumping zero branches to 1e-6, if the whole tree is
   still ~zero width, fall back to a uniform-branch dendrogram so it always
   renders. Trees with real branch lengths are left untouched.

5. **Download tree as plain Newick** in every tab that builds a tree.
   - Tree tab (`/tree`): a "Tree (Newick)" bar now shows in **all** viewer modes
     (previously only in static-PNG mode). Bug fixed: the ETE4-server and static
     branches now store `newick`/`tree_path` (via a `_built_newick()` helper),
     and `/api/job` reports `has_newick` for every viewer mode. Download filename
     cleaned to basename (was full path).
   - Profiling / High-res: per-Pfam "Download Newick" button; new endpoint
     `POST /api/highres/download-tree` returns the raw builder Newick.
   - The served Newick is the raw IQ-TREE/FastTree `.treefile` (clean standard
     Newick; not the polytomy-resolved variant).

**Logs updated:** `INTERFACE_CHANGES.md` (Round 5) and `CHANGES.md`.
**All verified** without a DB via the `--local_fasta` path + Flask test client.

### Local test artifacts / commands used
- Tools present on the Mac: mafft, einsi, clustalo, trimal, iqtree (fasttree missing → tested with `--ml iqtree`).
- The `--local_fasta` path builds trees with **no database**, which is how most
  of the above was verified end-to-end.

---

## 2) Deployment plan (Gunicorn + Nginx + systemd)

Target: **fork.cgenomicslab.org**, test purposes, **bare metal** on `nada-P8`
(confirmed: `systemd-detect-virt` = none). No VM (PI suggested one, but going
direct on the host for now).

### ⚠️ The single most important constraint
FORK keeps **job state in memory** (`_jobs`, live ETE4 viewer ports). So:
- **Run ONE Gunicorn worker with threads:**
  `--worker-class gthread --workers 1 --threads 8`. Multiple *worker processes*
  break job polling (a job started on one worker is invisible to another →
  intermittent "job not found", broken ETE4 viewers).
- **The 4 cores** go to the tree-building **subprocesses** (set the UI "CPU
  threads" field to 4), not to web workers.
- `--timeout 300` (launching an ETE4 viewer blocks up to 120 s) and **do NOT**
  set `--max-requests` (a mid-build worker restart would lose running jobs).
- WSGI entry point is `app:app`.

### Phases
1. **Provision / prep host** — Ubuntu, 4 cores for builds, DNS A record for
   fork.cgenomicslab.org → this server; ensure outbound internet + inbound 80/443.
2. **System + bioinformatics deps** — mafft, einsi, clustalo, trimal,
   iqtree/fasttree, headless-Qt libs for ETE4 (`QT_QPA_PLATFORM=offscreen` is
   already set in code).
3. **Code + conda env `bio_tools` + `.env`** — install requirements incl.
   **gunicorn**; point `.env` at the DB; pre-warm NCBI taxonomy
   (`python -c "from ete4 import NCBITaxa; NCBITaxa()"`).
4. **Gunicorn** — 1 gthread worker (config above); smoke-test on 127.0.0.1:8001.
5. **systemd** — service (WorkingDirectory = repo, conda python, Restart=always).
6. **Nginx** — reverse proxy fork.cgenomicslab.org → Gunicorn; raise
   `client_max_body_size` (FASTA/tblout uploads) and `proxy_read_timeout` (long
   polls + 120 s ETE4 launch). Check nothing already uses 80/443 first.
7. **TLS** — Let's Encrypt/certbot (HTTP-01 if reachable on 80; DNS-01 if behind
   NAT/proxy).
8. **End-to-end check** — build a tree, open ETE4 viewer, run combined
   Pfam + FASTA profile, download a Newick.

### Where we are right now — Phase 1 checkpoint (environment audit)
Bare-metal cautions: don't clobber existing services; check 80/443 first; cap
tree-build `--cpu` at 4 on this shared box.

Run and report:
```bash
# OS + resources
lsb_release -d; nproc; free -h; df -h ~
# does bio_tools have deps + tools?
which python gunicorn mafft trimal FastTree fasttree iqtree
python -c "import flask, ete4, mysql.connector; print('core deps OK')"
# gunicorn WSGI entry point importable?
python -c "import app; print('WSGI:', type(app.app))"
# DB reachable? (redact host if pasting)
grep -iE 'host|port' .env; nc -vz <db_host> 3306
# sudo + ports free?
sudo -v && echo sudo OK
sudo ss -tlnp | grep -E ':(80|443)\b' || echo "80/443 free"
# DNS
hostname -I; dig +short fork.cgenomicslab.org
```
Then decide certbot method: HTTP-01 (server directly reachable on 80/443) vs
DNS-01 (behind campus NAT/central proxy).

---

## 3) Immediate reminders for testing the app locally
- **Restart `python app.py`** after any `app.py` change (only *templates*
  auto-reload, not Python).
- **Hard-refresh** pages (`Cmd+Shift+R`) to pick up new template/JS.
- The Newick "Download" bar appears **on build completion**, in all viewer modes.
