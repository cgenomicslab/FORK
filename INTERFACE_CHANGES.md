# FORK — Green interface (test copy)

This folder is a **full copy** of `uniprot-lab-manager_copy` with a **reskinned interface**
plus a small, read-only data layer so the Overview dashboard shows **real** activity.

Created: 2026-07-06 · Updated: 2026-07-09
(real KPIs; Recent runs driven by the job registry; wider ETE4 explorers; **renamed to
FORK**; About tab; colormap + partition-lock + ETE4 domain-track fixes; **new CGLab tree
logo**; **click-to-select node paths**; **combine DB Pfams with an uploaded FASTA**)

---

## Round 5 — new logo, click-to-select node paths, DB + FASTA profiling (2026-07-09)

- **New logo: the CGLab green tree mark.** Replaced the wooden three-tined fork SVG (Round 4)
  with the uploaded circular green tree image (`static/logo.png`, 207×207). The
  `<symbol id="fork-logo">` in `base.html` now embeds the PNG via
  `<image href="/static/logo.png">` on a square `viewBox="0 0 48 48"`; every brand spot
  (sidebar, top bar, bottom-left wordmark, About hero) still references it through `<use>`, so
  the single symbol edit updates them all. The four `.fork-mark` usages were squared to
  `0 0 48 48` so the round logo isn't stretched.

- **Click-to-select node paths** (`highres.html`, `profiling.html`, `main.css`). In the
  **Node path** partition mode, the "Show node list" table is now interactive: clicking a row
  **adds or removes** that node's path in the "Node paths" textarea — no copy-paste. Selected
  rows highlight with a ✓, and the highlight also tracks manual edits to the textarea. Rows get
  real `addEventListener('click', …)` handlers (not inline `onclick`). New
  `.nodelist-row` / `.nl-check` / `.selected` styles in `main.css`. The root node (empty path)
  is shown as `(root)` and is non-clickable. Applied to **both** the profiling and high-res
  pages, which share this partition UI.

- **Combine DB Pfams with an uploaded FASTA** (`profiling.html`, `highres.html`, `app.py`,
  `tree_builder.py`, `get_reference_uniprot_set_lib.py`). The high-res build form has a new
  **External FASTA (optional)** upload. You can supply Pfams, a FASTA, or **both**: the FASTA
  is built into its own gene tree via `tree_from_db.py --local_fasta`, appears as its own
  partition block alongside the DB Pfams, and its subclades join the same profile matrix.
  - **Headers must be `{taxid}.{accession}`** (e.g. `>9606.P04637`) — that's how each sequence
    is attributed to a taxon. If none of the headers match, the build is rejected up front with
    a clear error instead of silently dropping sequences.
  - `tree_builder.py`: `build_one_tree` / `build_trees` accept `local_fasta(s)`; the FASTA
    **content hash** is folded into the cache key, so re-uploading different content rebuilds
    instead of reusing a stale cached tree.
  - `get_reference_uniprot_set_lib.py`: `get_highres_profile` / `fetch_highres_profile` accept
    an optional `acc_to_taxon_override`. Sequences absent from the DB are attributed to the
    taxon carried in the tree leaf name; **the DB lookup always wins**, so pure-DB results are
    unchanged. Added `_taxon_names_by_ids` so external taxa still get organism names.
  - `app.py`: `build-trees` accepts `combine_fasta` (validates headers, sanitizes a label from
    the filename, dedupes against Pfams, stores the raw FASTA for the download step);
    `compute-profile` builds the taxid override from the leaf names; `download-fasta` falls back
    to the uploaded sequences for accessions not in the DB.

**Verified.** New-logo markup served on every page and `/static/logo.png` returns
`200 image/png`. Click-to-select logic unit-tested (add / add / toggle-off with correct
highlight) and the served pages carry the real click handlers. FASTA path exercised
end-to-end through the Flask endpoints **with no database** (the `--local_fasta` path needs
none): bad headers → HTTP 400; a FASTA-only build completes, the entry is labelled and
partitions by node path. The profile merge was unit-tested with a stubbed DB — a combined
DB + FASTA matrix attributes every taxon correctly and `missing` empties out, while the
**no-override path is byte-for-byte identical to the old behaviour**. The one path needing a
live DB (`compute-profile` / `download-fasta` round-trip) was validated against stubbed DB
responses; run one real combined build against MySQL to confirm end-to-end.

---

## Round 4 — renamed to FORK (2026-07-08)

- **Name changed to FORK everywhere.** Every user-facing occurrence of the two previous names
  was replaced with **FORK** — page titles, the brand/wordmark in `base.html`, the footer, the
  About page, `README.md`, `API_REFERENCE.md`, docstrings/comments, and this log. Repo/file
  names and URLs (e.g. `uniprot-lab-manager.yml`, the GitHub remote) were left untouched.
- **New logo: a wooden three-tined fork** next to the name. Defined once as an SVG
  `<symbol id="fork-logo" viewBox="0 0 40 48">` at the top of `base.html` and reused in all
  three brand spots (sidebar, top bar, bottom-left) plus the About hero via `<use>`, styled by
  `.fork-mark` in `main.css` (wood brown `#8a5a2b`). Replaced the previous wave/"P" marks.

---

## Round 2 — rebrand + bug fixes (2026-07-07)

Requested batch of fixes, all in the green copy only:

1. **External colormap in the ETE4 tree** (`tree_from_db.py`). Branch colouring only did an
   *exact* taxid match, so a colormap that uses clade/group taxids (e.g. 9443 Primates)
   coloured nothing. Now each leaf taxid is resolved to a colour by walking its NCBI lineage
   (exact match first, then the first matching ancestor) — the same approach the high-res
   profile viewer already uses. Colormap parsing is also tolerant now (skips blank/comment/
   malformed lines). `evoltype` is copied into `node.props` after annotation so the
   **ortho/para** layer always has data (it's an on-demand toggle in ETE4's layer panel).
2. **Footer** (`base.html`). Removed **EggNOG**; **Pfam → InterPro** (the old `pfam.xfam.org`
   link was dead, which was the "broken click"); links now open reliably.
3. **Rebrand to FORK** (`base.html`, `main.css`). New name with a minimal wave mark,
   shown top-left (sidebar brand + top bar) and bottom-left (sidebar wordmark). Above SYSTEM
   sits the **real interactive CGLab pixel logo** — the `<canvas>` cursor-repel animation
   reproduced verbatim from cgenomicslab.org (pixels flee the mouse and spring back), linking
   to <https://cgenomicslab.org/>. The same logo appears on the About page. Page titles and
   footer updated to FORK. (The earlier static four-squares PNG was removed.)
4. **About tab** (`templates/about.html`, `/about` route). Describes the lab (from
   cgenomicslab.org) and the tool (from the README's module table).
5. **Subclade partition locking** (`profiling.html`). Locking now always recomputes the
   partition from the current inputs and returns the parts directly, instead of depending on a
   fragile cached preview — so node-path and auto-duplication partitions lock reliably with or
   without clicking Preview, and a failed lock now shows a clear message instead of silently
   doing nothing.
6. **ETE4 domain rendering — REVERTED.** Attempts to change the domain track caused an ETE4
   draw-time 500 ("when drawing"). The domain code (`_draw_domains`, aliases, and the MSA
   `column=4`) is now **byte-identical to the original** — one labelled `RectFace` per domain,
   unchanged. Verified by diff against the untouched copy and by hitting the ETE4 draw endpoint
   (HTTP 200). Nothing else in `tree_from_db.py` was reverted (the colormap/ortho-para fixes in
   item 1 remain; they run before drawing and don't affect it).
7. **Continuous profiling flow** (`profiling.html`) — see the interface-files table below.

---

## Round 3 — dashboard restore + robustness (2026-07-08)

- **Auto port selection** (`app.py`). The app tries 8080, then the next free port if it's
  taken, and prints the chosen URL — no file editing, multiple users can run side by side.
  `PORT=8090 python app.py` forces a base port. (The app-port picker is `_pick_app_port`, kept
  separate from the existing `_find_free_port` that assigns ETE4 viewer ports.)
- **"Profile" on a Recent run restores that build** (`index.html`, `profiling.html`,
  `recent_runs.py`). A finished **high-res** run now links to `/profiling?job=<id>`, which
  reveals the high-res stage and jumps straight to the **partition step** with the build's
  tree summary — the already-built trees are reused, no rebuild. Single-tree runs link to
  `/tree` ("Open"). Works while the run is still loaded in memory (same app session); after a
  restart the row remains but shows a clear "no longer loaded — re-run" message (the in-memory
  ETE tree objects aren't serialized). The recent-runs entry now carries the job `type` so the
  dashboard can route correctly.
- **Logo = the plain black "P" of FORK, with a cartoon wave wrapping around it.** No
  separate icon: the first letter of the **FORK** wordmark is the ordinary black serif P,
  and a blue cartoon wave (foam-crest curl top-left, water sweeping under and curling up the
  right, white foam + spray droplets) wraps around it; the rest of the word ("hyloWave") is
  normal text. Defined once as `<symbol id="pw-logo" viewBox="0 0 50 56">` (top of `base.html`)
  — the wave drawn first, then the black serif `<text>P</text>` on top — and reused inline via
  `<use>` in all three brand spots as `.pw-letter` (`height: 1.15em`, baseline-aligned). Wave
  blue `#2f86d6`, foam white; the P matches the wordmark (serif, weight 600, ink `#1A1C1A`).
- **Templates auto-reload** (`app.py`). `TEMPLATES_AUTO_RELOAD` is on, so edits to templates
  (e.g. the logo in `base.html`) appear on a normal browser reload without restarting the app.
- **ETE4 domains restored to the original** (`tree_from_db.py`). All domain-track experiments
  were rolled back — `_draw_domains`, the aliases, and the MSA `column=4` now match the
  untouched original exactly (one labelled RectFace per domain). This fixed the "when drawing"
  500 that the changed versions introduced.
- **Removed the top search box** (`base.html`). The `⌘K` field next to **Database** was
  visual-only, so it's gone; the Database button and DB panel are unchanged.
- **Live build progress on the Phylogenetic Tree page** (`app.py`, `tree.html`). The build
  subprocess's output now streams into the job log for **all** viewers — previously the ETE4
  and static viewers ran silently (`Popen` with no capture / a blocking `subprocess.run`), so
  a long build showed no feedback. The status headline now reads **"Building tree…"** and
  updates live with the latest progress line (fetching → aligning → building tree → starting
  ETE4). The d3 viewer already streamed; this brings ETE4 and static in line.

---

## How to run

Same as the original — from this folder:

```bash
python app.py
# it prints the URL it chose, e.g. http://localhost:8080
```

> **Auto port selection.** The app tries **8080** first; if that's taken (e.g. another
> user already has an instance up), it automatically moves to the next free port (8081,
> 8082, …) and prints the chosen URL — **no more editing the file**. To force a specific
> base port, set `PORT`: `PORT=8090 python app.py`.

It reads the same `.env` and talks to the same MySQL database as the original.

---

## What changed

### Interface files
| File | Change |
|------|--------|
| `static/css/main.css` | Re-themed from purple to a phylogenetic-green palette; replaced the dark top header with a **sidebar + top-bar app shell**; added serif/mono type roles and dashboard component styles (KPI tiles, module cards, runs table). Also **widened the ETE4 explorers**: content cap 1180 → 1600px, filter column 320 → 300px, page padding trimmed, and the ETE4 iframe height 800 → 900px (applies to tree / profiling / highres). |
| `templates/base.html` | Rebuilt the page frame as **left sidebar + top bar + main column**. All Jinja blocks (`title`, `page_header`, `content`, `scripts`) are unchanged, so every page still fills them the same way. The DB config panel markup and **all element IDs are identical** to the original. |
| `templates/index.html` | Home rebuilt as an **Overview dashboard**: KPI strip, module cards, and a **real** recent-runs table (fetches `/api/recent-runs`). |
| `templates/profiling.html` | **Continuous flow.** Replaced the two hide-each-other sub-tabs (Presence/Absence vs. High-Resolution) with one continuous page: Stage 1 (Presence/Absence) on top, then a **"Continue to High-Resolution Profiling →"** button that reveals Stage 2 below and **carries the typed selections forward** (version, Pfam, taxonomy IDs, e-value). No analysis behaviour changed — same forms, same routes, same element IDs; only the navigation between the two stages. |
| `static/js/app.js` | One line: the active-nav highlighter now also matches sidebar links (`.sidebar a.nav-item`). |

### Data layer for the dashboard
| File | Change |
|------|--------|
| `recent_runs.py` | **New file.** Turns the app's job registry into dashboard rows and mirrors them into a small history file (`.run_history.json`) so runs persist across restarts. Touches no analysis output. |
| `.run_history.json` | **New file, auto-created.** Tiny bookkeeping store of recent runs (max 100). Safe to delete; it just clears the Recent-runs list. |
| `app.py` | **Additive only:** `import recent_runs as rr`; a new `GET /api/recent-runs` route; a `created` timestamp + `meta` block on each job in `_new_job(...)`; and that `meta` (family, method, output path) is filled at the two build call sites (tree + high-res). No existing route or logic was modified. |

### Design tokens (the "green" identity)
- **Accent** `#1C6B54` (phylogenetic green), hover `#155744`, soft `#E7F0EB`
- **Paper** `#F6F7F4`, cards `#FFFFFF`, ink `#1A1C1A`, muted `#64685F`, hairline `#E4E5DF`
- **Semantic:** done = green, building = amber `#B9722E`, error = red (kept distinct from the accent)
- **Type:** serif (Iowan Old Style / Palatino) for titles & big numbers; system sans for UI;
  monospace for identifiers (Pfam IDs, leaf counts)

---

## What's real now (no more placeholders)

**KPI strip** — all four tiles are live:
- **Reference versions** — count from `/api/db-info`
- **Latest version** — newest version string from `/api/db-info`
- **Database** — real Online/Offline status
- **Completed runs** — real count of finished runs (from `/api/recent-runs`)

**Recent runs table** — real, and now driven by the **job registry** instead of a fixed
folder scan. Every analysis the app launches (Phylogenetic Tree *and* high-res builds) is
registered with its family, method, and output path, so it appears here **no matter which
output path you chose** — the moment it starts:
- **building** while it runs, then **done** (or **failed**) when it finishes.
- **Family** (Pfam), **Method** (aligner → tree method), **Built** (relative time), **Status**.
- **Leaves**: filled with the real count once the tree file is available (single trees show
  the tree's leaf count; high-res shows the combined total). Shows `—` while that isn't known
  yet — e.g. for the live ETE4 viewer, which doesn't hand back a tree-file path.

Runs are mirrored to `.run_history.json`, so they **survive an app restart** too.

### Why the list was empty before
The first version only scanned `/tmp/highres_runs`. The **Phylogenetic Tree** page writes to
*your chosen* output directory (e.g. `scripts/test_results`), so those builds were never in
the scanned folder and never showed up. The registry-based approach fixes that for any path.

---

## Still placeholder / not wired

- The **⌘K search box** in the top bar is visual only (no search backend yet).
- Recent runs deliberately does **not** show taxon input or partition settings: those aren't
  reliably known at the job level, so showing them would be guessing. Columns are limited to
  what's truthfully known (Family, Leaves, Method, Built, Status).

---

## What stayed exactly the same

- Analysis Python left untouched: `subclade_partition.py`, `viz_utils.py`, `utils.py`,
  `ete_profile.py`, `ete_highres_profile.py`, everything in `setup/`.
  **Exceptions from Round 5** (all additive / backward-compatible): `tree_builder.py` and
  `get_reference_uniprot_set_lib.py` gained the optional `local_fasta` / `acc_to_taxon_override`
  paths — with no override or FASTA, behaviour is identical to before — and `app.py` gained the
  `combine_fasta` handling in `build-trees` / `compute-profile` / `download-fasta`.
  Aside from those, `app.py` remains additive (the `/about` route plus the earlier recent-runs
  helpers). `tree_from_db.py` was edited **only** for the Round-2 fixes above (colormap lineage
  resolution + tolerant parsing, `evoltype`→props, and the binned ETE4 domain track); the
  tree-building / alignment / static-render paths are unchanged, and its existing `--local_fasta`
  flag (already present) is what Round 5 drives.
- Other templates (`tree.html`, `presence.html`, `utilities.html`) — not edited; they inherit
  the green look via shared CSS. `highres.html` and `profiling.html` got the Round-5
  click-to-select node list and External-FASTA upload (their existing analysis forms, routes,
  and element IDs are otherwise unchanged); `profiling.html` also carries the earlier
  continuous-flow navigation change.
- All JS-critical IDs (`#db-panel`, `#db-toggle-btn`, `#header-db-dot`, `#db-host` …), so the
  Database panel, connection check, and job polling keep working.

---

## Verified

- Booted `python app.py`; HTTP 200 on `/`, `/tree`, `/profiling`, `/utilities`,
  `/presence`, `/highres`; `/api/recent-runs` returns valid JSON.
- Recent-runs logic unit-tested with a synthetic job registry: a run shows as
  **building** while running, flips to **done** with a real leaf count on completion
  (single tree → 4 leaves; high-res → 650+1111 = 1761 combined), and still appears after
  a simulated restart (loaded from `.run_history.json`). Test data cleaned up.
- `/api/db-info` returns data whenever MySQL is reachable (same as the original).

---

## Reverting / comparing

Standalone copy — delete the whole `uniprot-lab-manager-green/` folder to discard it.
Your original `uniprot-lab-manager_copy/` is untouched.
