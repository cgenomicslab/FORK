#!/usr/bin/env bash
# =============================================================================
# update_db.sh — safely import a new UniProt Reference Proteomes version into
# the FORK database and publish it to the live app, WITHOUT disrupting the
# running tool.
#
# WHY THIS IS SAFE
#   FORK's tables are versioned. A new release is imported as a NEW `version`
#   value — it only ADDS rows, never touches the version the live app is
#   serving. The app keeps querying the currently-published version the whole
#   time. Only after the import is verified do we "publish" by flipping ONE
#   line in .env (FORK_DEFAULT_VERSION) and restarting the app. If ANY step
#   before publish fails, the old version stays published and users see no
#   change at all.
#
# WHAT IT DOES (in order)
#   0. Takes a single lock so two updates can never overlap.
#   1. Pre-flight: checks free disk, checks the version isn't already imported.
#   2. Import proteins        -> uniprot_sync_v7.py  -version <new>   (additive)
#   3. Populate HMM results   -> pyhmmer_hmmsearch.py --version <new> (resumable)
#   4. Verify: sane row counts for <new>; abort (no publish) if not.
#   5. Publish: set FORK_DEFAULT_VERSION=<new> in .env; restart the app.
#   6. Health check: confirm the app is up and serving.
#   7. Prune: if more than KEEP_VERSIONS exist, delete the oldest version's
#      rows from every version-keyed table (batched, so it never locks the DB).
#   Everything is logged; heavy steps run under nice/ionice to yield to the app.
#
# USAGE
#   ./update_db.sh <version>              # e.g. ./update_db.sh 2026_07
#   ./update_db.sh <version> --dry-run    # print what it would do, change nothing
#   ./update_db.sh <version> --no-prune   # import + publish, keep all versions
#   ./update_db.sh <version> --gc-sequences   # also GC orphan sequences (slow)
#
# Run it from an admin shell (or a systemd timer). It is NOT triggered by the
# web app and does not need the app stopped.
# =============================================================================

set -euo pipefail

# ----------------------------- configuration --------------------------------
FORK_DIR="${FORK_DIR:-/home/cglab/FORK}"          # where app.py + .env live
SETUP_DIR="${SETUP_DIR:-$FORK_DIR/setup}"         # where the sync scripts live
ENV_FILE="${ENV_FILE:-$FORK_DIR/.env}"
CONDA_ENV="${CONDA_ENV:-bio_tools}"
CONDA_SH="${CONDA_SH:-/opt/conda/etc/profile.d/conda.sh}"
FORK_SERVICE="${FORK_SERVICE:-fork}"              # systemd unit name
APP_URL="${APP_URL:-http://127.0.0.1:8080/}"      # local health-check URL

KEEP_VERSIONS="${KEEP_VERSIONS:-2}"               # how many versions to retain
MIN_FREE_GB="${MIN_FREE_GB:-250}"                 # abort import if less free
IMPORT_THREADS="${IMPORT_THREADS:-4}"             # cap for pyhmmer, etc.
DELETE_BATCH="${DELETE_BATCH:-50000}"             # rows per prune DELETE

# REQUIRED for the HMM step: the directory that contains hmm_profiles/Pfam-A.hmm
# (the same Pfam HMM database used when the DB was first built). Set it here or
# via the environment before running.
HMM_OUTPUT_DIR="${HMM_OUTPUT_DIR:-}"
SYNC_EXTRA_ARGS="${SYNC_EXTRA_ARGS:-}"            # extra flags for uniprot_sync_v7.py
PYHMMER_EXTRA_ARGS="${PYHMMER_EXTRA_ARGS:-}"      # extra flags for pyhmmer_hmmsearch.py
REUSE_PRIOR_HMM="${REUSE_PRIOR_HMM:-0}"           # 1 = copy HMM results for sequences
                                                  # unchanged from a prior version instead
                                                  # of re-searching (VALIDATE first!)

LOG_DIR="${LOG_DIR:-$FORK_DIR/logs}"
LOCK_FILE="/tmp/fork_update_db.lock"

# ----------------------------- args -----------------------------------------
NEW="${1:-}"
DRY_RUN=0; NO_PRUNE=0; GC_SEQ=0
shift || true
for arg in "$@"; do
  case "$arg" in
    --dry-run)      DRY_RUN=1 ;;
    --no-prune)     NO_PRUNE=1 ;;
    --gc-sequences) GC_SEQ=1 ;;
    *) echo "Unknown option: $arg" >&2; exit 2 ;;
  esac
done
if [[ -z "$NEW" || ! "$NEW" =~ ^[0-9]{4}_[0-9]{2}$ ]]; then
  echo "Usage: $0 <version>  (e.g. 2026_07)  [--dry-run] [--no-prune] [--gc-sequences]" >&2
  exit 2
fi

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/update_db_${NEW}_$(date +%Y%m%d_%H%M%S).log"

log()  { echo "[$(date '+%F %T')] $*" | tee -a "$LOG_FILE"; }
die()  { log "ERROR: $*"; log "ABORTED — the published version was NOT changed; the live tool is unaffected."; exit 1; }
run()  { if [[ $DRY_RUN -eq 1 ]]; then log "DRY-RUN would run: $*"; else log "+ $*"; eval "$@"; fi; }

# ----------------------------- DB helper ------------------------------------
# Read DB creds from .env and use a private defaults-file so the password never
# appears in `ps`/the process list.
[[ -f "$ENV_FILE" ]] || die ".env not found at $ENV_FILE"
# shellcheck disable=SC1090
DB_USER="$(grep -E '^DB_USER=' "$ENV_FILE" | cut -d= -f2-)"
DB_PASS="$(grep -E '^DB_PASSWORD=' "$ENV_FILE" | cut -d= -f2-)"
DB_NAME="$(grep -E '^DB_NAME=' "$ENV_FILE" | cut -d= -f2-)"
DB_HOST="$(grep -E '^DB_HOST=' "$ENV_FILE" | cut -d= -f2- || true)"; DB_HOST="${DB_HOST:-localhost}"

# The sync scripts read DB creds from the environment (they load .env from their
# OWN setup/ folder, which may not exist). Export what we parsed from the app's
# .env so they connect no matter where they look.
export DB_HOST DB_USER DB_NAME
export DB_PASSWORD="$DB_PASS"

MYCNF="$(mktemp)"; chmod 600 "$MYCNF"
cat > "$MYCNF" <<EOF
[client]
user=$DB_USER
password=$DB_PASS
host=$DB_HOST
EOF
cleanup() { rm -f "$MYCNF"; }
trap cleanup EXIT

sql() { mysql --defaults-extra-file="$MYCNF" -N -B "$DB_NAME" -e "$1"; }

# ----------------------------- lock -----------------------------------------
exec 9>"$LOCK_FILE"
flock -n 9 || die "Another update is already running (lock $LOCK_FILE held)."

log "=== FORK DB update -> version $NEW  (dry-run=$DRY_RUN) ==="

# ----------------------------- 1. pre-flight --------------------------------
FREE_GB=$(df -BG --output=avail "$(dirname "$(readlink -f /var/lib/mysql)")" 2>/dev/null | tail -1 | tr -dc '0-9')
FREE_GB="${FREE_GB:-0}"
log "Free disk where MySQL lives: ${FREE_GB} GB (need >= ${MIN_FREE_GB})"
[[ "$FREE_GB" -ge "$MIN_FREE_GB" ]] || die "Not enough free disk to import a new version safely."

# The HMM step needs the Pfam HMM database present.
[[ -n "$HMM_OUTPUT_DIR" ]] || die "HMM_OUTPUT_DIR is not set — point it at the dir containing hmm_profiles/Pfam-A.hmm (edit the config at the top, or export it)."
[[ -f "$HMM_OUTPUT_DIR/hmm_profiles/Pfam-A.hmm" ]] || die "Pfam-A.hmm not found at $HMM_OUTPUT_DIR/hmm_profiles/Pfam-A.hmm"

EXISTING=$(sql "SELECT COUNT(*) FROM proteins WHERE version='$NEW';" || echo 0)
if [[ "${EXISTING:-0}" -gt 0 ]]; then
  log "Version $NEW already has $EXISTING protein rows — assuming already imported. Skipping import."
  SKIP_IMPORT=1
else
  SKIP_IMPORT=0
fi

# Activate the conda env for the sync scripts.
if [[ $DRY_RUN -eq 0 ]]; then
  # shellcheck disable=SC1090
  source "$CONDA_SH" && conda activate "$CONDA_ENV"
fi
NICE="nice -n 19 ionice -c3"

# ----------------------------- 2. import proteins ---------------------------
if [[ "$SKIP_IMPORT" -eq 0 ]]; then
  log "Importing proteins for $NEW (additive; existing versions untouched)…"
  run "$NICE python '$SETUP_DIR/uniprot_sync_v7.py' -version '$NEW' $SYNC_EXTRA_ARGS"
fi

# ----------------------------- 3. HMM search --------------------------------
log "Running HMM search for $NEW (resumable via its checkpoint)…"
REUSE_FLAG=""
if [[ "$REUSE_PRIOR_HMM" -eq 1 ]]; then
  REUSE_FLAG="--reuse-prior"
  log "REUSE-PRIOR enabled — only genuinely new sequences are searched; unchanged ones are copied from the prior version."
fi
run "$NICE python '$SETUP_DIR/pyhmmer_hmmsearch.py' --version '$NEW' --output-dir '$HMM_OUTPUT_DIR' $REUSE_FLAG $PYHMMER_EXTRA_ARGS"

# ----------------------------- 4. verify ------------------------------------
if [[ $DRY_RUN -eq 0 ]]; then
  N_PROT=$(sql "SELECT COUNT(*) FROM proteins WHERE version='$NEW';")
  N_HMM=$(sql  "SELECT COUNT(*) FROM hmm_search_results WHERE version='$NEW';")
  log "Verify: version $NEW has $N_PROT proteins, $N_HMM HMM hits."
  [[ "${N_PROT:-0}" -gt 100000 ]] || die "Suspiciously few proteins for $NEW ($N_PROT) — refusing to publish."
  [[ "${N_HMM:-0}"  -gt 100000 ]] || die "Suspiciously few HMM hits for $NEW ($N_HMM) — refusing to publish."
else
  log "DRY-RUN: skipping verify."
fi

# ----------------------------- 5. publish -----------------------------------
log "Publishing: setting FORK_DEFAULT_VERSION=$NEW in $ENV_FILE and restarting '$FORK_SERVICE'…"
if [[ $DRY_RUN -eq 0 ]]; then
  cp "$ENV_FILE" "$ENV_FILE.bak.$(date +%s)"
  if grep -q '^FORK_DEFAULT_VERSION=' "$ENV_FILE"; then
    sed -i "s/^FORK_DEFAULT_VERSION=.*/FORK_DEFAULT_VERSION=$NEW/" "$ENV_FILE"
  else
    echo "FORK_DEFAULT_VERSION=$NEW" >> "$ENV_FILE"
  fi
  sudo systemctl restart "$FORK_SERVICE"
else
  log "DRY-RUN: would flip FORK_DEFAULT_VERSION and restart the service."
fi

# ----------------------------- 6. health check ------------------------------
if [[ $DRY_RUN -eq 0 ]]; then
  sleep 5
  if curl -fsS "$APP_URL" >/dev/null 2>&1; then
    log "Health check OK — app is up on $APP_URL."
  else
    die "App did not come back up after restart! Check: journalctl -u $FORK_SERVICE -n 50"
  fi
fi

# ----------------------------- 7. prune old versions ------------------------
if [[ "$NO_PRUNE" -eq 0 ]]; then
  mapfile -t VERSIONS < <(sql "SELECT DISTINCT version FROM proteins ORDER BY version DESC;")
  log "Versions present: ${VERSIONS[*]:-none}"
  if [[ "${#VERSIONS[@]}" -gt "$KEEP_VERSIONS" ]]; then
    # tables that actually have a `version` column (schema-agnostic)
    mapfile -t VTABLES < <(sql "SELECT table_name FROM information_schema.columns WHERE table_schema='$DB_NAME' AND column_name='version';")
    for (( idx=KEEP_VERSIONS; idx<${#VERSIONS[@]}; idx++ )); do
      OLD="${VERSIONS[$idx]}"
      log "Pruning old version $OLD from: ${VTABLES[*]}"
      for tbl in "${VTABLES[@]}"; do
        while :; do
          if [[ $DRY_RUN -eq 1 ]]; then log "DRY-RUN would batch-delete from $tbl WHERE version=$OLD"; break; fi
          DEL=$(sql "DELETE FROM \`$tbl\` WHERE version='$OLD' LIMIT $DELETE_BATCH; SELECT ROW_COUNT();" | tail -1)
          [[ "${DEL:-0}" -eq 0 ]] && break
          sleep 0.2   # let the live app breathe between batches
        done
      done
    done
    # Optional: reclaim sequences that no remaining version references.
    if [[ "$GC_SEQ" -eq 1 && $DRY_RUN -eq 0 ]]; then
      log "GC: deleting orphan sequences (no longer referenced by any version)…"
      sql "DELETE s FROM sequences s LEFT JOIN proteins p ON p.seq_id=s.seq_id WHERE p.seq_id IS NULL;"
    fi
  else
    log "Only ${#VERSIONS[@]} version(s) present (<= KEEP_VERSIONS=$KEEP_VERSIONS) — nothing to prune."
  fi
else
  log "Pruning disabled (--no-prune)."
fi

log "=== DONE — version $NEW imported, verified, and published. ==="
