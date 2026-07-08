"""
recent_runs.py — real "Recent runs" data for the Overview dashboard (green UI).

Every analysis the app launches is registered in the in-memory job registry
(`_jobs` in app.py) with a small `meta` block (family, method, output path) and a
`created` timestamp. This module turns those jobs into dashboard rows — so a run
shows up **no matter which output path you chose**, the moment it starts.

To survive an app restart, each observed run is also mirrored into a tiny history
file next to this module (`.run_history.json`), keyed by job id. Reads/writes here
never touch your analysis outputs — only this one bookkeeping file.
"""

import os
import json
import time
from pathlib import Path

HISTORY_FILE = Path(__file__).resolve().parent / ".run_history.json"
MAX_HISTORY = 100

# job status (app.py) -> label shown on the dashboard
_STATUS_MAP = {"running": "building", "done": "done", "error": "failed"}


def _count_leaves(tree_path):
    """Leaf count of a single Newick tree = commas + 1 (cheap, exact)."""
    try:
        with open(tree_path) as fh:
            text = fh.read()
        return text.count(",") + 1 if text.strip() else None
    except Exception:
        return None


def _leaves_from_job(job):
    """Best-effort real leaf count from a finished job's result."""
    r = job.get("result") or {}

    # high-res build: sum the per-Pfam leaf counts
    summary = r.get("summary")
    if isinstance(summary, dict) and summary:
        total, got = 0, False
        for info in summary.values():
            n = (info or {}).get("n_leaves")
            if isinstance(n, int):
                total += n
                got = True
        return total if got else None

    # single tree (d3 viewer stores the tree file path)
    tp = r.get("tree_path")
    if tp and os.path.isfile(tp):
        return _count_leaves(tp)

    return None


def _job_entry(job_id, job):
    meta = job.get("meta") or {}
    family = meta.get("pfam") or job.get("type", "run")
    return {
        "run_id": job_id,
        "family": family,
        "method": meta.get("method", "—"),
        "leaves": _leaves_from_job(job),
        "mtime": job.get("created") or time.time(),
        "status": _STATUS_MAP.get(job.get("status"), job.get("status") or "—"),
        "output": meta.get("output", ""),
        "type": job.get("type", ""),  # 'highres_trees' or 'tree' — lets the dashboard link back
    }


def _load_history():
    try:
        with open(HISTORY_FILE) as fh:
            data = json.load(fh)
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _save_history(hist):
    """Keep the most recent MAX_HISTORY runs; write atomically."""
    items = sorted(hist.values(), key=lambda e: e.get("mtime", 0), reverse=True)
    hist = {e["run_id"]: e for e in items[:MAX_HISTORY]}
    try:
        tmp = str(HISTORY_FILE) + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(hist, fh)
        os.replace(tmp, HISTORY_FILE)
    except Exception:
        pass
    return hist


def collect(jobs=None, limit=12):
    """
    Merge live jobs with the persisted history, newest first.

    A live job always wins over its stored copy, so a run flips from
    'building' to 'done' (and gains a real leaf count) as it progresses.
    """
    hist = _load_history()
    if jobs:
        for jid, job in list(jobs.items()):
            hist[jid] = _job_entry(jid, job)
    hist = _save_history(hist)

    entries = sorted(hist.values(), key=lambda e: e.get("mtime", 0), reverse=True)
    return entries[:limit]
