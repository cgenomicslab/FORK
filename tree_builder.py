"""
Per-Pfam gene tree building by calling tree_from_db.py as a
subprocess, with caching and parsing of the resulting
Newick into ete4 Tree objects.

The high-resolution phylogenetic profiling pipeline:

    build_trees(pfams, ...)
            │
            ▼  {pfam: {"tree": ete4.Tree, "leaves": [...], ...}}
    subclade_partition.partition_by_depth/by_mrca   per tree
            │
            ▼  {pfam: {"A": {accs}, "B": {accs}, ...}}
    library.fetch_highres_profile(...)
            │
            ▼  pandas matrix

Cache
----------------
Each (pfam, taxids, exclude_taxids, version, evalue, aln, ml, gt) combination
gets its own subdirectory under `output_root`, named after a short MD5 hash
of the parameters. Inside each subdirectory, tree_from_db.py is invoked
with --prefix={subdir}/{pfam}. Re-running with identical parameters reuses
the existing tree (tree_from_db.py itself checks for file existence at
every stage). Re-running with different parameters builds in a new
subdirectory — old runs are not deleted.

"""

import os
import subprocess
import hashlib
import json
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional

from ete4 import Tree

# -----------------------------------------------------------------------------
# Leaf-name
# -----------------------------------------------------------------------------


"""
tree_from_db.py writes FASTA headers as ">{taxon_id}.{accession}", so tree
leaves look like "9606.P04637". `parse_leaf_to_accession` and
`strip_leaf_prefix_in_subclades` convert these back to bare accessions
for the profile-assembly step.
"""


def parse_leaf_to_accession(leaf_name: str) -> str:
    """
    "9606.P04637"  ->  "P04637"

    Splits on the FIRST dot.
    """
    parts = leaf_name.split(".", 1)
    return parts[1] if len(parts) > 1 else leaf_name


def parse_leaf_to_taxid(leaf_name: str) -> Optional[int]:
    """ "9606.P04637" -> 9606 ; returns None if not parseable."""
    parts = leaf_name.split(".", 1)
    if len(parts) > 1:
        try:
            return int(parts[0])
        except ValueError:
            return None
    return None


def strip_leaf_prefix_in_subclades(
    subclades: Dict[str, set],
) -> Dict[str, set]:
    """
    Convert a {label: set(leaf_names)} dict (output of subclade_partition)
    into {label: set(accessions)} by stripping the "taxid." prefix from
    each leaf.
    """
    return {
        label: {parse_leaf_to_accession(leaf) for leaf in leaves}
        for label, leaves in subclades.items()
    }


# -----------------------------------------------------------------------------
# Cache key + path construction
# -----------------------------------------------------------------------------


def cache_key(
    pfam: str,
    taxids: Optional[List[int]],
    exclude_taxids: Optional[List[int]],
    version: str,
    evalue: Optional[float],
    aln: str,
    ml: str,
    gt: str,
    local_fasta_hash: Optional[str] = None,
) -> str:
    """
    Short ID for one tree-build configuration. Identical
    inputs always produce the same key - changing any parameter produces
    a different key.

    For entries built from an uploaded FASTA (rather than the DB), a hash of
    the file content is folded in so that re-uploading a *different* FASTA
    under the same label rebuilds instead of reusing a stale cached tree.
    """
    params = {
        "pfam": pfam,
        "tax": sorted(taxids) if taxids else None,
        "extax": sorted(exclude_taxids) if exclude_taxids else None,
        "ver": version,
        "evalue": evalue,
        "aln": aln,
        "ml": ml,
        "gt": gt,
        "local_fasta": local_fasta_hash,
    }
    blob = json.dumps(params, sort_keys=True).encode()
    return f"{pfam}_{hashlib.md5(blob).hexdigest()[:10]}"


def _tree_filename(prefix: str, aln: str, ml: str, gt: str) -> str:
    """
    Reproduce tree_from_db.py's output naming so we can locate the tree
    file without scanning the directory.

    {prefix}.{aln_ext}  ->  +.gt{gt_clean}  ->  +.lg.fasttree / +.treefile
    """
    aln_ext_map = {"mafft": ".mft", "einsi": ".einsi", "clustalo": ".clustalo"}
    if aln not in aln_ext_map:
        raise ValueError(f"unknown aln method: {aln}")
    aln_path = f"{prefix}{aln_ext_map[aln]}"

    gt_clean = gt.replace(
        ".", ""
    )  # matches tree_from_db.py args['gt'].replace(".", "")
    trimal_path = f"{aln_path}.gt{gt_clean}"

    if ml == "fasttree":
        return f"{trimal_path}.lg.fasttree"
    if ml == "iqtree":
        return f"{trimal_path}.treefile"
    raise ValueError(f"unknown ml method: {ml}")


# -----------------------------------------------------------------------------
# Building tree
# -----------------------------------------------------------------------------

ProgressCb = Optional[Callable[[str], None]]


def build_one_tree(
    pfam: str,
    output_root: str,
    version: str,
    taxids: Optional[List[int]] = None,
    exclude_taxids: Optional[List[int]] = None,
    evalue: Optional[float] = None,
    aln: str = "mafft",
    ml: str = "fasttree",
    gt: str = "0.01",
    cpu: int = 4,
    tree_from_db_path: str = "tree_from_db.py",
    python_exe: Optional[str] = None,
    force: bool = False,
    progress_callback: ProgressCb = None,
    local_fasta: Optional[str] = None,
) -> dict:
    """
    Build (or load from cache) the tree for ONE entry.

    Normally the sequences are pulled from the reference DB by `pfam`. If
    `local_fasta` is given, that file's sequences are used instead (headers
    must be "{taxid}.{accession}") and `pfam` is treated as a display label.
    """
    local_fasta_hash = None
    if local_fasta:
        with open(local_fasta, "rb") as fh:
            local_fasta_hash = hashlib.md5(fh.read()).hexdigest()
    key = cache_key(
        pfam, taxids, exclude_taxids, version, evalue, aln, ml, gt, local_fasta_hash
    )
    cache_dir = Path(output_root) / key
    cache_dir.mkdir(parents=True, exist_ok=True)
    prefix = str(cache_dir / pfam)

    tree_path = _tree_filename(prefix, aln, ml, gt)
    cached = os.path.isfile(tree_path) and not force

    result = {
        "pfam": pfam,
        "cache_key": key,
        "prefix": prefix,
        "tree_path": tree_path,
        "resolved_tree_path": None,
        "tree": None,
        "leaves": [],
        "cached": cached,
        "stderr": None,
        "error": None,
    }

    if not cached:
        if progress_callback:
            progress_callback(f"[{pfam}] Aligning sequences and inferring tree…")
        cmd = [
            python_exe or sys.executable,
            tree_from_db_path,
            "--pfam",
            pfam,
            "--version",
            version,
            "--prefix",
            prefix,
            "--aln",
            aln,
            "--ml",
            ml,
            "--gt",
            gt,
            "--cpu",
            str(cpu),
            "--no_ncbi",  # skip NCBI annotation — irrelevant for profiling
            "--no_explore",  # no ETE4 server
        ]
        if taxids:
            cmd += ["--taxids", ",".join(str(t) for t in taxids)]
        if exclude_taxids:
            cmd += ["--exclude_taxids", ",".join(str(t) for t in exclude_taxids)]
        if evalue is not None:
            cmd += ["--evalue", str(evalue)]
        if local_fasta:
            cmd += ["--local_fasta", local_fasta]

        proc = subprocess.run(cmd, capture_output=True, text=True)
        result["stderr"] = proc.stderr

        if proc.returncode != 0:
            result["error"] = (
                f"Tree build failed (exit code {proc.returncode})."
            )
            return result
    else:
        if progress_callback:
            progress_callback(f"[{pfam}] Reusing previously built tree…")

    if not os.path.isfile(tree_path):
        result["error"] = f"expected tree file missing: {tree_path}"
        return result

    try:
        with open(tree_path) as f:
            nwk = f.read().strip()
        tree = Tree(nwk)
        # FastTree writes an unrooted tree with a trifurcating root
        # (children 0,1,2). ETE's explorer displays it as binary. Resolve
        # the polytomy here so the node-path numbering the user reads in the
        # ETE viewer matches what partition_by_node_path resolves against.
        if len(tree.children) > 2:
            tree.resolve_polytomy()
        tree.ladderize()
        # Write the resolved tree so the ETE explorer and the node-path
        # partitioner read identical structure (identical 0/1 numbering).
        resolved_path = tree_path + ".resolved.nwk"
        tree.write(outfile=resolved_path)
        result["resolved_tree_path"] = resolved_path
    except Exception as e:
        result["error"] = f"failed to parse Newick: {e}"
        return result

    result["tree"] = tree
    result["leaves"] = [leaf.name for leaf in tree.leaves()]
    return result


# -----------------------------------------------------------------------------
# Building many trees
# -----------------------------------------------------------------------------


def build_trees(
    pfams: List[str],
    output_root: str,
    version: str,
    progress_callback: ProgressCb = None,
    local_fastas: Optional[Dict[str, str]] = None,
    **kwargs,
) -> Dict[str, dict]:
    """
    Build trees for multiple entries (serial). Returns {label: result_dict}.

    `pfams` are built from the reference DB. `local_fastas` is an optional
    {label: fasta_path} mapping — those entries are built from the uploaded
    FASTA instead, and appear in the results alongside the DB Pfams so the
    downstream partition/profile steps treat both uniformly.
    """
    results = {}
    entries = [(pfam, None) for pfam in pfams]
    entries += [(label, path) for label, path in (local_fastas or {}).items()]
    n = len(entries)
    for i, (label, lf) in enumerate(entries, 1):
        if progress_callback:
            src = "uploaded FASTA" if lf else "database"
            progress_callback(f"Gene tree {i} of {n} — {label} (from {src})")
        results[label] = build_one_tree(
            pfam=label,
            output_root=output_root,
            version=version,
            progress_callback=progress_callback,
            local_fasta=lf,
            **kwargs,
        )
        if progress_callback:
            r = results[label]
            ok = r["error"] is None
            status = "done" if ok else f"failed ({r['error']})"
            tag = " · reused from cache" if r["cached"] and ok else ""
            progress_callback(f"Gene tree {i} of {n} — {label}: {status}{tag}")
    return results
