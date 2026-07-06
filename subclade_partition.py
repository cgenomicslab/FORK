"""
Partition a gene tree into subclades for high-resolution phylogenetic profiling.

Two modes are supported, both returning the same data structure:

    {subclade_label: set(leaf_names)}

where leaf_names are the strings stored in `node.name` of the ETE4 tree
(UniProt accessions).

Mode 1 — Depth slider:
    Cut the tree at a chosen root-to-node distance (d). Every node whose
    parent's distance-from-root is < d and whose own distance-from-root
    is >= d becomes the root of a subclade. Leaves already separated from the rest before d,
    end up as singleton subclades.

Mode 2 — Manual MRCA picking:
    User supplies one or more groups of leaf names. For each group we
    compute the most recent common ancestor (MRCA) and call its leaves
    one subclade. Leaves not under any picked MRCA go into 'unassigned'
    (or are dropped, depending on `include_unassigned`).

Labels are assigned A, B, C, ... in ladderized left-to-right order so
that the column order in the final profile matrix is reproducible.

Quick usage
-----------
    from ete4 import Tree
    import subclade_partition as sp

    tree = Tree(open("PF00041.nwk").read())
    # Mode 1
    parts = sp.partition_by_depth(tree, threshold=0.8)
    # Mode 2
    parts = sp.partition_by_mrca(tree, [
        ["P12345_HUMAN", "Q67890_MOUSE"],   # -> MRCA of these two is subclade 1
        ["O11111_DROME"],
    ])
"""

from typing import Dict, Iterable, List, Set, Union

from ete4 import Tree

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------


def _dist_from_root(node) -> float:
    """
    Sum of branch lengths from the tree root to `node`..
    """
    d = 0.0
    cur = node
    while cur.up is not None:
        d += float(cur.dist or 0.0)
        cur = cur.up
    return d


def get_max_root_distance(tree: Tree) -> float:
    """
    Maximum root-to-leaf distance in the tree. We use this to set the upper
    bound of a depth slider in the UI.
    """
    return max(_dist_from_root(leaf) for leaf in tree.leaves())


def _leaf_names(node) -> Set[str]:
    """All leaf names below `node`, as a set."""
    return {leaf.name for leaf in node.leaves()}


def _ladderize_label(tree: Tree, subclade_roots: List) -> Dict[str, Set[str]]:
    """
    Given a list of subclade root nodes, return a dict labeled A, B, C, ...
    in ladderized left-to-right leaf order.

    "Ladderized order" = the order leaves appear when you draw the tree with
    smaller subtrees on top. We ladderize a copy of the tree (cheap), then
    record the order in which subclade roots are encountered while walking
    the tree in preorder.
    """
    # Ladderize in place (ETE4: sorts children by subtree size)
    tree.ladderize()

    roots_set = set(id(r) for r in subclade_roots)
    ordered_roots = []
    for node in tree.traverse("preorder"):
        if id(node) in roots_set and node not in ordered_roots:
            ordered_roots.append(node)

    # If any subclade roots weren't reached by traversal, append them.
    for r in subclade_roots:
        if r not in ordered_roots:
            ordered_roots.append(r)

    return {
        _index_to_label(i): _leaf_names(root) for i, root in enumerate(ordered_roots)
    }


def _index_to_label(i: int) -> str:
    """
    0 -> 'A', 1 -> 'B', ..., 25 -> 'Z', 26 -> 'AA', 27 -> 'AB', ...
    Excel-column-style. Trees with > 26 subclades happen.
    """
    label = ""
    n = i
    while True:
        label = chr(ord("A") + (n % 26)) + label
        n = n // 26 - 1
        if n < 0:
            return label


# -----------------------------------------------------------------------------
# Mode 1 — Depth slider
# -----------------------------------------------------------------------------


def partition_by_depth(tree: Tree, threshold: float) -> Dict[str, Set[str]]:
    """
    Partition the tree by cutting at root-to-node distance == threshold.

    A subclade root is any node whose parent's distance-from-root is < threshold
    and whose own distance-from-root is >= threshold. Walking preorder, as
    soon as we cross the threshold we record that node and prune the recursion
    (we don't descend further into already-claimed subtrees).

    - If `threshold` is shallower than every leaf, every direct child of root
      becomes a subclade .
    - If `threshold` is deeper than every internal node, every leaf is its own
      subclade .
    """
    if threshold <= 0:
        raise ValueError(f"threshold must be > 0 (got {threshold})")

    subclade_roots = []
    # Preorder, but we want to skip descendants of already-claimed nodes.
    stack = [tree]
    while stack:
        node = stack.pop()
        d = _dist_from_root(node)

        if node.up is not None:
            d_parent = _dist_from_root(node.up)
        else:
            d_parent = -1.0  # root: nothing above it

        if d_parent < threshold <= d:
            # This node is a subclade root. Don't descend further.
            subclade_roots.append(node)
        else:
            # Either we haven't reached the threshold yet (descend),
            # or this node is below the threshold but its parent already was

            for child in node.children:
                stack.append(child)

    return _ladderize_label(tree, subclade_roots)


# -----------------------------------------------------------------------------
# Mode 2 — Manual MRCA picking
# -----------------------------------------------------------------------------

MRCASpec = Union[str, Iterable[str]]


def partition_by_mrca(
    tree: Tree,
    mrca_specs: List[MRCASpec],
    include_unassigned: bool = True,
) -> Dict[str, Set[str]]:
    """
    Partition the tree by user-specified MRCAs.

    Each entry in `mrca_specs` defines one subclade:
      - a single leaf name (str)       -> singleton subclade {leaf}
      - an iterable of leaf names      -> subclade = all leaves under their MRCA

    Overlapping MRCAs are handled by *first-come-first-served*: if leaf L is
    under two picked MRCAs, it's assigned to whichever subclade was specified
    first. (We warn but don't raise.)

    Leaves under no picked MRCA become subclade 'unassigned' if
    `include_unassigned=True`; otherwise they are dropped.


    """
    all_leaves = set(leaf.name for leaf in tree.leaves())
    assigned: Set[str] = set()
    result: Dict[str, Set[str]] = {}

    for i, spec in enumerate(mrca_specs):
        label = _index_to_label(i)

        if isinstance(spec, str):
            leaves_here = {spec} if spec in all_leaves else set()
        else:
            spec_list = list(spec)
            spec_in_tree = [s for s in spec_list if s in all_leaves]

            if not spec_in_tree:
                leaves_here = set()
            elif len(spec_in_tree) == 1:
                leaves_here = set(spec_in_tree)
            else:
                mrca = tree.common_ancestor(spec_in_tree)
                leaves_here = _leaf_names(mrca)

        # First-come-first-served: drop leaves already assigned
        leaves_here = leaves_here - assigned
        assigned |= leaves_here
        result[label] = leaves_here

    if include_unassigned:
        leftover = all_leaves - assigned
        if leftover:
            result["unassigned"] = leftover

    return result


# -----------------------------------------------------------------------------
#     Partition the tree by explicit ETE4 node paths.
# -----------------------------------------------------------------------------


def partition_by_node_path(tree, paths):
    """

    Each path is a list of child indices from the root: [] = root,
    [0] = first child of root, [1, 1] = root's second child's second child.

    Unlike the depth slider, that takes every clade at one depth, this
    picks only the nodes the user names. Leaves not under any picked node are
    dropped entirely.

    Invalid paths are skipped.
    """
    result = {}
    label_i = 0
    for path in paths:
        node = tree
        ok = True
        for idx in path:
            if 0 <= idx < len(node.children):
                node = node.children[idx]
            else:
                ok = False
                break
        if not ok:
            continue
        result[_index_to_label(label_i)] = {leaf.name for leaf in node.leaves()}
        label_i += 1
    return result


def list_internal_nodes(tree, max_nodes=300):
    """
    List internal nodes with their ETE4 path and leaf count, to help find
    which paths to pick for partition_by_node_path.

    Returns list of dict: {"path": "0,1", "n_leaves": 42, "sample": "..."}
    Sorted by leaf count descending, capped at max_nodes.
    """
    rows = []
    stack = [(tree, [])]
    while stack:
        node, path = stack.pop()
        if not node.is_leaf:
            leaves = [leaf.name for leaf in node.leaves()]
            rows.append(
                {
                    "path": ",".join(str(i) for i in path),
                    "n_leaves": len(leaves),
                    "sample": ", ".join(sorted(leaves)[:3])
                    + (" ..." if len(leaves) > 3 else ""),
                }
            )
        for i, child in enumerate(node.children):
            stack.append((child, path + [i]))
    rows.sort(key=lambda r: r["n_leaves"], reverse=True)
    return rows[:max_nodes]


# -----------------------------------------------------------------------------
#     Partition the tree by duplication events based on taxonomic level
# -----------------------------------------------------------------------------
def _leaf_taxid(leaf_name):
    """Extract the NCBI taxid from a leaf name. Leaf names are 'taxid.accession'"""

    return leaf_name.split(".")[0]


def _target_species_set(tree, taxon_taxid, ncbi):
    """Return the set of taxids present as leaves in 'tree' whose NCBI lineage contains
    'taxon_taxid'. This way is defined which leaves belong to the taxonomic group we are
    partitioning within.
    """
    taxon_taxid = str(taxon_taxid)
    target = set()
    lineage_cache = {}

    for leaf in tree.leaves():
        taxid = _leaf_taxid(leaf.name)
        if taxid not in lineage_cache:
            try:
                lineage_cache[taxid] = {str(t) for t in (ncbi.get_lineage(int(taxid)) or [])}
            except Exception as e:
                print(f"[partition_by_duplication] lineage lookup failed for taxid={taxid!r}: {e}")
                lineage_cache[taxid] = set()
        if taxon_taxid in lineage_cache[taxid]:
            target.add(taxid)
    return target


def partition_by_duplication(tree, taxon_taxid, ncbi):
    target_species = _target_species_set(tree, taxon_taxid, ncbi)
    if not target_species:
        return {}
    # node -> set of target species found under that node
    node_species = {}
    duplication_nodes = []

    for node in tree.traverse(
        "postorder"
    ):  # leaves first, root last -- walking the tree bottom-up
        # we process the children before parents, so by the time we reach an internal node, we know what's below each of its children
        if node.is_leaf:
            taxid = _leaf_taxid(node.name)
            node_species[node] = {taxid} if taxid in target_species else set()
            continue

        child_sets = [node_species[c] for c in node.children]
        node_species[node] = set().union(*child_sets) if child_sets else set()

        is_duplication = False
        for i in range(len(child_sets)):
            for j in range(i + 1, len(child_sets)):
                if child_sets[i] & child_sets[j]:
                    is_duplication = True
                    break
            if is_duplication:
                break
        if is_duplication:
            duplication_nodes.append(node)
    if not duplication_nodes:
        return {}

    # Only split at the outermost duplication nodes.
    # Walk preorder (root → leaves); once a subtree is claimed by a duplication
    # node, skip any nested duplication nodes inside it to avoid overlapping subclades.
    duplication_set = set(id(n) for n in duplication_nodes)
    subclade_roots = []
    claimed = set()

    for node in tree.traverse("preorder"):
        if id(node) in claimed:
            continue
        if id(node) in duplication_set:
            for child in node.children:
                subclade_roots.append(child)
                for desc in child.traverse():
                    claimed.add(id(desc))

    if not subclade_roots:
        return {}

    return _ladderize_label(tree, subclade_roots)


# -----------------------------------------------------------------------------
# To test, run: python subclade_partition.py
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    # Toy tree, branch lengths chosen so depth-cutting has obvious answers.
    #
    #          /-A    (d=1.0)
    #     /---|
    #    |    \-B    (d=1.0)
    # ---|
    #    |    /-C    (d=1.0)
    #    |---|
    #        |    /-D    (d=1.0)
    #         \--|
    #            \-E    (d=1.0)
    #
    nwk = "((A:1.0, B:1.0):0.5, (C:1.0, (D:1.0, E:1.0):0.5):0.5);"
    t = Tree(nwk)

    print("=== Tree ===")
    print(t.to_str(show_internal=False))

    print(f"\nMax root-to-leaf distance: {get_max_root_distance(t):.2f}")

    print("\n=== Mode 1: partition_by_depth(threshold=0.5) ===")
    # At d=0.5 we cut on the two branches leaving root -> 2 subclades.
    parts = partition_by_depth(Tree(nwk), threshold=0.5)
    for label, leaves in parts.items():
        print(f"  {label}: {sorted(leaves)}")

    print("\n=== Mode 1: partition_by_depth(threshold=1.0) ===")
    # Deeper cut --> more subclades.
    parts = partition_by_depth(Tree(nwk), threshold=1.0)
    for label, leaves in parts.items():
        print(f"  {label}: {sorted(leaves)}")

    print("\n=== Mode 2: partition_by_mrca([[A, B], [D, E]]) ===")
    parts = partition_by_mrca(Tree(nwk), [["A", "B"], ["D", "E"]])
    for label, leaves in parts.items():
        print(f"  {label}: {sorted(leaves)}")

    print("\n=== Mode 2: partition_by_mrca([[A, B]], include_unassigned=False) ===")
    parts = partition_by_mrca(Tree(nwk), [["A", "B"]], include_unassigned=False)
    for label, leaves in parts.items():
        print(f"  {label}: {sorted(leaves)}")
