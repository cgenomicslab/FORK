"""
NCBI species tree coloured by comparison group — used by FORK's Comparative
("concept check") tab.

Two modes, both run as a subprocess by app.py so the Qt/ETE machinery never
touches the Flask request threads (Qt rendering off the main thread can crash
the whole process):

    # interactive ETE4 smartview explorer (default)
    python ete_species_tree.py -t taxids.txt -c colormap.txt -p 5001

    # static PNG render (Qt treeview) -> outfile
    python ete_species_tree.py -t taxids.txt -c colormap.txt --render out.png

    taxids.txt    one taxid per line
    colormap.txt  "<taxid>\\t<hexcolour>" per line
"""

import argparse
import sys
from ete4 import NCBITaxa

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxids", required=True, help="one taxid per line")
parser.add_argument("-c", "--colormap", required=True, help="taxid<tab>hexcolour per line")
parser.add_argument("-p", "--port", type=int, default=5001, help="port (interactive mode)")
parser.add_argument("--render", help="render a static PNG to this path and exit")
args = parser.parse_args()

# ---- inputs -------------------------------------------------------------
taxids = []
for line in open(args.taxids):
    line = line.strip()
    if line and not line.startswith("#"):
        taxids.append(line.split()[0])

taxid2color = {}
for line in open(args.colormap):
    parts = line.split()
    if len(parts) >= 2:
        taxid2color[parts[0]] = parts[1]

# ---- NCBI species tree --------------------------------------------------
ncbi = NCBITaxa()
tree = ncbi.get_topology(taxids, intermediate_nodes=False)
tree.annotate_ncbi_taxa(taxid_attr="taxid")


def _node_color(node):
    """Group colour for a node: a leaf takes its own; an internal node takes
    the shared colour only when all its leaves agree."""
    if node.is_leaf:
        return taxid2color.get(str(node.props.get("taxid", "")))
    cols = {taxid2color.get(str(l.props.get("taxid", ""))) for l in node.leaves()}
    cols.discard(None)
    return next(iter(cols)) if len(cols) == 1 else None


# Pre-assign the colour as a prop so both renderers can read it.
for _node in tree.traverse():
    _col = _node_color(_node)
    if _col:
        _node.add_prop("branch_color", _col)


# =========================================================================
# STATIC PNG (Qt treeview) — runs in this subprocess's main thread
# =========================================================================
if args.render:
    from ete4.treeview import TreeStyle, NodeStyle, TextFace as TVTextFace

    for _n in tree.traverse():
        _n.dist = 0.0
    tree.to_ultrametric(topological=True)

    for node in tree.traverse():
        ns = NodeStyle()
        ns["hz_line_width"] = 4
        ns["vt_line_width"] = 4
        ns["size"] = 0
        col = node.props.get("branch_color")
        if col:
            ns["hz_line_color"] = col
            ns["vt_line_color"] = col
        node.set_style(ns)
        if node.is_leaf:
            label = node.props.get("sci_name", node.name)
            node.add_face(
                TVTextFace(f"  {label} ({node.props.get('taxid', node.name)})"),
                column=0,
                position="branch-right",
            )

    ts = TreeStyle()
    ts.show_leaf_name = False
    tree.render(args.render, w=2000, units="px", tree_style=ts)
    sys.exit(0)


# =========================================================================
# INTERACTIVE ETE4 smartview explorer
# =========================================================================
from ete4.smartview import Layout, TextFace


def branch_col(node):
    col = node.props.get("branch_color", "black")
    return {
        "hz-line": {"stroke": col},
        "vt-line": {"stroke": col},
    }


def scientific_name_layout(node):
    if not node.is_leaf:
        return
    if "sci_name" in node.props:
        label = f"{node.props['sci_name']} ({node.props.get('taxid', node.name)})"
    else:
        label = str(node.name)
    return TextFace(label, position="right", fs_min=6, fs_max=25)


def clade_name_layout(node):
    if node.is_leaf or "sci_name" not in node.props:
        return
    return TextFace(str(node.props["sci_name"]), position="right", fs_min=6, fs_max=25)


my_tree_style = {
    "shape": "rectangular",
    "aligned-leaves": True,
    "show-popup-props": None,
}

# Ensure a drawable width (topology-only NCBI trees are zero-length) — same
# guard used elsewhere so SmartView never raises "Cannot draw tree with width 0".
for _n in tree.traverse():
    if _n.up is not None and (not _n.dist or _n.dist <= 0):
        _n.dist = 1e-6


def _root_dist(_leaf):
    _d, _cur = 0.0, _leaf
    while _cur.up is not None:
        _d += _cur.dist or 0.0
        _cur = _cur.up
    return _d


if max((_root_dist(_l) for _l in tree.leaves()), default=0.0) < 1e-5:
    for _n in tree.traverse():
        if _n.up is not None:
            _n.dist = 1.0

tree.explore(
    layouts=[
        Layout(name="example", draw_tree=my_tree_style),
        Layout(name="group colours", draw_node=branch_col),
        Layout(name="scientific names", draw_node=scientific_name_layout),
        Layout(name="clade names", active=False, draw_node=clade_name_layout),
    ],
    port=args.port,
    open_browser=False,
    keep_server=True,
)
