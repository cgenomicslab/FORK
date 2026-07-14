"""
Interactive ETE4 species tree coloured by comparison group.

Used by FORK's Comparative ("concept check") tab: given a list of taxids and a
taxid→colour map (group A vs group B vs both), build the NCBI topology and open
it in the ETE4 smartview explorer with branches coloured by group.

Modelled on ete_profile.py's smartview usage so it matches the installed ETE4
version's API.

    python ete_species_tree.py -t taxids.txt -c colormap.txt -p 5001

    taxids.txt   one taxid per line
    colormap.txt  "<taxid>\\t<hexcolour>" per line
"""

import argparse
from ete4 import NCBITaxa
from ete4.smartview import Layout, TextFace

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxids", required=True, help="one taxid per line")
parser.add_argument("-c", "--colormap", required=True, help="taxid<tab>hexcolour per line")
parser.add_argument("-p", "--port", type=int, default=5001)
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
tree = ncbi.get_topology(taxids)
tree.annotate_ncbi_taxa(taxid_attr="taxid")

# Pre-assign a branch colour to every node: a leaf takes its group's colour;
# an internal node takes the shared colour only if all its leaves agree.
for _node in tree.traverse():
    if _node.is_leaf:
        _col = taxid2color.get(str(_node.props.get("taxid", "")))
    else:
        _cols = {
            taxid2color.get(str(_l.props.get("taxid", ""))) for _l in _node.leaves()
        }
        _cols.discard(None)
        _col = next(iter(_cols)) if len(_cols) == 1 else None
    if _col:
        _node.add_prop("branch_color", _col)


def branch_col(node):
    col = node.props.get("branch_color", "black")
    return {
        "hz-line": {"stroke": col},
        "vt-line": {"stroke": col},
    }


def scientific_name_layout(node):
    if "sci_name" in node.props:
        label = f"{node.props['sci_name']} ({node.props.get('taxid', node.name)})"
    else:
        label = str(node.name)
    return TextFace(label, position="right", fs_min=6, fs_max=25)


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
    ],
    port=args.port,
    open_browser=False,
    keep_server=True,
)
