"""ETE4 interactive viewer for a High-Res Profile matrix on the NCBI taxonomy tree.

Reads the matrix CSV produced by the 'Assemble Profile' step, extracts taxids
from row labels, builds the NCBI species tree, and renders each matrix column
(subclade) as a coloured aligned rectangle — matching the style of ete_profile.py.

Called by Flask via subprocess:
    python ete_highres_profile.py -m matrix.csv [-c colormap.txt] [-p PORT]
"""

import argparse
import csv
import sys
import seaborn as sns
from ete4 import NCBITaxa
from ete4.smartview import Layout, TextFace, RectFace
import utils

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--matrix",   help="Profile matrix CSV from High-Res Profile step 3", required=True)
parser.add_argument("-c", "--colormap", help="Branch-colour map: taxid <TAB> hex-colour, one per line")
parser.add_argument("-p", "--port",     help="Port for ete4 explore server", type=int, default=5003)
args = parser.parse_args()

# ── Load matrix ─────────────────────────────────────────────────────────────

matrix  = {}   # {taxid_str: {col: float}}
columns = []
taxid_to_label = {}

with open(args.matrix, newline='') as fh:
    reader = csv.DictReader(fh)
    fieldnames = reader.fieldnames or []
    if not fieldnames:
        print("Empty or invalid matrix CSV.")
        sys.exit(1)
    row_col   = fieldnames[0]
    columns   = fieldnames[1:]

    for row in reader:
        label  = row[row_col]
        # Row labels are "taxid  sci_name" (two spaces) — extract taxid
        taxid  = label.strip().split()[0]
        if not taxid.isdigit():
            continue
        matrix[taxid] = {}
        for col in columns:
            try:
                matrix[taxid][col] = float(row[col])
            except (ValueError, TypeError):
                matrix[taxid][col] = 0.0
        taxid_to_label[taxid] = label.strip()

if len(matrix) < 2:
    print(f"Need at least 2 taxids in the matrix (found {len(matrix)}).")
    sys.exit(1)

taxids = list(matrix.keys())

# ── Optional colormap ────────────────────────────────────────────────────────

taxid2color: dict[str, str] = {}
if args.colormap:
    with open(args.colormap) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                taxid2color[parts[0]] = parts[1]

# ── NCBI taxonomy tree ───────────────────────────────────────────────────────

ncbi = NCBITaxa()
tree = ncbi.get_topology(taxids, intermediate_nodes=False)
tree.to_ultrametric()
tax2names, tax2lineages, tax2rank = tree.annotate_ncbi_taxa(taxid_attr='taxid')

# ── Colour palette per column ────────────────────────────────────────────────

col_palette = sns.color_palette("viridis", len(columns))
col_colors  = dict(zip(columns, col_palette))
col_max     = {
    col: max((matrix.get(t, {}).get(col, 0) for t in taxids), default=1) or 1
    for col in columns
}

# ── Tree style ───────────────────────────────────────────────────────────────

my_tree_style = {
    'shape': 'rectangular',
    'aligned-leaves': True,
    'show-popup-props': None,
    'hz-line': {'stroke-width': 2, 'stroke': '#333333'},
    'vt-line': {'stroke-width': 2, 'stroke': '#333333'},
}


def branch_col(node):
    taxid = str(node.props.get('taxid', ''))
    if taxid in taxid2color:
        c = taxid2color[taxid]
        return {'hz-line': {'stroke': c}, 'vt-line': {'stroke': c}}


def name_layout(node):
    if 'sci_name' in node.props:
        sci = str(node.props['sci_name'])
        return TextFace(sci, position='right', fs_min=6, fs_max=25)
    else:
        return TextFace(node.name, position='right', fs_min=6, fs_max=25)


def header_layout(tree):
    faces = []
    for n, col in enumerate(columns):
        short = col[:22] + ('…' if len(col) > 22 else '')
        faces.append(TextFace(short, fs_min=5, fs_max=15, position='header', column=n))
    return faces


def count_layout(node):
    if not node.is_leaf:
        return
    taxid = str(node.props.get('taxid', ''))
    faces = []
    for n, col in enumerate(columns):
        val  = matrix.get(taxid, {}).get(col, 0)
        cmax = col_max[col]
        if val > 0:
            ratio  = min(val / cmax, 1.0)
            colour = utils.color_gradient('white', col_colors[col], mix=ratio)
            text   = str(int(val)) if val == int(val) else f'{val:.2f}'
        else:
            colour = '#f0f0f0'
            text   = '–'
        faces.append(RectFace(
            wmax=160, hmax=80,
            style={'fill': colour},
            position='aligned', column=n, text=text,
        ))
    return faces


# ── Launch ───────────────────────────────────────────────────────────────────

base_layout   = Layout(name='base',             draw_tree=my_tree_style)
branch_layout = Layout(name='branch colours',   draw_node=branch_col)
label_layout  = Layout(name='scientific names', draw_node=name_layout)
hdr_layout    = Layout(name='column headers',   draw_tree=header_layout)
cnt_layout    = Layout(name='subclade counts',  draw_node=count_layout)

print(f"Launching ETE profile viewer on port {args.port} with {len(taxids)} taxa, {len(columns)} columns.")

tree.explore(
    layouts=[base_layout, branch_layout, label_layout, hdr_layout, cnt_layout],
    port=args.port, open_browser=False, keep_server=True,
)
