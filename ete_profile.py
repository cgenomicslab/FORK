import argparse
import sys
from collections import defaultdict
import seaborn as sns
from ete4 import NCBITaxa
from ete4.smartview import Layout, TextFace, RectFace
import utils

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--input", help="tblout output file from running hmmsearch", required=True
)
parser.add_argument(
    "-c", "--colormap", help="txt file w. colours for taxonomic groups", required=True
)
parser.add_argument(
    "-t",
    "--taxids",
    help="txt file w. taxids to be visualized in species tree (one per line)",
    required=True,
)
parser.add_argument("-m", "--max", help="max value of domain counts for visualization")
parser.add_argument(
    "-p", "--port", help="port selection for ete explore", type=int, default=5001
)
args = parser.parse_args()

if args.taxids:
    speciesfile = args.taxids
    taxids = []
    with open(speciesfile) as f:
        lines = f.readlines()
        for line in lines:
            taxid = line.strip()
            taxids.append(taxid)
else:
    print("No valid taxid file inputted (requires txt file with one taxid per line).")
    sys.exit()

if args.input:
    tbloutfile = args.input

    # dictionary to store the structure: {domain: {taxid: count}}
    domain_taxid_counts = defaultdict(lambda: defaultdict(int))

    # dictionary to store the structure: {domain: {seqid: count}}
    domain_seqid_counts = defaultdict(lambda: defaultdict(int))

    domains2screen = []

    with open(tbloutfile) as file:

        lines = file.readlines()

        for line in lines:

            if line.startswith("#"):
                continue

            fields = line.split()

            # ensure there are enough columns
            if len(fields) < 8:
                print("Not valid format of tblout file!")
                sys.exit()

            # domain E-value
            evalue = float(fields[7])

            # keep only significant hits
            if evalue >= 1e-5:
                continue

            domain = fields[3].split(".")[0]
            seqid = fields[0]
            taxid = seqid.split(".")[0]

            domain_taxid_counts[domain][taxid] += 1
            domain_seqid_counts[domain][seqid] += 1

            if domain not in domains2screen:
                domains2screen.append(domain)

else:
    print("No valid HMMsearch file inputted " "(requires tblout output file).")
    sys.exit()

# custom branch colors based on taxonomic groups
if args.colormap:
    txtfile = args.colormap
    taxid2color = {}
    with open(txtfile) as f:
        lines = f.readlines()
        for line in lines:
            if not line.startswith("#"):
                cols = line.strip().split("\t")
                if len(cols) >= 2:
                    taxon_group = cols[0]
                    hex_color = cols[1]
                    taxid2color[taxon_group] = hex_color
                else:
                    print("Not valid format of colormap file!")
                    sys.exit()
else:
    print(
        "Requires colormap to distinguish taxonomic groups (txt file, tab separated, group: colour)."
    )
    sys.exit()


if args.max:
    max_dom_value = int(args.max)
else:
    max_dom_value = 50
print(f"max colour capped on domain count: {max_dom_value}")

# NCBI taxonomy
ncbi = NCBITaxa()

# NCBI species tree
tree = ncbi.get_topology(taxids)
print(tree)
tree.to_ultrametric()

# add the rest of the props w. NCBI db (sci_name, taxid, named_lineage, lineage, rank)
tax2names, tax2lineages, tax2rank = tree.annotate_ncbi_taxa(taxid_attr="taxid")

# generate a list of colors for each domain
custom_palette = sns.color_palette("viridis", len(domains2screen))

# define domain colors for the gradient
domain_colors = {}
for domain in domains2screen:
    domain_colors[domain] = custom_palette.pop(0)

# generalized style of the tree
my_tree_style = {
    "shape": "rectangular",
    "aligned-leaves": True,
    "show-popup-props": None,
    "hz-line": {"stroke-width": 2, "stroke": "black"},
    "vt-line": {"stroke-width": 2, "stroke": "black"},
}

# modify branch colour based on taxonomic grouping
priority = [
    "9606",  # Human
    "936053",  # Rhizopus delemar
    "4827",  # Mucorales
    "4751",  # Fungi
    "33090",  # Plants
    "33213",  # Bilateria
]


def branch_col(node):

    tax_id = int(node.props["taxid"])
    lineage = tax2lineages[tax_id]

    for taxid in priority:

        if int(taxid) in lineage:

            return {
                "hz-line": {"stroke": taxid2color[taxid]},
                "vt-line": {"stroke": taxid2color[taxid]},
            }


# display their scientific names
def scientific_name_layout(node):
    if "sci_name" in node.props:
        scient_name = str(node.props["sci_name"])
        return TextFace(scient_name, position="right", fs_min=6, fs_max=25)
    else:
        return TextFace(node.name, position="right", fs_min=6, fs_max=25)


# domain visualization headers
def default_header(tree):
    faces = []
    for n, domain in enumerate(domains2screen):
        title = f"{domain}"
        faces.append(TextFace(title, fs_min=6, fs_max=25, position="header", column=n))
    return faces


# domain count visualization (per leaf node)
def default_counts(node):
    maxval = max_dom_value
    faces = []
    if node.is_leaf:
        taxid = str(node.props["taxid"])
        for n, domain in enumerate(domains2screen):
            if taxid in domain_taxid_counts[domain]:
                count = int(domain_taxid_counts[domain][taxid])
            else:
                count = 0

            if count > 0:
                cap = min(count, maxval)
                ratio = cap / maxval
                colour = utils.color_gradient("white", domain_colors[domain], mix=ratio)
                text = str(count)
            else:
                text = "-"
                colour = "red"
            faces.append(
                RectFace(
                    wmax=200,
                    hmax=100,
                    style={"fill": colour},
                    position="aligned",
                    column=n,
                    text=text,
                )
            )
        return faces


# layouts
base_layout = Layout(name="example", draw_tree=my_tree_style)
node_layout = Layout(name="branch colours", draw_node=branch_col)
name_layout = Layout(name="scientific names", draw_node=scientific_name_layout)
def_header_layout = Layout(name="heatmap header", draw_tree=default_header)
def_counts_layout = Layout(name="domain counts", draw_node=default_counts)


tree.explore(
    layouts=[
        base_layout,
        node_layout,
        name_layout,
        def_header_layout,
        def_counts_layout,
    ],
    port=args.port,
    open_browser=False,
    keep_server=True,
)
