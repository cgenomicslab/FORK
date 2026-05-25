# ------------------LIBRARIES----------------------------
from get_reference_uniprot_set_lib import (
    UniProtRetriever,
    get_db_config,
    fetch_domains_by_accession,
)
import argparse
import os
import sys
import re
import subprocess

# Safely import ETE4
try:
    from ete4 import PhyloTree
    from ete4.smartview import Layout, TextFace, SeqFace, BASIC_LAYOUT

    ETE_AVAILABLE = True
except ImportError as e:
    print(f"WARNING--Could not import ETE4 ({e}).")
    ETE_AVAILABLE = False
# --------------------------------------------------------


def get_fasta_from_db(
    version, taxids, pfams, evalue_cutoff=None, exclude_taxon_ids=None
):
    taxid_ints = [int(t) for t in taxids] if taxids else None
    seen = set()
    all_records = []

    with UniProtRetriever(get_db_config()) as db:
        for pfam in pfams:
            records = db.get_proteins_by_hmm_hit(
                version=version,
                hmm_query=pfam,
                evalue_cutoff=evalue_cutoff,
                taxon_ids=taxid_ints,
                exclude_taxon_ids=exclude_taxon_ids,
            )
            for r in records:
                if r["accession"] not in seen:
                    seen.add(r["accession"])
                    all_records.append(r)

    print("INFO--Retrieved %d unique sequences from local DB" % len(all_records))

    fasta_lines = []
    seqid2gene = {}
    for r in all_records:
        seqid = "%s.%s" % (r["taxon_id"], r["accession"])
        fasta_lines.append(">%s\n%s" % (seqid, r["sequence"]))
        seqid2gene[seqid] = r["name"]

    return "\n".join(fasta_lines) + "\n", seqid2gene


def get_seqs(fastafile):
    name2seq = {}
    seq = ""
    head = ""
    for line in open(fastafile):
        if line.startswith(">"):
            if seq:
                name2seq[head] = seq
                seq = ""
                head = line.lstrip(">").rstrip()
            else:
                head = line.lstrip(">").rstrip()
        else:
            seq += line.rstrip()
    if head:
        name2seq[head] = seq
    return name2seq


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Uniprot/Pfam-based protein family evolution analysis"
    )
    parser.add_argument("--pfam", required=True, type=str)
    parser.add_argument("--taxids", required=False, type=str)
    parser.add_argument("--cpu", type=str, default="32")
    parser.add_argument("--ml", default="fasttree", choices=["fasttree", "iqtree"])
    parser.add_argument(
        "--aln", default="mafft", choices=["mafft", "einsi", "clustalo"]
    )
    parser.add_argument("--gt", type=str, default="0.1")
    parser.add_argument("--colormap", required=False, type=str)
    parser.add_argument("--version", required=True, type=str)
    parser.add_argument("--evalue", required=False, type=float, default=None)
    parser.add_argument("--local_fasta", required=False, type=str)
    parser.add_argument("--MSA", action="store_true")
    parser.add_argument("--positions", required=False, type=str)
    parser.add_argument("--prefix", required=True, type=str)
    parser.add_argument("--port", required=False, type=int, default=5001)
    parser.add_argument("--no_ncbi", action="store_true")
    parser.add_argument("--no_explore", action="store_true")
    parser.add_argument("--exclude_taxids", required=False, type=str)

    args = vars(parser.parse_args())

    print("INFO--Processing Pfam domains:", args["pfam"])
    pfams = args["pfam"].split(",")

    if type(args["taxids"]) == str:
        if os.path.isfile(args["taxids"]):
            taxids = [
                line.strip().split()[0]
                for line in open(args["taxids"])
                if not line.startswith("#")
            ]
        else:
            try:
                taxids = args["taxids"].split(",")
                print("Taxids loaded from command line")
            except:
                raise ValueError("Taxids should be comma-separated or path")
    else:
        taxids = None

    colormap = {}
    if args.get("colormap"):
        try:
            with open(args["colormap"]) as cm_f:
                colormap = {
                    line.split()[0]: line.split()[1].strip()
                    for line in cm_f
                    if line.strip()
                }
        except Exception as e:
            print(f"WARNING--Failed parsing colormap file: {e}")

    seqid2gene = {}
    filename_fasta = "%s.fa" % args["prefix"]

    if os.path.isfile(filename_fasta):
        filename_seqid2name = filename_fasta.replace(".fa", ".seqid2gname.tab")
        if os.path.isfile(filename_seqid2name):
            for line in open(filename_seqid2name):
                f = line.strip().split("\t")
                if len(f) >= 2:
                    seqid2gene[f[0]] = f[1]
    else:
        exclude_taxids = (
            [t.strip() for t in args["exclude_taxids"].split(",")]
            if args.get("exclude_taxids")
            else None
        )
        fasta, seqid2gene = get_fasta_from_db(
            version=args["version"],
            taxids=taxids,
            pfams=pfams,
            evalue_cutoff=args.get("evalue"),
            exclude_taxon_ids=exclude_taxids,
        )
        with open(filename_fasta, "w") as out:
            out.write(fasta)
        filename_seqid2name = filename_fasta.replace(".fa", ".seqid2gname.tab")
        with open(filename_seqid2name, "w") as out:
            for seqid, gname in seqid2gene.items():
                print("%s\t%s" % (seqid, gname), file=out)

    aln_cpu = args["cpu"] if args["cpu"] != "AUTO" else 32
    trimal_available = os.system("which trimal > /dev/null 2>&1") == 0

    def run_trimal(filename_aln):
        filename_trimal = filename_aln + ".gt%s" % args["gt"].replace(".", "")
        if trimal_available:
            if not os.path.isfile(filename_trimal):
                os.system(
                    "trimal -in %s -out %s -gt %s"
                    % (filename_aln, filename_trimal, args["gt"])
                )
        else:
            filename_trimal = filename_aln
        return filename_trimal

    if args["aln"] == "mafft":
        filename_aln = filename_fasta.replace(".fa", ".mft")
        if not os.path.isfile(filename_aln):
            os.system(
                "mafft --quiet --thread %s %s > %s"
                % (aln_cpu, filename_fasta, filename_aln)
            )
        filename_trimal = run_trimal(filename_aln)
    elif args["aln"] == "einsi":
        filename_aln = filename_fasta.replace(".fa", ".einsi")
        if not os.path.isfile(filename_aln):
            os.system(
                "einsi --thread %s %s > %s" % (aln_cpu, filename_fasta, filename_aln)
            )
        filename_trimal = run_trimal(filename_aln)
    elif args["aln"] == "clustalo":
        filename_aln = filename_fasta.replace(".fa", ".clustalo")
        if not os.path.isfile(filename_aln):
            os.system(
                "clustalo --threads %s -i %s -o %s"
                % (aln_cpu, filename_fasta, filename_aln)
            )
        filename_trimal = run_trimal(filename_aln)

    if args["ml"] == "fasttree":
        filename_tree = filename_trimal + ".lg.fasttree"
        if not os.path.isfile(filename_tree):
            ret = os.system("fasttree -lg %s > %s" % (filename_trimal, filename_tree))
            if ret != 0:
                sys.exit(1)
    elif args["ml"] == "iqtree":
        filename_tree = filename_trimal + ".treefile"
        if not os.path.isfile(filename_tree):
            ret = os.system(
                "iqtree -s %s --prefix %s -mset LG -B 1000 -T %s"
                % (filename_trimal, filename_trimal, args["cpu"])
            )
            if ret != 0:
                sys.exit(1)

    print("INFO--Loading tree in %s" % filename_tree)

    def generate_itol_color_strip(filename_tree, colormap, seqid2gene):
        itol_file = filename_tree + ".itol_colors.txt"
        unique_taxids = set()
        if seqid2gene:
            for seqid in seqid2gene.keys():
                unique_taxids.add(seqid.split(".")[0])
        else:
            try:
                with open(filename_tree) as tf:
                    raw_tree = tf.read()
                    found_ids = set(re.findall(r"([\d]+\.[a-zA-Z\d_]+)", raw_tree))
                    for seqid in found_ids:
                        unique_taxids.add(seqid.split(".")[0])
            except Exception as e:
                pass

        distinct_palette = [
            "#e6194B",
            "#3cb44b",
            "#ffe119",
            "#4363d8",
            "#f58231",
            "#911eb4",
            "#42d4f4",
            "#f032e6",
            "#bfef45",
            "#fabed4",
            "#469990",
            "#dcbeff",
            "#9A6324",
            "#fffac8",
            "#800000",
            "#aaffc3",
            "#808000",
            "#ffd8b1",
            "#000075",
            "#a9a9a9",
        ]
        auto_colormap = {
            taxid: distinct_palette[i % len(distinct_palette)]
            for i, taxid in enumerate(sorted(list(unique_taxids)))
        }

        with open(itol_file, "w") as f:
            f.write(
                "DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tTaxon Color Map\nCOLOR\t#ff0000\nDATA\n"
            )
            if seqid2gene:
                for seqid, name in seqid2gene.items():
                    f.write(
                        f"{seqid}\t{auto_colormap.get(seqid.split('.')[0], '#bcc3d0')}\t{seqid.split('.')[0]}\n"
                    )
            else:
                try:
                    with open(filename_tree) as tf:
                        for seqid in set(
                            re.findall(r"([\d]+\.[a-zA-Z\d_]+)", tf.read())
                        ):
                            f.write(
                                f"{seqid}\t{auto_colormap.get(seqid.split('.')[0], '#bcc3d0')}\t{seqid.split('.')[0]}\n"
                            )
                except Exception:
                    pass

    generate_itol_color_strip(filename_tree, colormap, seqid2gene)

    # ---------------------------------------------------------
    # --- ETE4 VISUAL EXPLORER ---
    # ---------------------------------------------------------
    if not args["no_explore"] and ETE_AVAILABLE:
        print(f"INFO--Starting ETE4 server on port {args['port']}")

        def _draw_leaf(node):
            if node.is_leaf:
                name_parts = node.name.split(".")
                display = name_parts[1] if len(name_parts) > 1 else node.name
                gene = node.props.get("gene_name", "")
                sci = node.props.get("sci_name", "")
                label = f"{display}  {gene}  [{sci}]".strip(" []")
                node.add_face(
                    TextFace(label, color="black"), column=0, position="branch_right"
                )

        leaf_name_layout = Layout(name="Leaf names", active=True, draw_node=_draw_leaf)

        # ── OTHER ETE4 LAYOUT SYNTAX (Using node.add_face) ──
        # def node_names_style(node):
        #     if node.is_leaf:
        #         try:
        #             sci = node.props.get('sci_name', '')
        #             gene = node.props.get('gene_name', '-')
        #             name_parts = node.name.split(".")
        #             display_name = name_parts[1] if len(name_parts) > 1 else node.name

        #             node.add_face(TextFace(display_name), column=0, position='right')
        #             node.add_face(TextFace(f" ({gene}) "), column=1, position='right')
        #             node.add_face(TextFace(f" ({sci}) "), column=2, position='right')
        #         except Exception:
        #             node.add_face(TextFace(node.name), column=0, position='right')

        def layout_seqface(node):
            if node.is_leaf:
                seq = node.props.get("seq")
                if seq:
                    node.add_face(
                        SeqFace(seq, seqtype="aa"), column=0, position="aligned"
                    )

        try:
            with open(filename_tree, "r") as f:
                nwk_str = f.read().strip()
            if not nwk_str.endswith(";"):
                nwk_str += ";"

            t = PhyloTree(
                nwk_str, sp_naming_function=lambda node: node.name.split(".")[0]
            )
            t.set_outgroup(t.get_midpoint_outgroup())
            t.resolve_polytomy(descendants=True)

            if not args["no_ncbi"]:
                try:
                    t.annotate_ncbi_taxa()
                except Exception as e:
                    print(f"WARNING--NCBI taxonomy annotation failed: {e}")

            # Attach sequences if MSA is present
            name2seq = {}
            if args.get("MSA") and os.path.isfile(filename_aln):
                name2seq = get_seqs(filename_aln)

            # ── ATTACH LABELS DIRECTLY TO NODES (Bypasses Layout Bugs) ──
            for node in t.traverse():
                if colormap and "lineage" in node.props:
                    for taxid in node.props["lineage"][::-1]:
                        if str(taxid) in colormap:
                            node.add_prop("color", colormap[str(taxid)])
                            break
                if node.is_leaf:
                    node.add_prop("gene_name", seqid2gene.get(node.name, "-"))
                    if args.get("MSA") and node.name in name2seq:
                        node.add_prop("seq", name2seq[node.name])
                        node.add_face(
                            SeqFace(name2seq[node.name], seqtype="aa"),
                            column=3,
                            position="aligned",
                        )

            t.explore(
                layouts=[BASIC_LAYOUT],
                keep_server=True,
                quiet=True,
                port=args["port"],
                host="0.0.0.0",
                show_leaf_name=True,
            )

        except Exception as e:
            print(f"ERROR--Failed to launch ETE4 explorer: {e}")
