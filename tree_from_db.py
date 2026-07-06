# ===========================================================================================================================

#                                                       IMPORTS

# ===========================================================================================================================


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
import shutil
import hashlib

os.environ["QT_QPA_PLATFORM"] = "offscreen"


try:
    from ete4 import PhyloTree
    from ete4.smartview import CircleFace
    from ete4.treeview import TreeStyle, SeqMotifFace, TextFace, NodeStyle

    # from ete4.smartview import Layout, TextFace, SeqFace, RectFace, BASIC_LAYOUT
    ETE_AVAILABLE = True
except ImportError as e:
    print(f"WARNING--Could not import ETE4 ({e}).")
    ETE_AVAILABLE = False
# --------------------------------------------------------


def get_fasta_from_db(
    version,
    taxids,
    pfams,
    evalue_cutoff=None,
    exclude_taxon_ids=None,
    max_per_taxon=None,
    pfam_source="hmmsearch",
    pfam_logic="or",
):
    taxid_ints = [int(t) for t in taxids] if taxids else None
    seen = {}
    acc2pfams = {}
    all_records = []

    with UniProtRetriever(get_db_config()) as db:
        for pfam in pfams:
            # uniprot source has no evalue/exclude filters
            if pfam_source == "uniprot":
                records = db.get_proteins(
                    version=version,
                    taxon_ids=taxid_ints,
                    pfam_id=pfam,
                )
            else:
                records = db.get_proteins_by_hmm_hit(
                    version=version,
                    hmm_query=pfam,
                    evalue_cutoff=evalue_cutoff,
                    taxon_ids=taxid_ints,
                    exclude_taxon_ids=exclude_taxon_ids,
                )
            for r in records:
                acc = r["accession"]
                if acc not in seen:
                    seen[acc] = r
                acc2pfams.setdefault(acc, set()).add(pfam)

    if pfam_logic == "and" and len(pfams) > 1:
        wanted = set(pfams)
        all_records = [seen[a] for a in seen if acc2pfams[a] == wanted]
        print(
            "INFO--Multi-Pfam AND logic: %d proteins carry all %d Pfams %s"
            % (len(all_records), len(pfams), ",".join(pfams))
        )
    else:
        all_records = list(seen.values())

    if len(pfams) > 1:
        from collections import Counter

        combo_counts = Counter(
            "+".join(sorted(acc2pfams[r["accession"]])) for r in all_records
        )
        print("INFO--Pfam combination breakdown (proteins entering the tree):")
        for combo, cnt in combo_counts.most_common():
            print("       %-40s %d" % (combo, cnt))

    print("INFO--Retrieved %d unique sequences from local DB" % len(all_records))

    if max_per_taxon is not None:
        from collections import defaultdict

        per_taxon = defaultdict(list)
        for r in all_records:
            per_taxon[r["taxon_id"]].append(r)
        all_records = []
        for taxid, recs in per_taxon.items():
            all_records.extend(recs[:max_per_taxon])
        print(
            "INFO--After subsampling (%d per taxon): %d sequences"
            % (max_per_taxon, len(all_records))
        )

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


# Color generation
def _domain_color(name):
    h = int(hashlib.md5(name.encode()).hexdigest()[:6], 16)
    return "#{:02x}{:02x}{:02x}".format((h >> 16) & 0xFF, (h >> 8) & 0xFF, h & 0xFF)


# ===========================================================================================================================

#                                                       ARGUMENTS

# ===========================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Uniprot/Pfam-based protein family evolution analysis"
    )
    parser.add_argument("--pfam", required=True, type=str)
    parser.add_argument("--taxids", required=False, type=str)
    parser.add_argument("--cpu", type=str, default="4")
    parser.add_argument("--ml", default="fasttree", choices=["fasttree", "iqtree"])
    parser.add_argument(
        "--aln", default="mafft", choices=["mafft", "einsi", "clustalo"]
    )
    parser.add_argument("--gt", type=str, default="0.01")
    parser.add_argument("--colormap", required=False, type=str)
    parser.add_argument("--version", required=True, type=str)
    parser.add_argument("--evalue", required=False, type=float, default=None)
    parser.add_argument("--local_fasta", required=False, type=str)
    parser.add_argument(
        "--pfam_source",
        default="hmmsearch",
        choices=["hmmsearch", "uniprot"],
        help="Where to pull Pfam membership from: 'hmmsearch' = local Pfam-A "
        "hmmsearch results (default), 'uniprot' = UniProt's own protein_pfam "
        "assignments.",
    )
    parser.add_argument(
        "--pfam_logic",
        default="or",
        choices=["or", "and"],
        help="With >1 Pfam: 'or' = proteins matching any Pfam (default), "
        "'and' = only proteins matching all Pfams.",
    )
    parser.add_argument(
        "--color_by",
        default="taxon",
        choices=["taxon", "pfam"],
        help="Branch colouring: 'taxon' = by taxonomy/lineage (default), "
        "'pfam' = by which queried Pfam(s) each protein carries.",
    )
    parser.add_argument("--MSA", action="store_true")
    parser.add_argument("--positions", required=False, type=str)
    parser.add_argument("--prefix", required=True, type=str)
    parser.add_argument("--port", required=False, type=int, default=5001)
    parser.add_argument("--no_ncbi", action="store_true")
    parser.add_argument("--no_explore", action="store_true")
    parser.add_argument("--exclude_taxids", required=False, type=str)
    parser.add_argument("--output_dir", required=False, type=str, default=None)
    parser.add_argument(
        "--max_per_taxon",
        required=False,
        type=int,
        default=None,
        help="Keep at most N sequences per taxon before alignment. "
        "Useful for large gene families to keep the tree manageable.",
    )
    parser.add_argument("--render_ete_static", action="store_true")
    parser.add_argument(
        "--static_layers", required=False, type=str, default="names,domains,colors,gene"
    )
    parser.add_argument("--use_resolved", action="store_true")

    args = vars(parser.parse_args())

    if args.get("output_dir"):
        os.makedirs(args["output_dir"], exist_ok=True)
        args["prefix"] = os.path.join(args["output_dir"], args["prefix"])

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
    elif args.get("local_fasta"):
        # headers must be "{taxon_id}.{accession}"
        with open(args["local_fasta"]) as src, open(filename_fasta, "w") as out:
            out.write(src.read())
        for head in get_seqs(filename_fasta):
            seqid2gene[head] = "-"
        filename_seqid2name = filename_fasta.replace(".fa", ".seqid2gname.tab")
        with open(filename_seqid2name, "w") as out:
            for seqid, gname in seqid2gene.items():
                print("%s\t%s" % (seqid, gname), file=out)
        print(
            "INFO--Loaded %d sequences from local FASTA: %s"
            % (len(seqid2gene), args["local_fasta"])
        )
    else:
        if type(args["exclude_taxids"]) == str:
            if os.path.isfile(args["exclude_taxids"]):
                exclude_taxids = [
                    line.strip().split()[0]
                    for line in open(args["exclude_taxids"])
                    if not line.startswith("#")
                ]
            else:
                try:
                    exclude_taxids = args["exclude_taxids"].split(",")
                    print("Taxids to be excluded loaded from command line")
                except:
                    raise ValueError(
                        "Taxids to be excluded should be comma-separated or path"
                    )
        else:
            exclude_taxids = None
        fasta, seqid2gene = get_fasta_from_db(
            version=args["version"],
            taxids=taxids,
            pfams=pfams,
            evalue_cutoff=args.get("evalue"),
            exclude_taxon_ids=exclude_taxids,
            max_per_taxon=args.get("max_per_taxon"),
            pfam_source=args.get("pfam_source", "hmmsearch"),
            pfam_logic=args.get("pfam_logic", "or"),
        )
        with open(filename_fasta, "w") as out:
            out.write(fasta)
        filename_seqid2name = filename_fasta.replace(".fa", ".seqid2gname.tab")
        with open(filename_seqid2name, "w") as out:
            for seqid, gname in seqid2gene.items():
                print("%s\t%s" % (seqid, gname), file=out)

    aln_cpu = str(args["cpu"]) if args["cpu"] != "AUTO" else "4"
    trimal_available = shutil.which("trimal") is not None

    def run_trimal(filename_aln):
        filename_trimal = filename_aln + ".gt%s" % args["gt"].replace(".", "")
        if trimal_available:
            if not os.path.isfile(filename_trimal):
                subprocess.run(
                    [
                        "trimal",
                        "-in",
                        filename_aln,
                        "-out",
                        filename_trimal,
                        "-gt",
                        args["gt"],
                    ],
                    check=True,
                )
        else:
            filename_trimal = filename_aln
        return filename_trimal

    # ===========================================================================================================================

    #                                                       ALIGNMENT

    # ===========================================================================================================================

    if args["aln"] == "mafft":
        filename_aln = filename_fasta.replace(".fa", ".mft")
        if not os.path.isfile(filename_aln):
            try:
                with open(filename_aln, "w") as out:
                    subprocess.run(
                        ["mafft", "--quiet", "--thread", aln_cpu, filename_fasta],
                        stdout=out,
                        check=True,
                    )
            except subprocess.CalledProcessError:
                if os.path.isfile(filename_aln):
                    os.remove(filename_aln)
                sys.exit("ERROR--mafft alignment failed.")
        filename_trimal = run_trimal(filename_aln)
    elif args["aln"] == "einsi":
        filename_aln = filename_fasta.replace(".fa", ".einsi")
        if not os.path.isfile(filename_aln):
            try:
                with open(filename_aln, "w") as out:
                    subprocess.run(
                        ["einsi", "--thread", aln_cpu, filename_fasta],
                        stdout=out,
                        check=True,
                    )
            except subprocess.CalledProcessError:
                if os.path.isfile(filename_aln):
                    os.remove(filename_aln)
                sys.exit("ERROR--einsi alignment failed.")
        filename_trimal = run_trimal(filename_aln)
    elif args["aln"] == "clustalo":
        filename_aln = filename_fasta.replace(".fa", ".clustalo")
        if not os.path.isfile(filename_aln):
            try:
                subprocess.run(
                    [
                        "clustalo",
                        "--threads",
                        aln_cpu,
                        "-i",
                        filename_fasta,
                        "-o",
                        filename_aln,
                    ],
                    check=True,
                )
            except subprocess.CalledProcessError:
                if os.path.isfile(filename_aln):
                    os.remove(filename_aln)
                sys.exit("ERROR--clustalo alignment failed.")
        filename_trimal = run_trimal(filename_aln)

    # ===========================================================================================================================

    #                                                       TREE BUILDING

    # ===========================================================================================================================

    if args["ml"] == "fasttree":
        filename_tree = filename_trimal + ".lg.fasttree"
        if not os.path.isfile(filename_tree):
            try:
                with open(filename_tree, "w") as out:
                    subprocess.run(
                        ["fasttree", "-lg", filename_trimal],
                        stdout=out,
                        check=True,
                    )
            except subprocess.CalledProcessError:
                if os.path.isfile(filename_tree):
                    os.remove(filename_tree)
                sys.exit("ERROR--FastTree failed.")
    elif args["ml"] == "iqtree":
        filename_tree = filename_trimal + ".treefile"
        if not os.path.isfile(filename_tree):
            try:
                subprocess.run(
                    [
                        "iqtree",
                        "-s",
                        filename_trimal,
                        "--prefix",
                        filename_trimal,
                        "-mset",
                        "LG",
                        "-B",
                        "1000",
                        "-T",
                        aln_cpu,
                    ],
                    check=True,
                )
            except subprocess.CalledProcessError:
                sys.exit("ERROR--IQ-TREE failed.")

    print("INFO--Loading tree in %s" % filename_tree)

    # ===========================================================================================================================

    #                                                       ITOL COLORMAP

    # ===========================================================================================================================

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

        auto_colormap = {taxid: _domain_color(taxid) for taxid in sorted(unique_taxids)}

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

    # ===========================================================================================================================

    #                                                       ITOL DOMAINS

    # ===========================================================================================================================

    def generate_itol_domains(filename_tree, domain_dict, name2seq):
        itol_file = filename_tree + ".itol_domains.txt"

        names = sorted({d["hmm_name"] for hits in domain_dict.values() for d in hits})

        with open(itol_file, "w") as f:
            f.write(
                "DATASET_DOMAINS\nSEPARATOR TAB\nDATASET_LABEL\tPfam domains\nCOLOR\t#000000\n"
            )
            f.write("LEGEND_TITLE\tDomains\n")
            f.write("LEGEND_SHAPES\t" + "\t".join(["RE"] * len(names)) + "\n")
            f.write(
                "LEGEND_COLORS\t" + "\t".join(_domain_color(n) for n in names) + "\n"
            )
            f.write("LEGEND_LABELS\t" + "\t".join(names) + "\n")
            f.write("DATA\n")
            for seqid, seq in name2seq.items():
                acc = seqid.split(".", 1)[1] if "." in seqid else seqid
                hits = domain_dict.get(acc, [])
                if not hits:
                    continue
                fields = [seqid, str(len(seq))]
                for d in hits:
                    fields.append(
                        f"RE|{d['ali_from']}|{d['ali_to']}|{_domain_color(d['hmm_name'])}|{d['hmm_name']}"
                    )
                f.write("\t".join(fields) + "\n")

    # ===========================================================================================================================

    #                                                           ETE4

    # ===========================================================================================================================
    def get_species_name(node):
        return node.name.split(".")[0]

    if ETE_AVAILABLE:
        try:
            if args.get("use_resolved"):
                resolved_path = filename_tree + ".resolved.nwk"
                if os.path.exists(resolved_path):
                    filename_tree = resolved_path
            with open(filename_tree, "r") as f:
                nwk_str = f.read().strip()
            if not nwk_str.endswith(";"):
                nwk_str += ";"

            t = PhyloTree(nwk_str, sp_naming_function=get_species_name)
            if not args.get("use_resolved"):
                t.set_outgroup(t.get_midpoint_outgroup())
                t.resolve_polytomy(descendants=True)

            # Annotate internal nodes with evoltype ("D"=duplication, "S"=speciation)
            try:
                t.get_descendant_evol_events()
                print(
                    "INFO--Ortholog/paralog annotation complete (evoltype set on internal nodes)"
                )
            except Exception as e:
                print(f"WARNING--Could not annotate evol events: {e}")

            # Fetch Domain Data
            domain_dict = {}
            try:
                leaf_accs = [
                    n.name.split(".", 1)[1] for n in t.leaves() if "." in n.name
                ]
                if leaf_accs:
                    hits = fetch_domains_by_accession(args["version"], leaf_accs)
                    for d in hits:
                        domain_dict.setdefault(d["accession"], []).append(d)
                    for acc in domain_dict:
                        domain_dict[acc].sort(key=lambda x: x["ali_from"])
                print(f"INFO--Domain data loaded for {len(domain_dict)} proteins")
            except Exception as e:
                print(f"WARNING--Could not fetch domain data: {e}")

            if domain_dict:
                generate_itol_domains(
                    filename_tree, domain_dict, get_seqs(filename_fasta)
                )
                print(
                    f"INFO--iTOL domains file written: {filename_tree}.itol_domains.txt"
                )

            # optional MSA display window from --positions ("start:end")
            msa_start, msa_end = None, None
            if args.get("positions"):
                sep = ":" if ":" in args["positions"] else "-"
                parts = args["positions"].split(sep)
                if len(parts) == 2:
                    msa_start = int(parts[0]) if parts[0].strip() else None
                    msa_end = int(parts[1]) if parts[1].strip() else None

            # map each leaf to the queried Pfams it carries
            color_by = args.get("color_by", "taxon")
            acc2combo = {}
            for acc, hits in domain_dict.items():
                matched = []
                for p in pfams:
                    for d in hits:
                        if d["hmm_name"] == p or str(
                            d.get("hmm_accession", "")
                        ).startswith(p):
                            matched.append(p)
                            break
                acc2combo[acc] = "+".join(matched) if matched else "none"
            combo_colors = {
                combo: _domain_color(combo) for combo in set(acc2combo.values())
            }
            if len(pfams) > 1 and acc2combo:
                combo_file = filename_tree + ".pfam_combinations.tab"
                with open(combo_file, "w") as cf:
                    for n in t.leaves():
                        acc = n.name.split(".", 1)[1] if "." in n.name else n.name
                        cf.write(
                            "%s\t%s\t%s\n" % (n.name, acc, acc2combo.get(acc, "none"))
                        )
                print(f"INFO--Pfam combination table written: {combo_file}")

            # ---------------------------------------------------------------------------------------
            # MODE A: STATIC PNG GENERATOR WITH DOMAINS DISPLAYED
            # ---------------------------------------------------------------------------------------
            if args.get("render_ete_static"):
                print("INFO--Generating static tree image with custom domain shapes...")

                t_static = t.copy()
                # unaligned for domain coords, aligned for the MSA layer
                name2seq = get_seqs(filename_fasta)
                name2aln = (
                    get_seqs(filename_aln) if os.path.isfile(filename_aln) else {}
                )

                # Parse layer flags
                static_layers = set(
                    args.get("static_layers", "names,domains,colors,gene").split(",")
                )
                print(f"INFO--Static layers: {static_layers}")

                # Build auto-colormap from taxids (same palette as iTOL)
                auto_colormap = {}
                if "colors" in static_layers:
                    if not args.get("no_ncbi"):
                        try:
                            t_static.annotate_ncbi_taxa()
                        except Exception as e:
                            print(f"WARNING--NCBI annotation failed for static: {e}")
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
                        "#469990",
                        "#dcbeff",
                        "#9A6324",
                        "#800000",
                        "#aaffc3",
                        "#808000",
                        "#ffd8b1",
                        "#000075",
                        "#a9a9a9",
                    ]
                    unique_taxids = sorted(
                        {n.name.split(".")[0] for n in t_static.leaves()}
                    )
                    auto_colormap = {
                        tid: distinct_palette[i % len(distinct_palette)]
                        for i, tid in enumerate(unique_taxids)
                    }
                    # user-supplied colormap (taxid -> color) overrides the palette
                    if colormap:
                        auto_colormap.update(colormap)

                def _get_static_shape(domain_name):
                    shapes = ["[]", "()", "<>", "^", "v", "o"]
                    shape_idx = int(
                        hashlib.md5(domain_name.encode()).hexdigest()[-2:], 16
                    ) % len(shapes)
                    return shapes[shape_idx]

                for node in t_static.traverse():

                    nstyle = NodeStyle()
                    nstyle["hz_line_width"] = 4
                    nstyle["vt_line_width"] = 4
                    nstyle["size"] = 0

                    col = None
                    if color_by == "pfam" and combo_colors:
                        if node.is_leaf:
                            acc = (
                                node.name.split(".", 1)[1]
                                if "." in node.name
                                else node.name
                            )
                            key = acc2combo.get(acc)
                        else:
                            keys = {
                                acc2combo.get(
                                    l.name.split(".", 1)[1] if "." in l.name else l.name
                                )
                                for l in node.leaves()
                            }
                            key = next(iter(keys)) if len(keys) == 1 else None
                        col = combo_colors.get(key) if key else None
                    elif "colors" in static_layers and auto_colormap:
                        if node.is_leaf:
                            key = node.name.split(".")[0]
                        else:
                            keys = {l.name.split(".")[0] for l in node.leaves()}
                            key = next(iter(keys)) if len(keys) == 1 else None
                        col = auto_colormap.get(key) if key else None

                    if col:
                        nstyle["vt_line_color"] = col
                        nstyle["hz_line_color"] = col
                    node.set_style(nstyle)

                    if node.is_leaf:
                        name_parts = node.name.split(".")
                        accession = name_parts[1] if len(name_parts) > 1 else node.name

                        # --- Leaf name + gene label ---
                        if "names" in static_layers:
                            display = accession
                            gene = (
                                seqid2gene.get(node.name, "")
                                if "gene" in static_layers
                                else ""
                            )
                            node.add_face(
                                TextFace(f"{display}  {gene}  "),
                                column=0,
                                position="branch-right",
                            )

                        # --- Domain shapes ---
                        if "domains" in static_layers:
                            seq = name2seq.get(node.name, None)
                            domains = domain_dict.get(accession, [])
                            if seq and domains:
                                motifs = []
                                for d in domains:
                                    motifs.append(
                                        [
                                            d["ali_from"],
                                            d["ali_to"],
                                            _get_static_shape(d["hmm_name"]),
                                            None,
                                            14,
                                            "black",
                                            _domain_color(d["hmm_name"]),
                                            f"arial|8|black|{d['hmm_name']}",
                                        ]
                                    )
                                node.add_face(
                                    SeqMotifFace(seq, motifs=motifs, seq_format="-"),
                                    column=1,
                                    position="aligned",
                                )

                        # --- MSA aligned sequences ---
                        if "msa" in static_layers:
                            seq = name2aln.get(node.name, None)
                            if seq:
                                seq = seq[msa_start:msa_end]
                                node.add_face(
                                    SeqMotifFace(seq, motifs=[], seq_format="seq"),
                                    column=2,
                                    position="aligned",
                                )

                ts = TreeStyle()
                ts.show_leaf_name = False
                output_img = f"{args['prefix']}_tree_domains.png"

                n_leaves = len([n for n in t_static.traverse() if n.is_leaf])
                height = min(max(n_leaves * 12, 800), 60000)

                t_static.render(output_img, w=4000, units="px", tree_style=ts)
                print(f"Wrote file: {output_img}")
                print(f"INFO--Success. Static tree image saved to: {output_img}")

            # ---------------------------------------------------------------------------------------
            # MODE B: INTERACTIVE EXPLORER
            # ---------------------------------------------------------------------------------------
            elif not args.get("no_explore"):
                print(f"INFO--Starting ETE4 server on port {args['port']}")

                from ete4.smartview import (
                    Layout,
                    TextFace as SmartTextFace,
                    SeqFace,
                    RectFace,
                    BASIC_LAYOUT,
                )

                def _draw_leaf(node):
                    if not node.is_leaf:
                        return
                    name_parts = node.name.split(".")
                    display = name_parts[1] if len(name_parts) > 1 else node.name
                    gene = node.props.get("gene_name", "")
                    sci = node.props.get("sci_name", "")
                    label = f"{display}  {gene}  [{sci}]".strip(" []")
                    return SmartTextFace(
                        label, style="fill: black;", position="right", column=0
                    )

                leaf_name_layout = Layout(
                    name="Leaf names", active=True, draw_node=_draw_leaf
                )

                def layout_seqface(node):
                    if not node.is_leaf:
                        return
                    seq = node.props.get("seq")
                    if seq:
                        seq = seq[msa_start:msa_end]
                        return [
                            SeqFace(seq, seqtype="aa", position="aligned", column=4)
                        ]

                # Build aliases — one per unique domain name
                unique_domain_names = sorted(
                    {d["hmm_name"] for hits in domain_dict.values() for d in hits}
                )
                name_to_alias = {
                    n: f"dom-{i}" for i, n in enumerate(unique_domain_names)
                }
                domain_aliases = {
                    name_to_alias[n]: {
                        "fill": _domain_color(n),
                        "stroke": "black",
                        "stroke-width": "1",
                    }
                    for n in unique_domain_names
                }

                def _draw_domains(node):
                    if not node.is_leaf:
                        return
                    accession = (
                        node.name.split(".", 1)[1] if "." in node.name else node.name
                    )
                    domains = domain_dict.get(accession, [])
                    if not domains:
                        return
                    seq_len = max(d["ali_to"] for d in domains)
                    faces = []
                    for i, d in enumerate(domains):
                        width = max(
                            30, int((d["ali_to"] - d["ali_from"]) / seq_len * 200)
                        )
                        faces.append(
                            RectFace(
                                wmax=width,
                                hmax=18,
                                style=name_to_alias[d["hmm_name"]],
                                text=d["hmm_name"],
                                position="aligned",
                                column=i + 1,
                            )
                        )
                    return faces

                domain_layout = Layout(
                    name="Domains",
                    active=True,
                    draw_tree={"aliases": domain_aliases},
                    draw_node=_draw_domains,
                )

                if not args.get("no_ncbi"):
                    try:
                        t.annotate_ncbi_taxa()
                    except Exception as e:
                        print(f"WARNING--NCBI taxonomy annotation failed: {e}")

                # Attach sequences if MSA is present
                name2seq = {}
                if args.get("MSA") and os.path.isfile(filename_aln):
                    name2seq = get_seqs(filename_aln)
                seq_layout = Layout(
                    name="MSA",
                    active=bool(args.get("MSA")),
                    draw_node=layout_seqface,
                )

                def _draw_branch_color(node):
                    style = {"stroke-width": 5}
                    col = node.props.get("branch_color")
                    if col:
                        style["stroke"] = col
                    return {"hz-line": style, "vt-line": style}

                branch_color_layout = Layout(
                    name="Branch colors", active=True, draw_node=_draw_branch_color
                )

                # ORTHOLOGY / PARALOGY LAYOUT
                def _draw_ortho_para(node):
                    evol_type = node.props.get("evoltype", None)
                    if evol_type == "D":
                        return {
                            "dot": {
                                "shape": "circle",
                                "fill": "#c81e1e",
                                "stroke": "#960000",
                                "stroke-width": 1.5,
                                "radius": 7,
                            }
                        }
                    elif evol_type == "S":
                        return {
                            "dot": {
                                "shape": "circle",
                                "fill": "#1e1eb4",
                                "stroke": "#000082",
                                "stroke-width": 1.5,
                                "radius": 7,
                            }
                        }
                    return None

                ortho_para_layout = Layout(
                    name="ortho/para", active=False, draw_node=_draw_ortho_para
                )

                interactive_colormap = dict(colormap) if colormap else {}
                if color_by == "taxon" and not interactive_colormap:
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
                        "#469990",
                        "#dcbeff",
                        "#9A6324",
                        "#800000",
                        "#aaffc3",
                        "#808000",
                        "#ffd8b1",
                        "#000075",
                        "#a9a9a9",
                    ]
                    unique_taxids = sorted({n.name.split(".")[0] for n in t.leaves()})
                    interactive_colormap = {
                        tid: distinct_palette[i % len(distinct_palette)]
                        for i, tid in enumerate(unique_taxids)
                    }

                for node in t.traverse():
                    col = None
                    if color_by == "pfam" and combo_colors:
                        if node.is_leaf:
                            acc = (
                                node.name.split(".", 1)[1]
                                if "." in node.name
                                else node.name
                            )
                            key = acc2combo.get(acc)
                        else:
                            keys = {
                                acc2combo.get(
                                    l.name.split(".", 1)[1] if "." in l.name else l.name
                                )
                                for l in node.leaves()
                            }
                            key = next(iter(keys)) if len(keys) == 1 else None
                        col = combo_colors.get(key) if key else None
                    elif color_by == "taxon" and interactive_colormap:
                        if node.is_leaf:
                            key = node.name.split(".")[0]
                        else:
                            keys = {l.name.split(".")[0] for l in node.leaves()}
                            key = next(iter(keys)) if len(keys) == 1 else None
                        col = interactive_colormap.get(key) if key else None
                    if col:
                        node.add_prop("branch_color", col)
                    if node.is_leaf:
                        parts = node.name.split(".", 1)
                        acc = parts[1] if len(parts) > 1 else node.name
                        node.add_prop("taxid", parts[0])
                        node.add_prop("accession", acc)
                        node.add_prop("gene_name", seqid2gene.get(node.name, "-"))
                        if acc in acc2combo:
                            node.add_prop("pfam_combo", acc2combo[acc])
                        if domain_dict.get(acc):
                            node.add_prop(
                                "domains",
                                ", ".join(d["hmm_name"] for d in domain_dict[acc]),
                            )
                        if args.get("MSA") and node.name in name2seq:
                            node.add_prop("seq", name2seq[node.name])

                # ATTACH NODE PATHS FOR GUI VERIFICATION
                def attach_paths(node, current_path=""):
                    node.add_prop("node_path", current_path if current_path else "root")
                    for i, child in enumerate(node.children):
                        next_path = f"{current_path},{i}" if current_path else str(i)
                        attach_paths(child, next_path)

                attach_paths(t)

                # ETE4 SmartView raises "Cannot draw tree with width 0" when
                # any branch has dist=0 (e.g. polytomy-resolution nodes added
                # by tree_builder.py). Replace zeros with a tiny value.
                for _n in t.traverse():
                    if _n.up is not None and (not _n.dist or _n.dist <= 0):
                        _n.dist = 1e-6

                t.explore(
                    layouts=[
                        BASIC_LAYOUT,
                        leaf_name_layout,
                        domain_layout,
                        seq_layout,
                        branch_color_layout,
                        ortho_para_layout,
                    ],
                    keep_server=True,
                    quiet=True,
                    open_browser=False,
                    port=args["port"],
                    host="0.0.0.0",
                    show_leaf_name=True,
                    show_popup_props=[
                        "name",
                        "sci_name",
                        "taxid",
                        "accession",
                        "gene_name",
                        "pfam_combo",
                        "domains",
                        "rank",
                        "dist",
                        "support",
                        "evoltype",
                    ],
                )

        except Exception as e:
            print(f"ERROR--Failed running ETE4 Graphics Engine: {e}")
