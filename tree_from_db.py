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
import hashlib
os.environ["QT_QPA_PLATFORM"] = "offscreen"


try:
    from ete4 import PhyloTree
    from ete4.treeview import TreeStyle, SeqMotifFace, TextFace, NodeStyle

    # from ete4.smartview import Layout, TextFace, SeqFace, RectFace, BASIC_LAYOUT
    ETE_AVAILABLE = True
except ImportError as e:
    print(f"WARNING--Could not import ETE4 ({e}).")
    ETE_AVAILABLE = False
# --------------------------------------------------------

def get_fasta_from_db(version, taxids, pfams, evalue_cutoff=None, exclude_taxon_ids=None):
    taxid_ints = [int(t) for t in taxids] if taxids else None
    seen        = set()
    all_records = []

    with UniProtRetriever(get_db_config()) as db:
        for pfam in pfams:
            records = db.get_proteins_by_hmm_hit(
                version       = version,
                hmm_query     = pfam,
                evalue_cutoff = evalue_cutoff,
                taxon_ids     = taxid_ints,
                exclude_taxon_ids = exclude_taxon_ids,
            )
            for r in records:
                if r["accession"] not in seen:
                    seen.add(r["accession"])
                    all_records.append(r)

    print("INFO--Retrieved %d unique sequences from local DB" % len(all_records))

    fasta_lines = []
    seqid2gene  = {}
    for r in all_records:
        seqid = "%s.%s" % (r["taxon_id"], r["accession"])
        fasta_lines.append(">%s\n%s" % (seqid, r["sequence"]))
        seqid2gene[seqid] = r["name"]

    return "\n".join(fasta_lines) + "\n", seqid2gene

def get_seqs(fastafile):
    name2seq = {}
    seq = ''
    head = ''
    for line in open(fastafile):
        if line.startswith('>'):
            if seq:
                name2seq[head] = seq
                seq = ''
                head = line.lstrip('>').rstrip()
            else:
                head = line.lstrip('>').rstrip()
        else:
            seq += line.rstrip()
    if head:
        name2seq[head] = seq
    return name2seq

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Uniprot/Pfam-based protein family evolution analysis')
    parser.add_argument('--pfam', required=True, type=str)
    parser.add_argument('--taxids', required=False, type=str)
    parser.add_argument('--cpu', type=str, default='4')
    parser.add_argument('--ml', default="fasttree", choices=['fasttree','iqtree'])
    parser.add_argument('--aln', default="mafft", choices=['mafft','einsi','clustalo'])
    parser.add_argument('--gt', type=str, default="0.01")
    parser.add_argument('--colormap', required=False, type=str)
    parser.add_argument('--version', required=True, type=str)
    parser.add_argument('--evalue', required=False, type=float, default=None)
    parser.add_argument('--local_fasta', required=False, type=str)
    parser.add_argument('--MSA', action='store_true')
    parser.add_argument('--positions', required=False, type=str)
    parser.add_argument('--prefix', required=True, type=str)
    parser.add_argument('--port', required=False, type=int, default=5001)
    parser.add_argument('--no_ncbi', action='store_true')
    parser.add_argument('--no_explore', action='store_true')
    parser.add_argument('--exclude_taxids', required=False, type=str)
    parser.add_argument('--output_dir', required=False, type=str, default=None)
    parser.add_argument('--render_ete_static', action='store_true')
    parser.add_argument('--static_layers', required=False, type=str,
                    default="names,domains,colors,gene")
    

    
    args = vars(parser.parse_args())
    
    if args.get("output_dir"):
        os.makedirs(args["output_dir"], exist_ok=True)
        args["prefix"] = os.path.join(args["output_dir"], args["prefix"])

    print("INFO--Processing Pfam domains:", args['pfam'])
    pfams = args['pfam'].split(',')
    
    if type(args['taxids']) == str:
        if os.path.isfile(args['taxids']):
            taxids = [line.strip().split()[0] for line in open(args['taxids']) if not line.startswith('#')]
        else:
            try:
                taxids = args['taxids'].split(',')
                print('Taxids loaded from command line')
            except:
                raise ValueError('Taxids should be comma-separated or path')
    else:
        taxids = None

    colormap = {}
    if args.get('colormap'):
        try:
            with open(args['colormap']) as cm_f:
                colormap = {line.split()[0]: line.split()[1].strip() for line in cm_f if line.strip()}
        except Exception as e:
            print(f"WARNING--Failed parsing colormap file: {e}")

    seqid2gene = {}
    filename_fasta = "%s.fa" % args['prefix']

    if os.path.isfile(filename_fasta):
        filename_seqid2name = filename_fasta.replace(".fa", ".seqid2gname.tab")
        if os.path.isfile(filename_seqid2name):
            for line in open(filename_seqid2name):
                f = line.strip().split('\t')
                if len(f) >= 2:
                    seqid2gene[f[0]] = f[1]
    else:
        exclude_taxids = [t.strip() for t in args['exclude_taxids'].split(',')] if args.get('exclude_taxids') else None
        fasta, seqid2gene = get_fasta_from_db(
            version       = args['version'],
            taxids        = taxids,
            pfams         = pfams,
            evalue_cutoff = args.get('evalue'),
            exclude_taxon_ids = exclude_taxids,
        )
        with open(filename_fasta, "w") as out:
            out.write(fasta)
        filename_seqid2name = filename_fasta.replace(".fa", ".seqid2gname.tab")
        with open(filename_seqid2name, 'w') as out:
            for seqid, gname in seqid2gene.items():
                print('%s\t%s' % (seqid, gname), file=out)

    aln_cpu = args['cpu'] if args['cpu'] != 'AUTO' else 4
    trimal_available = os.system("which trimal > /dev/null 2>&1") == 0

    def run_trimal(filename_aln):
        filename_trimal = filename_aln + ".gt%s" % args['gt'].replace(".", "")
        if trimal_available:
            if not os.path.isfile(filename_trimal):
                os.system("trimal -in %s -out %s -gt %s" % (filename_aln, filename_trimal, args['gt']))
        else:
            filename_trimal = filename_aln
        return filename_trimal


# ===========================================================================================================================

#                                                       ALIGNMENT

# ===========================================================================================================================

    if args['aln'] == 'mafft':
        filename_aln = filename_fasta.replace(".fa", ".mft")
        if not os.path.isfile(filename_aln):
            os.system("mafft --quiet --thread %s %s > %s" % (aln_cpu, filename_fasta, filename_aln))
        filename_trimal = run_trimal(filename_aln)
    elif args['aln'] == 'einsi':
        filename_aln = filename_fasta.replace(".fa", ".einsi")
        if not os.path.isfile(filename_aln):
            os.system("einsi --thread %s %s > %s" % (aln_cpu, filename_fasta, filename_aln))
        filename_trimal = run_trimal(filename_aln)
    elif args['aln'] == 'clustalo':
        filename_aln = filename_fasta.replace(".fa", ".clustalo")
        if not os.path.isfile(filename_aln):
            os.system("clustalo --threads %s -i %s -o %s" % (aln_cpu, filename_fasta, filename_aln))
        filename_trimal = run_trimal(filename_aln)

    if args['ml'] == 'fasttree':
        filename_tree = filename_trimal + ".lg.fasttree"
        if not os.path.isfile(filename_tree):
            ret = os.system("fasttree -lg %s > %s" % (filename_trimal, filename_tree))
            if ret != 0: sys.exit(1)
    elif args['ml'] == 'iqtree':
        filename_tree = filename_trimal + ".treefile"
        if not os.path.isfile(filename_tree):
            ret = os.system("iqtree -s %s --prefix %s -mset LG -B 1000 -T %s" % (filename_trimal, filename_trimal, args['cpu']))
            if ret != 0: sys.exit(1)

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
                    found_ids = set(re.findall(r'([\d]+\.[a-zA-Z\d_]+)', raw_tree))
                    for seqid in found_ids:
                        unique_taxids.add(seqid.split(".")[0])
            except Exception as e:
                pass

        distinct_palette = ["#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990", "#dcbeff", "#9A6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9"]
        auto_colormap = {taxid: distinct_palette[i % len(distinct_palette)] for i, taxid in enumerate(sorted(list(unique_taxids)))}

        with open(itol_file, "w") as f:
            f.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tTaxon Color Map\nCOLOR\t#ff0000\nDATA\n")
            if seqid2gene:
                for seqid, name in seqid2gene.items():
                    f.write(f"{seqid}\t{auto_colormap.get(seqid.split('.')[0], '#bcc3d0')}\t{seqid.split('.')[0]}\n")
            else:
                try:
                    with open(filename_tree) as tf:
                        for seqid in set(re.findall(r'([\d]+\.[a-zA-Z\d_]+)', tf.read())):
                            f.write(f"{seqid}\t{auto_colormap.get(seqid.split('.')[0], '#bcc3d0')}\t{seqid.split('.')[0]}\n")
                except Exception:
                    pass

    generate_itol_color_strip(filename_tree, colormap, seqid2gene)
    
# ===========================================================================================================================

#                                                       ITOL DOMAINS

# ===========================================================================================================================
    
    def generate_itol_domains(filename_tree, domain_dict, name2seq):
        itol_file = filename_tree + ".itol_domains.txt"
        palette = ["#e6194B","#3cb44b","#4363d8","#f58231","#911eb4",
                "#42d4f4","#f032e6","#469990","#9A6324","#800000"]
        names = sorted({d["hmm_name"] for hits in domain_dict.values() for d in hits})
        color_of = {n: palette[i % len(palette)] for i, n in enumerate(names)}

        with open(itol_file, "w") as f:
            f.write("DATASET_DOMAINS\nSEPARATOR TAB\nDATASET_LABEL\tPfam domains\nCOLOR\t#000000\n")
            f.write("LEGEND_TITLE\tDomains\n")
            f.write("LEGEND_SHAPES\t" + "\t".join(["RE"] * len(names)) + "\n")
            f.write("LEGEND_COLORS\t" + "\t".join(color_of[n] for n in names) + "\n")
            f.write("LEGEND_LABELS\t" + "\t".join(names) + "\n")
            f.write("DATA\n")
            for seqid, seq in name2seq.items():
                acc = seqid.split(".", 1)[1] if "." in seqid else seqid
                hits = domain_dict.get(acc, [])
                if not hits:
                    continue
                fields = [seqid, str(len(seq))]
                for d in hits:
                    fields.append(f"RE|{d['ali_from']}|{d['ali_to']}|{color_of[d['hmm_name']]}|{d['hmm_name']}")
                f.write("\t".join(fields) + "\n")

# ===========================================================================================================================

#                                                           ETE4 

# ===========================================================================================================================
    def get_species_name(node):
        return node.name.split('.')[0]
    
    if ETE_AVAILABLE:
        try:
            with open(filename_tree, 'r') as f:
                nwk_str = f.read().strip()
            if not nwk_str.endswith(';'):
                nwk_str += ';'
                
            t = PhyloTree(nwk_str, sp_naming_function=get_species_name)
            t.set_outgroup(t.get_midpoint_outgroup())
            t.resolve_polytomy(descendants=True)
            
            # Fetch Domain Data
            domain_dict = {}
            try:
                leaf_accs = [n.name.split(".", 1)[1] for n in t.leaves() if "." in n.name]
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
                generate_itol_domains(filename_tree, domain_dict, get_seqs(filename_fasta))
                print(f"INFO--iTOL domains file written: {filename_tree}.itol_domains.txt")

            # Color generation 
            def _domain_color(name):
                h = int(hashlib.md5(name.encode()).hexdigest()[:6], 16)
                return "#{:02x}{:02x}{:02x}".format((h >> 16) & 0xFF, (h >> 8) & 0xFF, h & 0xFF)
            
            
            
            
            # ---------------------------------------------------------------------------------------
            # MODE A: STATIC PNG GENERATOR WITH DOMAINS DISPLAYED
            # ---------------------------------------------------------------------------------------
            if args.get("render_ete_static"):
                print("INFO--Generating static tree image with custom domain shapes...")
                
                t_static = t.copy()
                name2seq = get_seqs(filename_aln) if (args.get('MSA') and os.path.isfile(filename_aln)) else get_seqs(filename_fasta)
                
                # Parse layer flags
                static_layers = set(args.get("static_layers", "names,domains,colors,gene").split(","))
                print(f"INFO--Static layers: {static_layers}")

                # Build auto-colormap from taxids (same palette as iTOL)
                auto_colormap = {}
                if "colors" in static_layers:
                    if not args.get('no_ncbi'):
                        try:
                            t_static.annotate_ncbi_taxa()
                        except Exception as e:
                            print(f"WARNING--NCBI annotation failed for static: {e}")
                    distinct_palette = [
                        "#e6194B","#3cb44b","#ffe119","#4363d8","#f58231",
                        "#911eb4","#42d4f4","#f032e6","#bfef45","#469990",
                        "#dcbeff","#9A6324","#800000","#aaffc3","#808000",
                        "#ffd8b1","#000075","#a9a9a9"
                    ]
                    unique_taxids = sorted({n.name.split(".")[0] for n in t_static.leaves()})
                    auto_colormap = {tid: distinct_palette[i % len(distinct_palette)]
                                    for i, tid in enumerate(unique_taxids)}

                def _get_static_shape(domain_name):
                    shapes = ['[]', '()', '<>', '^', 'v', 'o']
                    shape_idx = int(hashlib.md5(domain_name.encode()).hexdigest()[-2:], 16) % len(shapes)
                    return shapes[shape_idx]

                for node in t_static.traverse():
                    
                    # --- Branch colouring ---
                    if "colors" in static_layers and auto_colormap:
                        taxid = node.name.split(".")[0] if node.is_leaf else None
                        if not taxid and 'lineage' in node.props:
                            for tid in node.props['lineage'][::-1]:
                                if str(tid) in auto_colormap:
                                    taxid = str(tid)
                                    break
                        if taxid and taxid in auto_colormap:
                            nstyle = NodeStyle()
                            nstyle["vt_line_color"] = auto_colormap[taxid]
                            nstyle["hz_line_color"] = auto_colormap[taxid]
                            nstyle["fgcolor"] = auto_colormap[taxid]
                            node.set_style(nstyle)

                    if node.is_leaf:
                        name_parts = node.name.split(".")
                        accession = name_parts[1] if len(name_parts) > 1 else node.name

                        # --- Leaf name + gene label ---
                        if "names" in static_layers:
                            display = accession
                            gene = seqid2gene.get(node.name, "") if "gene" in static_layers else ""
                            node.add_face(TextFace(f"{display}  {gene}  "), column=0, position='branch-right')

                        # --- Domain shapes ---
                        if "domains" in static_layers:
                            seq = name2seq.get(node.name, None)
                            domains = domain_dict.get(accession, [])
                            if seq and domains:
                                motifs = []
                                for d in domains:
                                    motifs.append([
                                        d["ali_from"],
                                        d["ali_to"],
                                        _get_static_shape(d["hmm_name"]),
                                        None,
                                        14,
                                        "black",
                                        _domain_color(d["hmm_name"]),
                                        f"arial|8|black|{d['hmm_name']}"
                                    ])
                                node.add_face(SeqMotifFace(seq, motifs=motifs, seq_format='-'), column=1, position='aligned')

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
            elif not args.get('no_explore'):
                print(f"INFO--Starting ETE4 server on port {args['port']}")
                
                from ete4.smartview import Layout, TextFace as SmartTextFace, SeqFace, RectFace, BASIC_LAYOUT
                
                def _draw_leaf(node):
                    if not node.is_leaf:
                        return
                    name_parts = node.name.split(".")
                    display = name_parts[1] if len(name_parts) > 1 else node.name
                    gene = node.props.get("gene_name", "")
                    sci = node.props.get("sci_name", "")
                    label = f"{display}  {gene}  [{sci}]".strip(" []")
                    return SmartTextFace(label, style="fill: black;", position="right", column=0)

                leaf_name_layout = Layout(name="Leaf names", active=True, draw_node=_draw_leaf)
                
                def layout_seqface(node):
                    if node.is_leaf:
                        seq = node.props.get('seq')
                        if seq:
                            node.add_face(SeqFace(seq, seqtype='aa'), column=0, position='aligned')
                
                # Build aliases — one per unique domain name
                unique_domain_names = sorted({d["hmm_name"]
                                            for hits in domain_dict.values()
                                            for d in hits})
                name_to_alias = {n: f"dom-{i}" for i, n in enumerate(unique_domain_names)}
                domain_aliases = {
                    name_to_alias[n]: {"fill": _domain_color(n),
                                    "stroke": "black",
                                    "stroke-width": "1"}
                    for n in unique_domain_names
                }

                def _draw_domains(node):
                    if not node.is_leaf:
                        return
                    accession = node.name.split(".", 1)[1] if "." in node.name else node.name
                    domains = domain_dict.get(accession, [])
                    if not domains:
                        return
                    seq_len = max(d["ali_to"] for d in domains)
                    faces = []
                    for i, d in enumerate(domains):
                        width = max(30, int((d["ali_to"] - d["ali_from"]) / seq_len * 200))
                        faces.append(RectFace(
                            wmax=width,
                            hmax=18,
                            style=name_to_alias[d["hmm_name"]],
                            text=d["hmm_name"],
                            position="aligned",
                            column=i + 1,
                        ))
                    return faces

                domain_layout = Layout(
                    name="Domains",
                    active=True,
                    draw_tree={"aliases": domain_aliases},   
                    draw_node=_draw_domains,
                )
                
                if not args.get('no_ncbi'):
                    try:
                        t.annotate_ncbi_taxa()
                    except Exception as e:
                        print(f"WARNING--NCBI taxonomy annotation failed: {e}")

                # Attach sequences if MSA is present
                name2seq = {}
                if args.get('MSA') and os.path.isfile(filename_aln):
                    name2seq = get_seqs(filename_aln)

                # ATTACH LABELS DIRECTLY TO NODES 
                

                for node in t.traverse():
                    if colormap and 'lineage' in node.props:
                        for taxid in node.props['lineage'][::-1]:
                            if str(taxid) in colormap:
                                node.add_prop("color", colormap[str(taxid)])
                                break
                    if node.is_leaf:
                        node.add_prop("gene_name", seqid2gene.get(node.name, "-"))
                        if args.get('MSA') and node.name in name2seq:
                            node.add_prop('seq', name2seq[node.name])
                            node.add_face(SeqFace(name2seq[node.name], seqtype='aa'),
                                        column=3, position='aligned')

                t.explore(
                    layouts=[BASIC_LAYOUT, leaf_name_layout, domain_layout],
                    keep_server=True,
                    quiet=True,
                    port=args['port'],
                    host='0.0.0.0',
                    show_leaf_name=True,
                )
                            

        except Exception as e:
             print(f"ERROR--Failed running ETE4 Graphics Engine: {e}")