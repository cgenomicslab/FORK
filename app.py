import os
import sys
import uuid
import time
import base64
import socket
import signal
import threading
import subprocess
import tempfile
from pathlib import Path
from ete4 import NCBITaxa

_ncbi = NCBITaxa()

# Load .env from the same directory as this script, before any DB imports
from dotenv import load_dotenv

_HERE = Path(__file__).resolve().parent
load_dotenv(_HERE / ".env", override=False)

from flask import (
    Flask,
    render_template,
    request,
    jsonify,
    session,
    send_file,
    Response,
    stream_with_context,
)
import requests as _req

import get_reference_uniprot_set_lib as uni
import viz_utils as viz
import subclade_partition as sp
import tree_builder as tb
import recent_runs as rr  # green UI: read-only dashboard helper
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

app = Flask(__name__)
app.secret_key = os.getenv("FLASK_SECRET_KEY") or os.urandom(24)
# Pick up template edits (e.g. the logo in base.html) on a plain page reload,
# without needing to restart the app.
app.config["TEMPLATES_AUTO_RELOAD"] = True
app.jinja_env.auto_reload = True


def _find_free_port(start=5001, end=5050):
    """Return the first TCP port in [start, end) that is not already bound."""
    for port in range(start, end):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("127.0.0.1", port))
                return port
            except OSError:
                continue
    raise RuntimeError(f"No free port found in range {start}–{end}")


# Tracks every ETE4 port launched this session so /static/images/ requests
# from PixiJS web workers can be forwarded to the right server.
_ete4_ports: set = set()


# In-memory store for background jobs (tree builds, etc.)
_jobs = {}
_job_lock = threading.Lock()


def _get_config():
    cfg = uni.get_db_config(
        host=session.get("db_host") or None,
        user=session.get("db_user") or None,
        password=session.get("db_password") or None,
        database=session.get("db_name") or None,
    )
    port = session.get("db_port")
    if port:
        cfg["port"] = int(port)
    return cfg


def _new_job(job_type, meta=None):
    job_id = uuid.uuid4().hex[:12]
    with _job_lock:
        _jobs[job_id] = {
            "type": job_type,
            "status": "running",
            "log": [],
            "result": None,
            "error": None,
            "created": time.time(),   # green UI: for the Recent runs dashboard
            "meta": meta or {},       # green UI: family / method / output path
        }
    return job_id


# ==============================================================================================================
#                                           PAGE ROUTES
# ==============================================================================================================


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/tree")
def tree_page():
    return render_template("tree.html")


@app.route("/presence")
def presence_page():
    return render_template("presence.html")


@app.route("/highres")
def highres_page():
    return render_template("highres.html")


@app.route("/profiling")
def profiling_page():
    return render_template("profiling.html")


@app.route("/utilities")
def utilities_page():
    return render_template("utilities.html")


@app.route("/about")
def about_page():
    return render_template("about.html")


# ==============================================================================================================
#                                           DB CONFIG
# ==============================================================================================================


@app.route("/api/db-config", methods=["POST"])
def api_db_config():
    data = request.json or {}
    session["db_host"] = data.get("host") or ""
    session["db_user"] = data.get("user") or ""
    session["db_password"] = data.get("password") or ""
    session["db_name"] = data.get("database") or ""
    session["db_port"] = data.get("port") or ""
    return jsonify({"ok": True})


@app.route("/api/db-defaults")
def api_db_defaults():
    """Return current env-based DB defaults and seed Flask session so connection works immediately."""
    host = os.getenv("DB_HOST", "localhost")
    user = os.getenv("DB_USER", "")
    database = os.getenv("DB_NAME", "")
    port = os.getenv("DB_PORT", "3306")
    password = os.getenv("DB_PASSWORD", "")

    # Seed session from env vars on first load so _get_config() works without
    # requiring the user to manually open and submit the DB panel.
    if not session.get("db_host"):
        session["db_host"] = host
        session["db_user"] = user
        session["db_password"] = password
        session["db_name"] = database
        session["db_port"] = port

    return jsonify(
        {
            "host": host,
            "user": user,
            "database": database,
            "port": port,
        }
    )


@app.route("/api/db-info")
def api_db_info():
    try:
        config = _get_config()
        with uni.UniProtRetriever(config) as db:
            versions = db.list_available_versions()
        return jsonify({"versions": versions})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/recent-runs")
def api_recent_runs():
    """Recent activity for the Overview dashboard.

    Reads the job registry and mirrors it into a small history file so runs
    persist across restarts. Does not touch any analysis output.
    """
    try:
        with _job_lock:
            runs = rr.collect(_jobs)
        return jsonify({"runs": runs})
    except Exception as e:
        return jsonify({"runs": [], "error": str(e)})


# ==============================================================================================================
#                                           GENERIC JOB POLLING
# ==============================================================================================================


@app.route("/api/job/<job_id>")
def api_job(job_id):
    with _job_lock:
        job = _jobs.get(job_id)
    if not job:
        return jsonify({"error": "Not found"}), 404

    out = {
        "status": job["status"],
        "log": job["log"][-20:],
        "error": job["error"],
    }

    if job["status"] == "done" and job["result"]:
        r = job["result"]
        out["viewer"] = r.get("viewer", "")
        out["prefix"] = r.get("prefix", "")
        if r.get("viewer") == "ete4":
            out["port"] = r.get("port")
        elif r.get("viewer") == "ete4_static":
            out["img_available"] = os.path.isfile(r.get("img_path", ""))
        elif r.get("viewer") == "d3":
            out["has_newick"] = bool(r.get("newick"))
            out["tree_base"] = r.get("tree_base", "")
            out["has_itol_colors"] = bool(r.get("itol_colors"))
            out["has_itol_domains"] = bool(r.get("itol_domains"))

    return jsonify(out)


# ==============================================================================================================
#                                           PHYLOGENETIC TREE
# ==============================================================================================================


@app.route("/api/run-tree", methods=["POST"])
def api_run_tree():
    f = request.form
    files = request.files

    pfam = f.get("pfam", "").strip()
    ver = f.get("version", "2026_01")
    prefix = f.get("prefix", "").strip()
    output_dir = f.get("output_dir", "").strip()
    aln = f.get("aln", "mafft")
    ml = f.get("ml", "fasttree")
    trimal_th = f.get("trimal_th", "0.01")
    cpu = f.get("cpu", "4")
    evalue = f.get("evalue", "").strip()
    no_ncbi = f.get("no_ncbi") == "true"
    pfam_source = f.get("pfam_source", "hmmsearch")
    pfam_logic = f.get("pfam_logic", "or")
    color_by = f.get("color_by", "taxon")
    msa_range = f.get("msa_range", "").strip()
    local_fasta = f.get("local_fasta", "").strip()
    viewer = f.get("viewer", "d3")
    attach_msa = f.get("attach_msa") == "true"
    static_layers = f.get("static_layers", "names,domains,colors,gene")

    if not pfam or not prefix:
        return jsonify({"error": "Pfam ID and output prefix are required."}), 400

    # Resolve taxids from text or uploaded file
    tax_raw = f.get("taxids", "").strip()
    if "tax_file" in files and files["tax_file"].filename:
        tax_raw = files["tax_file"].read().decode()
    tax_str = ",".join(
        l.strip() for l in tax_raw.replace(",", "\n").splitlines() if l.strip()
    )

    excl_raw = f.get("exclude_taxids", "").strip()
    if "excl_file" in files and files["excl_file"].filename:
        excl_raw = files["excl_file"].read().decode()
    excl_str = ",".join(
        l.strip() for l in excl_raw.replace(",", "\n").splitlines() if l.strip()
    )

    colormap_path = None
    if "colormap" in files and files["colormap"].filename:
        tmp = tempfile.NamedTemporaryFile(
            mode="wb", suffix=".colormap.txt", delete=False
        )
        tmp.write(files["colormap"].read())
        tmp.close()
        colormap_path = tmp.name

    job_id = _new_job(
        "tree",
        meta={
            "pfam": pfam,
            "method": f"{aln} → {ml}",
            "output": output_dir or prefix,
        },
    )

    def run():
        try:
            cmd = [
                sys.executable,
                "tree_from_db.py",
                "--pfam",
                pfam.replace(" ", ""),
                "--version",
                ver,
                "--prefix",
                prefix,
                "--aln",
                aln,
                "--ml",
                ml,
                "--cpu",
                cpu,
                "--gt",
                trimal_th,
            ]
            if output_dir:
                cmd += ["--output_dir", output_dir]
            if tax_str:
                cmd += ["--taxids", tax_str]
            if excl_str:
                cmd += ["--exclude_taxids", excl_str]
            if evalue:
                cmd += ["--evalue", evalue]
            if no_ncbi:
                cmd += ["--no_ncbi"]
            cmd += [
                "--pfam_source",
                pfam_source,
                "--pfam_logic",
                pfam_logic,
                "--color_by",
                color_by,
            ]
            if local_fasta:
                cmd += ["--local_fasta", local_fasta]
            if msa_range:
                cmd += ["--positions", msa_range]
            if colormap_path:
                cmd += ["--colormap", colormap_path]

            env = os.environ.copy()
            env["QT_QPA_PLATFORM"] = "offscreen"

            if viewer == "ete4":
                ete4_port = _find_free_port()
                _ete4_ports.add(ete4_port)
                cmd += ["--port", str(ete4_port)]
                if attach_msa:
                    cmd += ["--MSA"]
                os.system(f"fuser -k {ete4_port}/tcp >/dev/null 2>&1")
                proc = subprocess.Popen(
                    cmd,
                    env=env,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )

                # Stream the build's progress into the job log while it runs, so
                # the tree page shows live messages (fetching, aligning, building
                # tree, starting the ETE4 server…). Runs in a thread because the
                # ETE4 server keeps printing after the tree is built.
                def _pump_log(p=proc):
                    try:
                        for line in p.stdout:
                            with _job_lock:
                                _jobs[job_id]["log"].append(line.rstrip())
                    except Exception:
                        pass

                threading.Thread(target=_pump_log, daemon=True).start()

                deadline = time.time() + 600
                while time.time() < deadline:
                    try:
                        with socket.create_connection(
                            ("localhost", ete4_port), timeout=2
                        ):
                            break
                    except OSError:
                        time.sleep(3)
                else:
                    with _job_lock:
                        _jobs[job_id].update(
                            {
                                "status": "error",
                                "error": "ETE4 server did not start within 10 minutes.",
                            }
                        )
                    return

                with _job_lock:
                    _jobs[job_id].update(
                        {
                            "status": "done",
                            "result": {
                                "viewer": "ete4",
                                "port": ete4_port,
                                "prefix": prefix,
                            },
                        }
                    )

            elif viewer == "ete4_static":
                ete4_port = _find_free_port()
                cmd += ["--port", str(ete4_port), "--render_ete_static", "--no_explore"]
                if static_layers:
                    cmd += ["--static_layers", static_layers]

                # Stream progress into the job log so the tree page shows live
                # build messages while the static image is being rendered.
                proc = subprocess.Popen(
                    cmd,
                    env=env,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                for line in proc.stdout:
                    with _job_lock:
                        _jobs[job_id]["log"].append(line.rstrip())
                proc.wait()
                if proc.returncode != 0:
                    with _job_lock:
                        _jobs[job_id].update(
                            {
                                "status": "error",
                                "error": f"Pipeline exited with code {proc.returncode} (see log).",
                            }
                        )
                    return

                img_path = f"{prefix}_tree_domains.png"
                if output_dir:
                    img_path = os.path.join(output_dir, img_path)

                with _job_lock:
                    _jobs[job_id].update(
                        {
                            "status": "done",
                            "result": {
                                "viewer": "ete4_static",
                                "img_path": img_path,
                                "prefix": prefix,
                            },
                        }
                    )

            else:
                cmd += ["--no_explore"]
                proc = subprocess.Popen(
                    cmd,
                    env=env,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )

                for line in proc.stdout:
                    with _job_lock:
                        _jobs[job_id]["log"].append(line.rstrip())

                proc.wait()
                if proc.returncode != 0:
                    with _job_lock:
                        _jobs[job_id].update(
                            {
                                "status": "error",
                                "error": f"Pipeline exited with code {proc.returncode}.",
                            }
                        )
                    return

                aln_ext = ".mft" if aln == "mafft" else f".{aln}"
                trim_ext = f'.gt{trimal_th.replace(".", "")}'
                tree_base = (
                    f"{prefix}{aln_ext}{trim_ext}.lg.fasttree"
                    if ml == "fasttree"
                    else f"{prefix}{aln_ext}{trim_ext}.treefile"
                )
                tree_path = (
                    os.path.join(output_dir, tree_base) if output_dir else tree_base
                )

                if not os.path.isfile(tree_path):
                    with _job_lock:
                        _jobs[job_id].update(
                            {
                                "status": "error",
                                "error": f"Expected output not found: {tree_path}",
                            }
                        )
                    return

                newick = Path(tree_path).read_text()

                itol_colors = None
                cp = tree_path + ".itol_colors.txt"
                if os.path.isfile(cp):
                    itol_colors = Path(cp).read_text()

                itol_domains = None
                dp = tree_path + ".itol_domains.txt"
                if os.path.isfile(dp):
                    itol_domains = Path(dp).read_text()

                with _job_lock:
                    _jobs[job_id].update(
                        {
                            "status": "done",
                            "result": {
                                "viewer": "d3",
                                "newick": newick,
                                "tree_path": tree_path,
                                "tree_base": tree_base,
                                "itol_colors": itol_colors,
                                "itol_domains": itol_domains,
                                "prefix": prefix,
                            },
                        }
                    )

        except Exception as e:
            with _job_lock:
                _jobs[job_id].update({"status": "error", "error": str(e)})

    threading.Thread(target=run, daemon=True).start()
    return jsonify({"job_id": job_id})


@app.route("/api/tree-static-png/<job_id>")
def api_tree_static_png(job_id):
    with _job_lock:
        job = _jobs.get(job_id)
    if not job or job["status"] != "done":
        return jsonify({"error": "Not ready"}), 404
    r = job["result"]

    if r.get("viewer") == "d3":
        buf = viz.render_tree(r["newick"])
        b64 = base64.b64encode(buf.read()).decode()
        return jsonify({"png": b64})

    if r.get("viewer") == "ete4_static":
        p = r.get("img_path", "")
        if os.path.isfile(p):
            b64 = base64.b64encode(Path(p).read_bytes()).decode()
            return jsonify({"png": b64})

    return jsonify({"error": "No image available"}), 404


@app.route("/api/download/tree/<job_id>/<filetype>")
def api_download_tree(job_id, filetype):
    with _job_lock:
        job = _jobs.get(job_id)
    if not job or job["status"] != "done":
        return "Not ready", 404
    r = job["result"]
    prefix = r.get("prefix", "tree")

    if filetype == "newick":
        return Response(
            r.get("newick", ""),
            mimetype="text/plain",
            headers={"Content-Disposition": f"attachment; filename={prefix}.nwk"},
        )
    if filetype == "itol_colors":
        return Response(
            r.get("itol_colors", ""),
            mimetype="text/plain",
            headers={
                "Content-Disposition": f"attachment; filename={prefix}_itol_colors.txt"
            },
        )
    if filetype == "itol_domains":
        return Response(
            r.get("itol_domains", ""),
            mimetype="text/plain",
            headers={
                "Content-Disposition": f"attachment; filename={prefix}_itol_domains.txt"
            },
        )
    if filetype == "png":
        p = r.get("img_path", "")
        if p and os.path.isfile(p):
            return send_file(
                p,
                mimetype="image/png",
                as_attachment=True,
                download_name=os.path.basename(p),
            )

    return "Not found", 404


# ==============================================================================================================
#                                           PRESENCE / ABSENCE
# ==============================================================================================================


@app.route("/api/presence-matrix", methods=["POST"])
def api_presence_matrix():
    data = request.json or {}
    ver = data.get("version", "2026_01")
    pfam_queries = data.get("pfam_queries", [])
    tax_ids = data.get("tax_ids") or None
    evalue = data.get("evalue") or None

    if not pfam_queries:
        return jsonify({"error": "At least one Pfam query is required."}), 400

    try:
        import pandas as pd

        config = _get_config()
        rows = uni.fetch_presence_absence_matrix(
            ver, pfam_queries, tax_ids, evalue, db_config=config
        )

        if not rows:
            return jsonify(
                {
                    "rows": [],
                    "columns": [],
                    "matrix": [],
                    "heatmap": None,
                    "taxon_map": {},
                    "profiles": [],
                }
            )

        df = pd.DataFrame(rows)
        df["taxon_label"] = df.apply(
            lambda r: (
                f"{r['taxon_id']} · {r['scientific_name']}"
                if r.get("scientific_name")
                else str(r["taxon_id"])
            ),
            axis=1,
        )
        matrix = df.pivot_table(
            index="taxon_label",
            columns="hmm_name",
            values="protein_count",
            aggfunc="sum",
            fill_value=0,
        )
        matrix.index.name = "Organism"
        matrix.columns.name = None

        heatmap_b64 = None
        if matrix.shape[0] >= 2:
            buf = viz.draw_presence_absence_heatmap(
                matrix,
                title="Presence / Absence",
                cluster=(matrix.shape[1] >= 2),
                cmap="viridis",
            )
            buf.seek(0)
            heatmap_b64 = base64.b64encode(buf.read()).decode()

        taxon_map = (
            df[["taxon_label", "taxon_id"]]
            .drop_duplicates()
            .set_index("taxon_label")["taxon_id"]
            .to_dict()
        )

        return jsonify(
            {
                "columns": ["Organism"] + list(matrix.columns),
                "matrix": matrix.reset_index().to_dict(orient="records"),
                "heatmap": heatmap_b64,
                "taxon_map": {str(k): int(v) for k, v in taxon_map.items()},
                "profiles": list(matrix.columns),
                "n_rows": len(matrix),
            }
        )

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/drill-down", methods=["POST"])
def api_drill_down():
    data = request.json or {}
    ver = data.get("version", "2026_01")
    profile = data.get("profile", "")
    taxon_id = data.get("taxon_id")
    evalue = data.get("evalue") or None

    try:
        config = _get_config()
        cell_records = uni.fetch_accessions_for_cell(
            ver, profile, taxon_id, evalue, db_config=config
        )
        accessions = [r["accession"] for r in cell_records]

        subprofiles = uni.fetch_subprofile_hits(
            ver, accessions, evalue, exclude_queries=[profile], db_config=config
        )
        archs = uni.fetch_domain_architectures(
            ver, accessions, evalue, collapse_repeats=True, db_config=config
        )

        return jsonify(
            {
                "n_proteins": len(accessions),
                "subprofiles": subprofiles,
                "architectures": archs,
            }
        )

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/draw-architecture", methods=["POST"])
def api_draw_architecture():
    data = request.json or {}
    ver = data.get("version", "2026_01")
    accessions = data.get("accessions", [])
    evalue = data.get("evalue") or None
    title = data.get("title", "Domain Architecture")

    try:
        config = _get_config()
        domains = uni.fetch_domains_by_accession(
            ver, accessions, evalue, db_config=config
        )
        if not domains:
            return jsonify({"error": "No domain data found."}), 404

        buf = viz.draw_domain_architecture(domains, title=title)
        buf.seek(0)
        return jsonify({"png": base64.b64encode(buf.read()).decode()})

    except Exception as e:
        return jsonify({"error": str(e)}), 500


# ==============================================================================================================
#                                           HIGH RESOLUTION PROFILE
# ==============================================================================================================


@app.route("/api/highres/build-trees", methods=["POST"])
def api_highres_build_trees():
    f = request.form
    files = request.files

    pfam_text = f.get("pfams", "")
    ver = f.get("version", "2026_01")
    output_root = f.get("output_root", "/tmp/highres_runs").strip()
    aln = f.get("aln", "mafft")
    ml = f.get("ml", "fasttree")
    gt = f.get("gt", "0.01")
    cpu = int(f.get("cpu", 8))
    evalue = f.get("evalue", "").strip() or None

    tax_raw = f.get("taxids", "")
    if "tax_file" in files and files["tax_file"].filename:
        tax_raw += "\n" + files["tax_file"].read().decode()

    excl_raw = f.get("exclude_taxids", "")
    if "excl_file" in files and files["excl_file"].filename:
        excl_raw += "\n" + files["excl_file"].read().decode()

    pfams, seen = [], set()
    for tok in pfam_text.replace(",", " ").split():
        if tok and tok not in seen:
            seen.add(tok)
            pfams.append(tok)

    def parse_ids(raw):
        ids = []
        for line in raw.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            for tok in line.replace(",", " ").split():
                try:
                    ids.append(int(tok))
                except ValueError:
                    pass
        return sorted(set(ids)) or None

    taxids = parse_ids(tax_raw)
    exclude_taxids = parse_ids(excl_raw)

    if not pfams:
        return jsonify({"error": "At least one Pfam ID is required."}), 400

    job_id = _new_job(
        "highres_trees",
        meta={
            "pfam": ", ".join(pfams),
            "method": f"{aln} → {ml}",
            "output": output_root,
        },
    )

    def run():
        log = []

        def cb(msg):
            log.append(msg)
            with _job_lock:
                _jobs[job_id]["log"] = log[-25:]

        try:
            results = tb.build_trees(
                pfams=pfams,
                output_root=output_root,
                version=ver,
                taxids=taxids,
                exclude_taxids=exclude_taxids,
                evalue=float(evalue) if evalue else None,
                aln=aln,
                ml=ml,
                gt=gt,
                cpu=cpu,
                tree_from_db_path="tree_from_db.py",
                progress_callback=cb,
            )

            summary = {}
            for pfam, r in results.items():
                summary[pfam] = {
                    "error": r["error"],
                    "cached": r["cached"],
                    "n_leaves": len(r["leaves"]),
                    "leaves": r["leaves"],
                    "tree_path": r["tree_path"] or "",
                    "prefix": r.get("prefix", ""),
                    "cache_key": r.get("cache_key", ""),
                    "stderr": (r.get("stderr") or "")[:2000],
                    "max_depth": (
                        sp.get_max_root_distance(r["tree"]) if r["tree"] else 0.0
                    ),
                }

            with _job_lock:
                _jobs[job_id].update(
                    {
                        "status": "done",
                        "result": {
                            "summary": summary,
                            "tree_objects": results,
                            "version": ver,
                            "taxids": taxids,
                            "output_root": output_root,
                        },
                    }
                )

        except Exception as e:
            with _job_lock:
                _jobs[job_id].update({"status": "error", "error": str(e)})

    threading.Thread(target=run, daemon=True).start()
    return jsonify({"job_id": job_id})


@app.route("/api/highres/job/<job_id>")
def api_highres_job(job_id):
    with _job_lock:
        job = _jobs.get(job_id)
    if not job:
        return jsonify({"error": "Not found"}), 404

    out = {"status": job["status"], "log": job["log"][-25:], "error": job["error"]}
    if job["status"] == "done" and job["result"]:
        out["summary"] = job["result"]["summary"]
    return jsonify(out)


@app.route("/api/highres/partition", methods=["POST"])
def api_highres_partition():
    data = request.json or {}
    job_id = data.get("job_id")
    pfam = data.get("pfam")
    mode = data.get("mode", "depth")
    taxon_level = data.get("taxon_level")

    with _job_lock:
        job = _jobs.get(job_id)
    if not job or job["status"] != "done":
        return jsonify({"error": "Job not ready."}), 404

    tree_objects = job["result"].get("tree_objects", {})
    r = tree_objects.get(pfam, {})
    if not r or r.get("error") or not r.get("tree"):
        return jsonify({"error": f"No valid tree for {pfam}."}), 400

    tree = r["tree"]

    try:
        if mode == "depth":
            parts = sp.partition_by_depth(tree, float(data.get("threshold", 0.5)))
        elif mode == "mrca":
            groups = data.get("groups", [])
            parts = sp.partition_by_mrca(
                tree, groups, include_unassigned=data.get("include_unassigned", True)
            )
        elif mode == "node_path":
            paths = [[int(x) for x in p] for p in data.get("paths", [])]
            parts = sp.partition_by_node_path(tree, paths)
        elif (
            mode == "auto_duplication"
        ):  # auto split base on duplication events, by taxonomic level
            parts = sp.partition_by_duplication(tree, taxon_level, NCBITaxa())
        else:
            return jsonify({"error": "Unknown mode."}), 400

        parts_json = {label: sorted(members) for label, members in parts.items()}
        summary = [
            {"label": l, "n_leaves": len(m), "sample": sorted(m)[:3]}
            for l, m in parts.items()
        ]
        return jsonify({"parts": parts_json, "summary": summary})

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/highres/list-nodes", methods=["POST"])
def api_highres_list_nodes():
    data = request.json or {}
    with _job_lock:
        job = _jobs.get(data.get("job_id"))
    if not job or job["status"] != "done":
        return jsonify({"error": "Not ready"}), 404

    pfam = data.get("pfam")
    tree = job["result"]["tree_objects"].get(pfam, {}).get("tree")
    if not tree:
        return jsonify({"error": f"No tree for {pfam}"}), 400

    return jsonify({"nodes": sp.list_internal_nodes(tree)})


@app.route("/api/highres/launch-ete4", methods=["POST"])
def api_highres_launch_ete4():
    data = request.json or {}
    job_id = data.get("job_id")
    pfam = data.get("pfam")

    with _job_lock:
        job = _jobs.get(job_id)
    if not job or job["status"] != "done":
        return jsonify({"error": "Not ready"}), 404

    r = job["result"]["tree_objects"].get(pfam, {})
    ver = job["result"].get("version", "2026_01")

    try:
        port = _find_free_port()
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 500

    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    cmd = [
        sys.executable,
        "tree_from_db.py",
        "--pfam",
        pfam,
        "--version",
        ver,
        "--prefix",
        r["prefix"],
        "--port",
        str(port),
        "--use_resolved",
    ]
    proc = subprocess.Popen(cmd, env=env)

    deadline = time.time() + 120
    connected = False
    while time.time() < deadline:
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=2):
                connected = True
                break
        except OSError:
            time.sleep(2)

    if not connected:
        return jsonify({"error": "ETE4 server did not start in time."}), 500

    _ete4_ports.add(port)
    ete_key = f"ete4_{job_id}_{pfam}"
    with _job_lock:
        _jobs[ete_key] = {"pid": proc.pid, "port": port}

    return jsonify({"ok": True, "pid": proc.pid, "port": port, "ete_key": ete_key})


@app.route("/api/highres/stop-ete4", methods=["POST"])
def api_highres_stop_ete4():
    data = request.json or {}
    ete_key = data.get("ete_key")
    with _job_lock:
        ete = _jobs.pop(ete_key, None)
    if ete:
        try:
            os.kill(ete["pid"], signal.SIGKILL)
        except ProcessLookupError:
            pass
        os.system(f"fuser -k {ete['port']}/tcp >/dev/null 2>&1")
    return jsonify({"ok": True})


@app.route("/api/highres/compute-profile", methods=["POST"])
def api_highres_compute_profile():
    data = request.json or {}
    job_id = data.get("job_id")
    partitions = data.get("partitions", {})  # {pfam: {label: [leaf_names]}}
    binary = data.get("binary", False)
    log_scale = data.get("log_scale", False)
    cluster_cols = data.get("cluster_cols", False)

    with _job_lock:
        job = _jobs.get(job_id)
    if not job or job["status"] != "done":
        return jsonify({"error": "Job not ready."}), 404

    ver = job["result"].get("version", "2026_01")
    taxids = job["result"].get("taxids")

    # Strip taxid prefix from leaf names to get bare accessions
    pfam_subclade_map = {
        pfam: tb.strip_leaf_prefix_in_subclades(
            {label: set(leaves) for label, leaves in subclades.items()}
        )
        for pfam, subclades in partitions.items()
    }

    try:
        config = _get_config()
        out = uni.fetch_highres_profile(
            version=ver,
            pfam_subclade_map=pfam_subclade_map,
            taxon_ids=taxids or None,
            binary=binary,
            db_config=config,
        )

        matrix = out["matrix"]
        buf = viz.draw_highres_profile_heatmap(
            matrix,
            column_origin=out["column_origin"],
            taxon_names=out["taxon_names"],
            missing_accessions=out["missing_accessions"],
            binary=binary,
            log_scale=log_scale,
            cluster_cols=cluster_cols,
        )
        buf.seek(0)
        heatmap_b64 = base64.b64encode(buf.read()).decode()

        tn = out["taxon_names"]
        display = matrix.copy()
        display.index = [f"{tx}  {tn.get(tx, '')}" for tx in display.index]

        return jsonify(
            {
                "heatmap": heatmap_b64,
                "matrix": display.reset_index().to_dict(orient="records"),
                "matrix_csv": matrix.to_csv(),
                "missing": sorted(out["missing_accessions"]),
                "n_missing": len(out["missing_accessions"]),
            }
        )

    except Exception as e:
        return jsonify({"error": str(e)}), 500


# ==============================================================================================================
#                                           UTILITIES
# ==============================================================================================================


@app.route("/api/fetch-sequences", methods=["POST"])
def api_fetch_sequences():
    data = request.json or {}
    ver = data.get("version", "2026_01")
    tax_ids = data.get("tax_ids") or None
    proteome = data.get("proteome") or None
    go_id = data.get("go_id") or None
    pfam_id = data.get("pfam_id") or None

    try:
        config = _get_config()
        records = uni.fetch_sequences(
            ver, tax_ids, proteome, go_id, pfam_id, db_config=config
        )
        if not records:
            return jsonify({"records": [], "fasta": "", "count": 0})

        fasta = uni.fetch_fasta_string(
            ver, tax_ids, proteome, go_id, pfam_id, db_config=config
        )
        display = [{k: v for k, v in r.items() if k != "sequence"} for r in records]
        return jsonify({"records": display, "fasta": fasta, "count": len(records)})

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/hmm-search", methods=["POST"])
def api_hmm_search():
    data = request.json or {}
    ver = data.get("version", "2026_01")
    query = data.get("hmm_query", "")
    evalue = data.get("evalue") or None
    tax_ids = data.get("tax_ids") or None

    try:
        config = _get_config()
        records = uni.fetch_sequences_by_hmm_hit(
            ver, query, evalue, tax_ids, db_config=config
        )
        if not records:
            return jsonify({"records": [], "fasta": "", "count": 0})

        fasta = uni.fetch_fasta_string_by_hmm_hit(
            ver, query, evalue, tax_ids, db_config=config
        )
        display = [{k: v for k, v in r.items() if k != "sequence"} for r in records]
        return jsonify({"records": display, "fasta": fasta, "count": len(records)})

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/accession-lookup", methods=["POST"])
def api_accession_lookup():
    data = request.json or {}
    ver = data.get("version", "2026_01")
    accs = data.get("accessions", [])

    try:
        config = _get_config()
        records = uni.fetch_sequences_by_accession(ver, accs, db_config=config)
        if not records:
            return jsonify({"records": [], "fasta": "", "count": 0})

        with uni.UniProtRetriever(config) as db:
            fasta = db.to_fasta_string(records)

        return jsonify({"records": records, "fasta": fasta, "count": len(records)})

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/domain-lookup", methods=["POST"])
def api_domain_lookup():
    data = request.json or {}
    ver = data.get("version", "2026_01")
    accs = data.get("accessions", [])
    evalue = data.get("evalue") or None
    draw = data.get("draw", False)

    try:
        config = _get_config()
        domains = uni.fetch_domains_by_accession(ver, accs, evalue, db_config=config)
        if not domains:
            return jsonify({"domains": [], "count": 0})

        result = {"domains": domains, "count": len(domains)}
        if draw:
            title = f"Domain Architecture — {', '.join(accs[:3])}" + (
                " …" if len(accs) > 3 else ""
            )
            buf = viz.draw_domain_architecture(domains, title=title)
            buf.seek(0)
            result["png"] = base64.b64encode(buf.read()).decode()

        return jsonify(result)

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/go-domains", methods=["POST"])
def api_go_domains():
    data = request.json or {}
    ver = data.get("version", "2026_01")
    go_term = data.get("go_term", "")
    evalue = data.get("evalue") or None

    try:
        config = _get_config()
        results = uni.fetch_domains_by_go(ver, go_term, evalue, db_config=config)
        return jsonify({"results": results or [], "count": len(results or [])})

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/extract-branch", methods=["POST"])
def api_extract_branch():
    ver = request.form.get("version", "2026_01")
    if "branch_file" not in request.files:
        return jsonify({"error": "No file uploaded."}), 400

    from ete4 import PhyloTree

    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".nwk") as tmp:
            tmp.write(request.files["branch_file"].read())
            tmp_path = tmp.name

        t = PhyloTree(tmp_path)
        accs = []
        for leaf in t.leaves():
            parts = leaf.name.split(".")
            accs.append(parts[1] if len(parts) > 1 else leaf.name)

        if not accs:
            return jsonify({"error": "No leaves found in branch file."}), 400

        config = _get_config()
        records = uni.fetch_sequences_by_accession(ver, accs, db_config=config)
        fasta = ""
        if records:
            with uni.UniProtRetriever(config) as db:
                fasta = db.to_fasta_string(records)

        return jsonify(
            {
                "accessions": sorted(accs),
                "n_accessions": len(accs),
                "n_sequences": len(records),
                "fasta": fasta,
            }
        )

    except Exception as e:
        return jsonify({"error": str(e)}), 500
    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)


# ---------------------------------------------------------------------------
# High-Res: FASTA download for all locked subclades
# ---------------------------------------------------------------------------


@app.route("/api/highres/download-fasta", methods=["POST"])
def api_highres_download_fasta():
    data = request.json or {}
    job_id = data.get("job_id")
    partitions = data.get("partitions", {})  # {pfam: {label: [leaf_names]}}

    with _job_lock:
        job = _jobs.get(job_id)
    if not job or job["status"] != "done":
        return jsonify({"error": "Job not ready."}), 404

    ver = job["result"].get("version", "2026_01")

    pfam_subclade_map = {
        pfam: tb.strip_leaf_prefix_in_subclades(
            {label: set(leaves) for label, leaves in subclades.items()}
        )
        for pfam, subclades in partitions.items()
    }

    try:
        config = _get_config()
        fasta_parts = []

        for pfam, subclades in pfam_subclade_map.items():
            for label, accessions in subclades.items():
                accs = sorted(accessions)
                if not accs:
                    continue
                records = uni.fetch_sequences_by_accession(ver, accs, db_config=config)
                if not records:
                    continue
                with uni.UniProtRetriever(config) as db:
                    fasta_block = db.to_fasta_string(records)
                # Annotate each header with pfam + subclade label
                annotated = []
                for line in fasta_block.splitlines():
                    if line.startswith(">"):
                        annotated.append(f"{line} [pfam={pfam}|subclade={label}]")
                    else:
                        annotated.append(line)
                fasta_parts.append("\n".join(annotated))

        fasta_text = "\n".join(fasta_parts) + "\n"
        return Response(
            fasta_text,
            mimetype="text/plain",
            headers={
                "Content-Disposition": "attachment; filename=highres_subclades.fasta"
            },
        )
    except Exception as e:
        return jsonify({"error": str(e)}), 500


# ---------------------------------------------------------------------------
# High-Res: NCBI static species tree (re-implements Streamlit "2b" section)
# ---------------------------------------------------------------------------


@app.route("/api/highres/species-tree", methods=["POST"])
def api_highres_species_tree():
    job_id = request.form.get("job_id", "")
    files = request.files

    with _job_lock:
        job = _jobs.get(job_id) if job_id else None

    # Collect taxids from gene-tree leaves
    sp_taxids: list[str] = []
    if job and job["status"] == "done":
        for r in job["result"].get("tree_objects", {}).values():
            for leaf in r.get("leaves", []):
                t = leaf.split(".")[0]
                if t.isdigit():
                    sp_taxids.append(t)
    sp_taxids = sorted(set(sp_taxids))

    if len(sp_taxids) < 2:
        return (
            jsonify(
                {
                    "error": "Need at least 2 valid taxids (run tree-build first, or the trees had no leaves)."
                }
            ),
            400,
        )

    # Optional colormap file
    sp_colormap: dict[str, str] = {}
    if "colormap" in files and files["colormap"].filename:
        for line in files["colormap"].read().decode().splitlines():
            parts = line.split()
            if len(parts) >= 2:
                sp_colormap[parts[0]] = parts[1]

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    try:
        from ete4 import NCBITaxa
        from ete4.treeview import TreeStyle, NodeStyle, TextFace as TVTextFace

        ncbi = NCBITaxa()
        sp_tree = ncbi.get_topology(sp_taxids, intermediate_nodes=False)
        for n in sp_tree.traverse():
            n.dist = 0.0
        sp_tree.to_ultrametric(topological=True)

        _palette = [
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
        leaf_taxids = sorted({l.name for l in sp_tree.leaves()})
        auto_cmap = {
            tid: _palette[i % len(_palette)] for i, tid in enumerate(leaf_taxids)
        }
        auto_cmap.update(sp_colormap)

        sp_tree.annotate_ncbi_taxa(taxid_attr="name")

        for node in sp_tree.traverse():
            ns = NodeStyle()
            ns["hz_line_width"] = 4
            ns["vt_line_width"] = 4
            ns["size"] = 0
            if node.is_leaf:
                col = auto_cmap.get(node.name)
            else:
                leaf_cols = {auto_cmap.get(l.name) for l in node.leaves()}
                col = next(iter(leaf_cols)) if len(leaf_cols) == 1 else None
            if col:
                ns["hz_line_color"] = col
                ns["vt_line_color"] = col
            node.set_style(ns)
            if node.is_leaf:
                label = node.props.get("sci_name", node.name)
                node.add_face(
                    TVTextFace(f"  {label} ({node.name})"),
                    column=0,
                    position="branch-right",
                )

        ts = TreeStyle()
        ts.show_leaf_name = False

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = tmp.name
        sp_tree.render(tmp_path, w=2000, units="px", tree_style=ts)
        b64 = base64.b64encode(Path(tmp_path).read_bytes()).decode()
        os.unlink(tmp_path)

        return jsonify({"png": b64, "n_taxa": len(sp_taxids)})

    except Exception as e:
        return jsonify({"error": str(e)}), 500


# ---------------------------------------------------------------------------
# High-Res: ETE4 interactive view of the computed profile matrix on NCBI tree
# ---------------------------------------------------------------------------


@app.route("/api/highres/ete-matrix-viz", methods=["POST"])
def api_highres_ete_matrix_viz():
    matrix_csv = request.form.get("matrix_csv", "")
    files = request.files

    if not matrix_csv:
        return jsonify({"error": "No matrix data provided."}), 400

    try:
        port = _find_free_port()
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 500

    mat_tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False)
    mat_tmp.write(matrix_csv)
    mat_tmp.close()

    cmap_path = None
    if "colormap" in files and files["colormap"].filename:
        cmap_tmp = tempfile.NamedTemporaryFile(
            mode="wb", suffix=".colormap.txt", delete=False
        )
        cmap_tmp.write(files["colormap"].read())
        cmap_tmp.close()
        cmap_path = cmap_tmp.name

    cmd = [
        sys.executable,
        "ete_highres_profile.py",
        "-m",
        mat_tmp.name,
        "-p",
        str(port),
    ]
    if cmap_path:
        cmd += ["-c", cmap_path]

    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    os.system(f"fuser -k {port}/tcp >/dev/null 2>&1")
    proc = subprocess.Popen(cmd, env=env)

    deadline = time.time() + 120
    connected = False
    while time.time() < deadline:
        try:
            with socket.create_connection(("localhost", port), timeout=2):
                connected = True
                break
        except OSError:
            time.sleep(3)

    if not connected:
        proc.kill()
        return jsonify({"error": "ETE profile viewer did not start in time."}), 500

    _ete4_ports.add(port)
    return jsonify({"ok": True, "port": port, "pid": proc.pid})


# ---------------------------------------------------------------------------
# ETE Species Tree Profile — launch ete_profile.py as a subprocess
# ---------------------------------------------------------------------------


@app.route("/api/ete-profile/launch", methods=["POST"])
def api_ete_profile_launch():
    files = request.files
    form = request.form

    if "tblout" not in files or not files["tblout"].filename:
        return jsonify({"error": "tblout file is required."}), 400
    if "colormap" not in files or not files["colormap"].filename:
        return jsonify({"error": "colormap file is required."}), 400
    if "taxids" not in files or not files["taxids"].filename:
        return jsonify({"error": "taxids file is required."}), 400

    max_val = form.get("max", "").strip()

    try:
        port = _find_free_port()
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 500

    tblout_tmp = tempfile.NamedTemporaryFile(mode="wb", suffix=".tblout", delete=False)
    colormap_tmp = tempfile.NamedTemporaryFile(
        mode="wb", suffix=".colormap.txt", delete=False
    )
    taxids_tmp = tempfile.NamedTemporaryFile(
        mode="wb", suffix=".taxids.txt", delete=False
    )

    for tmp, key in [
        (tblout_tmp, "tblout"),
        (colormap_tmp, "colormap"),
        (taxids_tmp, "taxids"),
    ]:
        tmp.write(files[key].read())
        tmp.close()

    cmd = [
        sys.executable,
        "ete_profile.py",
        "-i",
        tblout_tmp.name,
        "-c",
        colormap_tmp.name,
        "-t",
        taxids_tmp.name,
        "-p",
        str(port),
    ]
    if max_val:
        cmd += ["-m", max_val]

    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"

    os.system(f"fuser -k {port}/tcp >/dev/null 2>&1")
    proc = subprocess.Popen(cmd, env=env)

    deadline = time.time() + 120
    connected = False
    while time.time() < deadline:
        try:
            with socket.create_connection(("localhost", port), timeout=2):
                connected = True
                break
        except OSError:
            time.sleep(3)

    if not connected:
        proc.kill()
        return jsonify({"error": "ETE profile server did not start in time."}), 500

    _ete4_ports.add(port)
    return jsonify({"ok": True, "port": port, "pid": proc.pid})


# ==============================================================================================================
#                                           ETE4 PROXY ROUTES
# Routes that forward the browser's requests to a locally-running ETE4 bottle server,
# so the browser never needs a direct SSH tunnel to the ETE4 port.
# ==============================================================================================================


@app.route("/ete4/gui")
def ete4_gui():
    """Serve ETE4's gui.html with asset paths and fetch() calls rewritten to go through Flask."""
    port = request.args.get("port", type=int)
    if not port:
        return "Missing ?port parameter", 400
    try:
        resp = _req.get(f"http://127.0.0.1:{port}/static/gui.html", timeout=5)
    except Exception as e:
        return f"ETE4 server on port {port} is not reachable: {e}", 502

    html = resp.text

    # Inject <base> so that relative paths in gui.html (e.g. href="gui.css",
    # href="images/icon.png") resolve through our proxy instead of to /ete4/*.
    html = html.replace("<head>", f'<head><base href="/ete4/{port}/static/">')

    # Rewrite root-relative /static/ asset references so they go through our proxy.
    # (<base> only handles relative URLs; root-relative ones still need rewrites.)
    html = html.replace('href="/static/', f'href="/ete4/{port}/static/')
    html = html.replace('src="/static/', f'src="/ete4/{port}/static/')

    # Inject a tiny fetch() interceptor before any other script runs.
    # ETE4's JS calls fetch("/trees/...") and fetch("/load") with root-relative paths.
    # We redirect those through /ete4/<port>/... so Flask can proxy them to the right server.
    intercept = f"""<script>
(function(){{
  var _f = window.fetch.bind(window);
  window.fetch = function(url, opts){{
    if (typeof url === 'string') {{
      if (url.startsWith('/trees') || url.startsWith('/load') || url.startsWith('/static/')) {{
        url = '/ete4/{port}' + url;
      }} else {{
        var _o = window.location.origin;
        if (url.startsWith(_o + '/static/') || url.startsWith(_o + '/trees') || url.startsWith(_o + '/load')) {{
          url = _o + '/ete4/{port}' + url.slice(_o.length);
        }}
      }}
    }}
    return _f(url, opts);
  }};
  // Remove 'port' from the visible URL so ETE4's set_query_string_values()
  // doesn't flag it as an unknown parameter and show the Swal warning.
  (function(){{
    var p = new URLSearchParams(location.search);
    p.delete('port');
    var qs = p.toString();
    history.replaceState(null, '', location.pathname + (qs ? '?' + qs : ''));
  }})();
  // Call reset_view() once the page is fully loaded and the iframe has
  // settled to its actual height. ETE4 initialises with a small offsetHeight
  // (whatever the iframe had before layout), so the initial zoom is wrong.
  // We re-fit the tree after the load event + a short pause so the async
  // init (init_trees, init_pixi) has completed and tree_size is populated.
  window.addEventListener('load', function() {{
    var attempts = 0;
    var t = setInterval(function() {{
      attempts++;
      if (typeof window.__ete4_reset_view === 'function') {{
        window.__ete4_reset_view();
        clearInterval(t);
      }} else if (attempts >= 30) {{
        clearInterval(t);
      }}
    }}, 200);
  }});
}})();
</script>
"""
    html = html.replace("</head>", intercept + "</head>")
    return html


@app.route("/ete4/<int:port>/static/<path:path>")
def ete4_static_proxy(port, path):
    """Forward ETE4 static file requests (JS, CSS, images) to the ETE4 bottle server."""
    try:
        resp = _req.get(f"http://127.0.0.1:{port}/static/{path}", timeout=10)
    except Exception:
        return "ETE4 server unavailable", 502
    # Expose reset_view globally so our injected script can call it after
    # the iframe reaches its final height.
    if path == "js/gui.js":
        patched = (
            resp.text
            + "\nwindow.__ete4_reset_view = function() { try { reset_view(); } catch(e) {} };\n"
        )
        return Response(patched, content_type="application/javascript")
    return Response(
        resp.content,
        content_type=resp.headers.get("Content-Type", "application/octet-stream"),
    )


@app.route("/static/images/<path:filename>")
def ete4_static_images_fallback(filename):
    """Catch /static/images/* requests from PixiJS web workers.

    PixiJS resolves the spritesheet PNG relative to the URL it stored before
    our fetch interceptor rewrote it, so the worker ends up requesting
    /static/images/spritesheet.png (no port in the path). We try every ETE4
    port launched this session and return the first successful response.
    """
    for port in list(_ete4_ports):
        try:
            resp = _req.get(
                f"http://127.0.0.1:{port}/static/images/{filename}", timeout=3
            )
            if resp.status_code == 200:
                return Response(
                    resp.content,
                    content_type=resp.headers.get(
                        "Content-Type", "application/octet-stream"
                    ),
                )
        except Exception:
            continue
    return "Not found", 404


@app.route("/ete4/<int:port>/trees", methods=["GET", "POST", "PUT", "DELETE"])
@app.route("/ete4/<int:port>/trees/<path:p>", methods=["GET", "POST", "PUT", "DELETE"])
def ete4_trees_proxy(port, p=""):
    """Forward all /trees/... API calls from ETE4's JS to the correct ETE4 server instance."""
    target = f"http://127.0.0.1:{port}/trees"
    if p:
        target += f"/{p}"
    # Strip Accept-Encoding so ETE4 does not brotli-compress the response;
    # this avoids having to decompress-then-re-compress in the proxy layer.
    fwd_headers = {
        k: v
        for k, v in request.headers
        if k.lower() not in ("host", "content-length", "accept-encoding")
    }
    try:
        resp = _req.request(
            method=request.method,
            url=target,
            params=request.args,
            headers=fwd_headers,
            data=request.get_data(),
            timeout=30,
            stream=True,
        )
        # Forward response headers except hop-by-hop ones Flask cannot re-send.
        skip = {"content-encoding", "content-length", "transfer-encoding", "connection"}
        out_headers = {k: v for k, v in resp.headers.items() if k.lower() not in skip}
        return Response(
            stream_with_context(resp.iter_content(chunk_size=4096)),
            status=resp.status_code,
            headers=out_headers,
        )
    except Exception as e:
        return jsonify({"error": str(e)}), 502


@app.route("/ete4/<int:port>/load", methods=["POST"])
def ete4_load_proxy(port):
    """Forward /load calls (used by ETE4's upload UI) to the ETE4 server."""
    try:
        resp = _req.post(
            f"http://127.0.0.1:{port}/load",
            json=request.get_json(),
            timeout=10,
        )
        return Response(
            resp.content,
            status=resp.status_code,
            content_type=resp.headers.get("Content-Type", "application/json"),
        )
    except Exception as e:
        return jsonify({"error": str(e)}), 502


def _pick_app_port(preferred, host="0.0.0.0", span=100):
    """Return `preferred` if it's free, otherwise the next free port above it.

    Lets several people run PhyloWave on the same machine without editing the
    file: if the port is taken, the app just moves to the next open one.
    (Named differently from `_find_free_port`, which picks ETE4 viewer ports.)
    """
    for port in range(preferred, preferred + span):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                s.bind((host, port))
                return port
            except OSError:
                continue  # in use — try the next one
    return preferred  # give up scanning; app.run will report the real error


if __name__ == "__main__":
    host = "0.0.0.0"
    preferred = int(os.environ.get("PORT", "8080"))
    port = _pick_app_port(preferred, host)
    if port != preferred:
        print(f"[PhyloWave] Port {preferred} is in use — using {port} instead.")
    print(f"[PhyloWave] Serving on http://localhost:{port}  (Ctrl-C to stop)")
    app.run(debug=False, host=host, port=port)
