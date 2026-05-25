import json
import os
import re
import hashlib


def parse_itol_colors(itol_colors_path: str) -> dict:
    colors = {}
    in_data = False
    try:
        with open(itol_colors_path) as f:
            for line in f:
                line = line.strip()
                if line == "DATA":
                    in_data = True
                    continue
                if in_data and line and not line.startswith("#"):
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        colors[parts[0]] = parts[1]
    except Exception:
        pass
    return colors


def auto_leaf_colors(newick_str: str) -> dict:
    leaf_names = re.findall(r"[\(,]([^\(\),;:]+)(?::[0-9.e+-]+)?", newick_str)
    leaf_names = [
        n.strip() for n in leaf_names if n.strip() and not n.strip().startswith("(")
    ]
    palette = [
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
        "#fabed4",
        "#fffac8",
    ]
    taxon_color = {}
    for name in leaf_names:
        taxon = name.split(".")[0]
        if taxon not in taxon_color:
            idx = int(hashlib.md5(taxon.encode()).hexdigest()[:4], 16) % len(palette)
            taxon_color[taxon] = palette[idx]
    return {name: taxon_color.get(name.split(".")[0], "#cccccc") for name in leaf_names}


def build_tree_html(
    newick_str: str,
    leaf_colors: dict = None,
    heatmap_data: dict = None,
    heatmap_profiles: list = None,
    title: str = "Phylogenetic Tree",
    height: int = 700,
) -> str:

    if leaf_colors is None:
        leaf_colors = auto_leaf_colors(newick_str)

    newick_json = json.dumps(newick_str)
    colors_json = json.dumps(leaf_colors)
    heatmap_json = json.dumps(heatmap_data or {})
    profiles_json = json.dumps(heatmap_profiles or [])
    title_json = json.dumps(title)

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ background: #f9f9f9; font-family: 'Segoe UI', Arial, sans-serif; user-select: none; overflow: auto; }}
  #hdr {{ display: flex; align-items: center; gap: 8px; padding: 5px 10px; background: white; border-bottom: 1px solid #e0e0e0; height: 34px; flex-shrink: 0; }}
  #hdr h3 {{ font-size: 12px; color: #333; margin-right: 4px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; max-width: 300px; }}
  .btn {{ padding: 2px 8px; font-size: 10px; border: 1px solid #ccc; border-radius: 3px; background: white; cursor: pointer; white-space: nowrap; transition: background 0.2s; }}
  .btn:hover {{ background: #f0f0f0; }}
  .btn.active {{ background: #d0e8ff; border-color: #80bdff; color: #004085; font-weight: bold; }}
  #leaf-count {{ font-size: 10px; color: #999; margin-left: auto; white-space: nowrap; }}
  #wrap {{ width: 100%; min-height: 500px; }}
  svg {{ display: block; }}
  .link {{ fill: none; stroke: #888; stroke-width: 0.8px; }}
  .link-dash {{ fill: none; stroke: #ccc; stroke-width: 0.5px; stroke-dasharray: 2,2; }}
  .node-circle {{ stroke: white; stroke-width: 1.2px; cursor: pointer; }}
  .leaf-label {{ font-size: 9px; fill: #333; dominant-baseline: middle; cursor: default; }}
  .bootstrap {{ font-size: 7px; fill: #c0392b; dominant-baseline: auto; pointer-events: none; }}
  .heat-cell {{ stroke: white; stroke-width: 0.4px; }}
  .col-header {{ font-size: 8px; fill: #555; }}
  #tip {{ position: fixed; background: rgba(20,20,20,0.88); color: #eee; padding: 5px 9px; border-radius: 4px; font-size: 10px; pointer-events: none; display: none; z-index: 999; max-width: 220px; line-height: 1.6; white-space: pre-line; }}
  #error-box {{ color: #d32f2f; background: #ffebee; padding: 15px; margin: 15px; border-radius: 5px; font-family: monospace; display: none; border: 1px solid #ffcdd2; }}
</style>
</head>
<body>

<div id="hdr">
  <h3 id="title-text">Loading...</h3>
  <button class="btn" id="btn-reset">⌖ Reset</button>
  <button class="btn" id="btn-ladder">⇅ Ladderize</button>
  <button class="btn" id="btn-align">≈ Align Leaves</button>
  <button class="btn" id="btn-exp">↓ SVG</button>
  <span id="leaf-count"></span>
</div>
<div id="error-box"></div>
<div id="wrap"><svg id="svg"></svg></div>
<div id="tip"></div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"></script>
<script>
"use strict";

try {{
    const NEWICK   = {newick_json};
    const COLORS   = {colors_json};
    const HEATMAP  = {heatmap_json};
    const PROFILES = {profiles_json};
    const TITLE    = {title_json};

    document.getElementById("title-text").textContent = TITLE || "Phylogenetic Tree";

    const LABEL_W   = 160;
    const STRIP_W   = 10;
    const STRIP_GAP = 4;
    const CELL_W    = 13;
    const CELL_H    = 12;
    const COL_HDR_H = 60;
    const LEG_H     = 30;
    const MARGIN    = {{ top: COL_HDR_H + 10, right: 20, bottom: LEG_H + 10, left: 20 }};

    const HEAT_ABSENT  = "#e8e8e8";
    const HEAT_PRESENT = "#27ae60";

    function parseNewick(s) {{
      if (!s) throw new Error("The Newick string is empty.");
      s = s.trim().replace(/;$/, "");
      const anc = [];
      let node = {{}};
      const toks = s.split(/\\s*([(),:])\\s*/);
      
      for (let i = 0; i < toks.length; i++) {{
        const tok = toks[i];
        if (!tok && tok !== "0") continue;
        const prev = toks[i - 1] || "";
        
        switch (tok) {{
          case "(":
            const c0 = {{}};
            if(!node.children) node.children = [];
            node.children.push(c0);
            anc.push(node);
            node = c0;
            break;
          case ",":
            const sib = {{}};
            if(anc.length === 0) throw new Error("Malformed Newick");
            anc[anc.length - 1].children.push(sib);
            node = sib;
            break;
          case ")":
            if(anc.length === 0) throw new Error("Malformed Newick");
            node = anc.pop();
            break;
          case ":":
            break;
          default:
            if (prev === ")" || prev === "(" || prev === ",") {{
              node.name = tok;
            }} else if (prev === ":") {{
              node.branchLength = parseFloat(tok) || 0;
            }}
        }}
      }}
      return node;
    }}

    function ladderize(node) {{
      const kids = node.children || node._children;
      if (!kids) return 1;
      const sizes = kids.map(ladderize);
      const order = sizes.map((s, i) => [s, i]).sort((a, b) => b[0] - a[0]).map(x => x[1]);
      const sorted = order.map(i => kids[i]);
      if (node.children) node.children = sorted;
      if (node._children) node._children = sorted;
      return sizes.reduce((a, b) => a + b, 0);
    }}

    // Modified to read true underlying tree data so scale stays consistent when collapsed
    function maxTreeDepth(n) {{
      const bl = n.branchLength || 0;
      const kids = n.children || n._children;
      if (!kids || kids.length === 0) return bl;
      return bl + Math.max(...kids.map(maxTreeDepth));
    }}

    let treeData = parseNewick(NEWICK);
    let isLadderized = false;
    let isAligned = false; 

    const svg   = d3.select("#svg");
    const wrap  = document.getElementById("wrap");
    const tip   = document.getElementById("tip");

    let zoomState = d3.zoomIdentity;

    function render() {{
      svg.selectAll("*").remove();
      // Rebuild hierarchy. If a node has _children instead of children, D3 treats it as a leaf!
      let root = d3.hierarchy(treeData, d => d.children);

      const nLeaves = root.leaves().length;
      if (nLeaves === 0) throw new Error("Tree parsed successfully but has 0 leaves.");
      
      const ROW_H   = Math.max(CELL_H, Math.min(24, Math.floor(500 / nLeaves)));
      const treeH   = nLeaves * ROW_H;
      const HEAT_W  = PROFILES.length * CELL_W;
      const svgW    = Math.max(600, wrap.clientWidth || 800);
      const treeW   = svgW - MARGIN.left - MARGIN.right - LABEL_W - STRIP_GAP - STRIP_W - CELL_W - HEAT_W;
      const svgH    = treeH + MARGIN.top + MARGIN.bottom;

      svg.attr("width", svgW).attr("height", svgH);
      wrap.style.height = svgH + "px";
      document.getElementById("leaf-count").textContent =
        nLeaves + " visible leaves" + (PROFILES.length ? " · " + PROFILES.length + " profiles" : "");

      d3.cluster().size([treeH, treeW])(root);

      if (!isAligned) {{
          // Calculate max depth from ORIGINAL data so the tree doesn't jump when collapsed
          const totalDepth = maxTreeDepth(treeData);
          const scale      = totalDepth > 0 ? treeW / totalDepth : 1;

          function setY(node, parentY) {{
            node.y = parentY + (node.data.branchLength || 0) * scale;
            if (node.children) node.children.forEach(c => setY(c, node.y));
          }}
          setY(root, 0);
      }}

      const g = svg.append("g").attr("class", "root-g");

      const zoom = d3.zoom()
        .scaleExtent([0.2, 20])
        .on("zoom", ev => {{
          zoomState = ev.transform;
          g.attr("transform", `translate(${{ev.transform.x + MARGIN.left}},${{ev.transform.y + MARGIN.top}}) scale(${{ev.transform.k}})`);
        }});

      svg.call(zoom).call(zoom.transform, zoomState.translate(0, 0));
      g.attr("transform", `translate(${{zoomState.x + MARGIN.left}},${{zoomState.y + MARGIN.top}}) scale(${{zoomState.k}})`);

      const linkG = g.append("g");
      root.links().forEach(link => {{
        const s = link.source, t = link.target;
        linkG.append("path").attr("class", "link").attr("d", `M${{s.y}},${{s.x}} H${{t.y}} V${{t.x}}`);
      }});

      if (!isAligned) {{
          root.leaves().forEach(leaf => {{
            if (leaf.y < treeW - 1) {{
              // Do not draw dashed lines for collapsed parent nodes
              if (!leaf.data._children) {{
                  g.append("line").attr("class", "link-dash").attr("x1", leaf.y).attr("y1", leaf.x).attr("x2", treeW).attr("y2", leaf.x);
              }}
            }}
          }});
      }}

      const nodeG = g.append("g");
      const nodeEnter = nodeG.selectAll("g.nd")
        .data(root.descendants())
        .enter().append("g")
        .attr("class", "nd")
        .attr("transform", d => `translate(${{d.y}},${{d.x}})`);

      nodeEnter.append("circle")
        .attr("class", "node-circle")
        .attr("r", d => d.data.children || d.data._children ? (d.data._children ? 4 : 2.5) : 1.8)
        .attr("fill", d => d.data.children || d.data._children ? (d.data._children ? "#e74c3c" : "#4a90d9") : "#999")
        .on("mouseover", (ev, d) => showTip(ev, d))
        .on("mouseout", () => hideTip())
        .on("click", (ev, d) => {{
          // Mutate the original raw data, not the D3 hierarchy!
          const dataNode = d.data;
          if (!dataNode.children && !dataNode._children) return; 
          ev.stopPropagation();
          
          if (dataNode._children) {{ 
            dataNode.children = dataNode._children; 
            dataNode._children = null; 
          }} else {{ 
            dataNode._children = dataNode.children; 
            dataNode.children = null; 
          }}
          render();
        }});

      nodeEnter.filter(d => d.children && d.data.name && !isNaN(+d.data.name))
        .append("text").attr("class", "bootstrap").attr("dx", -3).attr("dy", -3).attr("text-anchor", "end").text(d => (+d.data.name).toFixed(0));

      const leaves = root.leaves();
      leaves.forEach(leaf => {{
        // Don't draw labels/heatmap for collapsed groups
        if (leaf.data._children) return; 

        const lx = treeW, ly = leaf.x, name = leaf.data.name || "";
        const displayName = name.includes(".") ? name.split(".").slice(1).join(".") : name;
        
        g.append("text").attr("class", "leaf-label").attr("x", lx + 4).attr("y", ly).text(displayName)
         .on("mouseover", ev => showTipText(ev, formatTooltip(name, leaf))).on("mouseout", () => hideTip());

        const stripX = lx + LABEL_W;
        const color  = COLORS[name] || "#cccccc";
        g.append("rect").attr("x", stripX).attr("y", ly - STRIP_W / 2).attr("width", STRIP_W).attr("height", STRIP_W).attr("fill", color).attr("rx", 1)
         .on("mouseover", ev => showTipText(ev, "Taxon: " + name.split(".")[0])).on("mouseout", () => hideTip());

        PROFILES.forEach((prof, pi) => {{
          const cellX = stripX + STRIP_W + CELL_W / 2 + pi * CELL_W;
          const val   = (HEATMAP[name] && HEATMAP[name][prof]) ? 1 : 0;
          g.append("rect").attr("class", "heat-cell").attr("x", cellX - CELL_W / 2 + 0.5).attr("y", ly - CELL_H / 2)
           .attr("width", CELL_W - 1).attr("height", CELL_H - 1).attr("fill", val ? HEAT_PRESENT : HEAT_ABSENT)
           .on("mouseover", ev => showTipText(ev, prof + "\\n" + name + "\\n" + (val ? "✓ present" : "✗ absent"))).on("mouseout", () => hideTip());
        }});
      }});

      if (PROFILES.length > 0) {{
        const stripX = treeW + LABEL_W;
        PROFILES.forEach((prof, pi) => {{
          const cx = stripX + STRIP_W + CELL_W / 2 + pi * CELL_W;
          g.append("text").attr("class", "col-header").attr("transform", `translate(${{cx}}, -4) rotate(-55)`).attr("text-anchor", "start").text(prof);
        }});
      }}

      const taxonColors = {{}};
      root.leaves().forEach(leaf => {{
        if (leaf.data._children) return;
        const name = leaf.data.name || "", taxon = name.split(".")[0];
        if (!taxonColors[taxon]) taxonColors[taxon] = COLORS[name] || "#cccccc";
      }});

      const legendY = treeH + 12;
      let legendX = 0;
      Object.keys(taxonColors).forEach(taxon => {{
        const entry = g.append("g").attr("transform", `translate(${{legendX}}, ${{legendY}})`);
        entry.append("rect").attr("width", 8).attr("height", 8).attr("rx", 1).attr("fill", taxonColors[taxon]);
        entry.append("text").attr("x", 11).attr("y", 8).attr("font-size", "8px").attr("fill", "#444").text(taxon);
        legendX += Math.max(40, taxon.length * 6 + 18);
      }});
    }}

    function formatTooltip(name, leaf) {{
      const parts = name.split(".");
      const taxon = parts[0], acc = parts.slice(1).join(".");
      let txt = (acc || name);
      if (taxon && acc) txt += "\\nTaxon: " + taxon;
      if (leaf.data.branchLength != null) txt += "\\nBranch length: " + leaf.data.branchLength.toFixed(6);
      return txt;
    }}

    function showTipText(ev, text) {{ tip.style.display = "block"; tip.textContent = text; moveTip(ev); }}
    
    function showTip(ev, d) {{
      const name = d.data.name || "";
      let txt = "";
      
      // Use raw data to count total leaves hidden inside a collapsed node
      function countLeaves(n) {{
         if (!n.children && !n._children) return 1;
         let sum = 0;
         const arr = n.children || n._children || [];
         for(let i=0; i<arr.length; i++) sum += countLeaves(arr[i]);
         return sum;
      }}

      if (d.data.children || d.data._children) {{
        txt = "Internal node\\n";
        if (name && !isNaN(+name)) txt += "Bootstrap: " + (+name).toFixed(0) + "\\n";
        txt += "Subtree leaves: " + countLeaves(d.data);
        txt += d.data._children ? "\\n(collapsed — click to expand)" : "\\n(click to collapse)";
      }} else {{ txt = formatTooltip(name, d); }}
      
      if (!txt) return;
      tip.style.display = "block"; tip.textContent = txt; moveTip(ev);
    }}

    function hideTip() {{ tip.style.display = "none"; }}
    function moveTip(ev) {{ tip.style.left = (ev.clientX + 12) + "px"; tip.style.top = (ev.clientY - 8) + "px"; }}
    svg.on("mousemove", moveTip);

    document.getElementById("btn-reset").onclick = () => {{
      zoomState = d3.zoomIdentity;
      svg.transition().duration(350).call(d3.zoom().transform, d3.zoomIdentity);
      render();
    }};

    document.getElementById("btn-ladder").onclick = () => {{ isLadderized = !isLadderized; if (isLadderized) ladderize(treeData); render(); }};
    
    document.getElementById("btn-align").onclick = () => {{ 
      isAligned = !isAligned; 
      const btn = document.getElementById("btn-align");
      if(isAligned) btn.classList.add("active");
      else btn.classList.remove("active");
      render(); 
    }};

    document.getElementById("btn-exp").onclick = () => {{
      const a = document.createElement("a");
      a.href = URL.createObjectURL(new Blob([svg.node().outerHTML], {{type: "image/svg+xml"}}));
      a.download = "tree.svg";
      a.click();
    }};

    render();
    window.addEventListener("resize", () => {{ zoomState = d3.zoomIdentity; render(); }});

}} catch (error) {{
    console.error(error);
    const errBox = document.getElementById("error-box");
    errBox.style.display = "block";
    errBox.innerHTML = "<strong>ERROR:</strong><br><br>" + error.message;
}}
</script>
</body>
</html>"""
