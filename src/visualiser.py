"""
src/visualiser.py
-----------------
Generates interactive visualisations from the NetworkX graph.

  1. Pyvis  → self-contained interactive HTML (physics-driven force layout)
  2. Plotly → static / embeddable charts (degree distribution, heatmap, bar)
"""

from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import networkx as nx
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from pyvis.network import Network

from data.hpo_data import CATEGORY_COLORS, EFFECT_COLORS


# ── 1. Pyvis interactive network ──────────────────────────────────────────────

def render_pyvis(
    G: nx.Graph,
    output_path: str | Path,
    height: str = "750px",
    show_phenotypes: bool = False,
    physics: bool = True,
) -> Path:
    """
    Render the knowledge graph as a self-contained interactive HTML file.

    Args:
        G:               NetworkX graph from graph_builder.build_graph()
        output_path:     Where to write the .html file
        height:          Iframe height string (e.g. "750px")
        show_phenotypes: If False, hide HPO phenotype nodes to reduce clutter
        physics:         Enable force-directed layout simulation

    Returns:
        Path to generated HTML file
    """
    net = Network(
        height=height,
        width="100%",
        bgcolor="#0F172A",
        font_color="#E2E8F0",
        directed=False,
    )

    if physics:
        net.set_options("""
        {
          "physics": {
            "enabled": true,
            "barnesHut": {
              "gravitationalConstant": -8000,
              "centralGravity": 0.3,
              "springLength": 120,
              "springConstant": 0.04,
              "damping": 0.12,
              "avoidOverlap": 0.8
            },
            "maxVelocity": 50,
            "minVelocity": 0.1,
            "solver": "barnesHut",
            "stabilization": { "iterations": 200 }
          },
          "nodes": {
            "font": { "size": 12, "face": "Arial" },
            "borderWidth": 2,
            "shadow": { "enabled": true, "size": 6, "color": "rgba(0,0,0,0.5)" }
          },
          "edges": {
            "smooth": { "type": "continuous" },
            "shadow": false
          },
          "interaction": {
            "hover": true,
            "tooltipDelay": 150,
            "hideEdgesOnDrag": true
          }
        }
        """)

    for node_id, attrs in G.nodes(data=True):
        ntype = attrs.get("node_type", "gene")

        # optionally skip phenotype nodes
        if not show_phenotypes and ntype == "phenotype":
            continue

        label   = attrs.get("label", str(node_id))
        color   = attrs.get("color", "#64748B")
        size    = attrs.get("size", 14)
        shape   = attrs.get("shape", "dot")

        # build tooltip
        if ntype == "disease":
            diseases_str = attrs.get("description","")
            tip = (
                f"<b>{label}</b><br>"
                f"Type: Disease<br>"
                f"Category: {attrs.get('category','')}<br>"
                f"HPO: {attrs.get('hpo_id','')}<br>"
                f"OMIM: {attrs.get('omim','')}<br>"
                f"<i>{diseases_str[:80]}</i>"
            )
        elif ntype == "gene":
            dis_list = ", ".join(attrs.get("diseases", []))
            tip = (
                f"<b>{label}</b><br>"
                f"Type: Gene<br>"
                f"Diseases: {dis_list or '—'}<br>"
                f"Degree: {G.degree(node_id)}"
            )
        else:
            tip = f"<b>{label}</b><br>HPO phenotype term"

        net.add_node(
            node_id,
            label   = label if ntype != "phenotype" else "",
            title   = tip,
            color   = {"background": color, "border": "#1E293B",
                       "highlight": {"background": "#FBBF24", "border": "#F59E0B"}},
            size    = size,
            shape   = shape,
            font    = {"size": 11 if ntype == "phenotype" else 13,
                       "color": "#CBD5E1"},
        )

    for u, v, attrs in G.edges(data=True):
        ntype_u = G.nodes[u].get("node_type", "")
        ntype_v = G.nodes[v].get("node_type", "")

        # skip edges involving hidden phenotype nodes
        if not show_phenotypes and ("phenotype" in (ntype_u, ntype_v)):
            continue

        effect  = attrs.get("effect", "")
        pubmed  = attrs.get("pubmed")
        tip     = f"Effect: {effect}"
        if pubmed:
            tip += f"<br>PubMed: {pubmed}"

        net.add_edge(
            u, v,
            title = tip,
            color = {"color": attrs.get("color","#475569"),
                     "highlight": "#FBBF24"},
            width = attrs.get("width", 1),
        )

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    net.save_graph(str(out))

    # inject a legend + title banner into the HTML
    _inject_legend(out)
    return out


def _inject_legend(html_path: Path):
    """Inject a legend and title bar into the Pyvis HTML output."""
    legend_html = """
<style>
  #graph-header {
    position: fixed; top: 0; left: 0; right: 0; z-index: 9999;
    background: linear-gradient(135deg, #0F172A 0%, #1E3A5F 100%);
    padding: 10px 20px; display: flex; align-items: center; gap: 20px;
    border-bottom: 1px solid #334155; font-family: Arial, sans-serif;
  }
  #graph-header h1 { color: #F1F5F9; font-size: 16px; margin: 0; flex: 1; }
  .legend-group { display: flex; align-items: center; gap: 14px; flex-wrap: wrap; }
  .legend-item { display: flex; align-items: center; gap: 5px;
                 font-size: 11px; color: #CBD5E1; }
  .legend-dot  { width: 12px; height: 12px; border-radius: 50%; flex-shrink: 0; }
  .legend-diamond { width: 12px; height: 12px; transform: rotate(45deg); flex-shrink: 0; }
  #mynetwork   { margin-top: 52px !important; }
  #graph-footer {
    position: fixed; bottom: 0; left: 0; right: 0; z-index: 9999;
    background: #0F172A; padding: 4px 16px;
    font-size: 11px; color: #64748B; font-family: Arial, sans-serif;
    border-top: 1px solid #1E293B;
  }
</style>
<div id="graph-header">
  <h1>🧬 HPO Gene–Phenotype Knowledge Graph</h1>
  <div class="legend-group">
    <span style="font-size:11px;color:#94A3B8;font-weight:600">NODE TYPE</span>
    <div class="legend-item">
      <div class="legend-diamond" style="background:#3B82F6"></div> Disease
    </div>
    <div class="legend-item">
      <div class="legend-dot" style="background:#64748B"></div> Gene (single disease)
    </div>
    <div class="legend-item">
      <div class="legend-dot" style="background:#F97316"></div> Gene (pleiotropic)
    </div>
  </div>
  <div class="legend-group">
    <span style="font-size:11px;color:#94A3B8;font-weight:600">EDGE EFFECT</span>
    <div class="legend-item"><div class="legend-dot" style="background:#DC2626;border-radius:0;height:3px;width:20px"></div> Causal</div>
    <div class="legend-item"><div class="legend-dot" style="background:#F97316;border-radius:0;height:3px;width:20px"></div> Risk</div>
    <div class="legend-item"><div class="legend-dot" style="background:#22C55E;border-radius:0;height:3px;width:20px"></div> Protective</div>
  </div>
</div>
<div id="graph-footer">
  Data: HPO · OMIM · ClinVar · GWAS Catalog &nbsp;|&nbsp;
  Drag nodes to rearrange · Scroll to zoom · Hover for details &nbsp;|&nbsp;
  Built by Rayala Madhu Bhanu Varma · GenepoweRx
</div>
"""
    content = html_path.read_text(encoding="utf-8")
    # inject after <body>
    content = content.replace("<body>", "<body>\n" + legend_html, 1)
    html_path.write_text(content, encoding="utf-8")


# ── 2. Plotly charts ─────────────────────────────────────────────────────────

def build_degree_distribution(G: nx.Graph) -> go.Figure:
    """Bar chart of gene node degrees (connectivity)."""
    gene_degrees = sorted(
        [(n, G.degree(n), G.nodes[n].get("diseases", []))
         for n, d in G.nodes(data=True) if d.get("node_type") == "gene"],
        key=lambda x: -x[1]
    )[:25]  # top 25

    genes      = [x[0] for x in gene_degrees]
    degrees    = [x[1] for x in gene_degrees]
    is_pleo    = [len(x[2]) > 1 for x in gene_degrees]
    colors     = ["#F97316" if p else "#3B82F6" for p in is_pleo]
    hover      = [f"<b>{g}</b><br>Degree: {d}<br>Diseases: {', '.join(dis)}"
                  for g, d, dis in gene_degrees
                  for dis in [G.nodes[g].get("diseases", [])]]

    fig = go.Figure(go.Bar(
        x=genes, y=degrees, marker_color=colors,
        customdata=[[", ".join(G.nodes[g].get("diseases", []))] for g in genes],
        hovertemplate="<b>%{x}</b><br>Graph degree: %{y}<br>Diseases: %{customdata[0]}<extra></extra>",
        text=degrees, textposition="outside",
    ))
    fig.update_layout(
        title=dict(text="Top 25 Genes by Graph Degree<br><sup>Orange = pleiotropic (linked to >1 disease)</sup>",
                   font=dict(size=15), x=0.5),
        xaxis=dict(title="", tickangle=-40, tickfont=dict(size=10)),
        yaxis=dict(title="Graph Degree (connections)", gridcolor="#f1f5f9"),
        plot_bgcolor="white", paper_bgcolor="white",
        height=420, margin=dict(l=60, r=20, t=80, b=100),
    )
    return fig


def build_disease_gene_heatmap(G: nx.Graph) -> go.Figure:
    """Heatmap: disease categories × gene count."""
    from data.hpo_data import HPO_DATA, CATEGORY_COLORS
    from collections import defaultdict

    cat_genes: dict[str, set] = defaultdict(set)
    cat_diseases: dict[str, list] = defaultdict(list)

    for disease, info in HPO_DATA.items():
        cat = info.get("category", "Other")
        cat_diseases[cat].append(disease)
        for g in info["genes"]:
            cat_genes[cat].add(g)

    categories = sorted(cat_genes.keys())
    data_rows = []
    for cat in categories:
        data_rows.append({
            "category":     cat,
            "n_genes":      len(cat_genes[cat]),
            "n_diseases":   len(cat_diseases[cat]),
            "color":        CATEGORY_COLORS.get(cat, "#94A3B8"),
        })
    df = pd.DataFrame(data_rows).sort_values("n_genes", ascending=True)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=df["category"], x=df["n_genes"],
        orientation="h",
        marker_color=df["color"].tolist(),
        text=df["n_genes"],
        textposition="outside",
        customdata=df[["n_diseases"]].values,
        hovertemplate="<b>%{y}</b><br>Genes: %{x}<br>Diseases: %{customdata[0]}<extra></extra>",
        name="Unique genes",
    ))
    fig.update_layout(
        title=dict(text="Unique Genes per Disease Category", font=dict(size=15), x=0.5),
        xaxis=dict(title="Number of unique genes", gridcolor="#f1f5f9"),
        yaxis=dict(title=""),
        plot_bgcolor="white", paper_bgcolor="white",
        height=400, margin=dict(l=140, r=60, t=70, b=50),
    )
    return fig


def build_pleiotropy_chart(G: nx.Graph) -> go.Figure:
    """Horizontal bar of pleiotropic genes and which diseases they bridge."""
    rows = []
    for n, d in G.nodes(data=True):
        if d.get("node_type") != "gene":
            continue
        diseases = list(set(d.get("diseases", [])))
        if len(diseases) > 1:
            rows.append({"gene": n, "n_diseases": len(diseases),
                         "diseases": ", ".join(sorted(diseases))})
    if not rows:
        fig = go.Figure()
        fig.add_annotation(text="No pleiotropic genes in current selection",
                           showarrow=False, font=dict(size=14))
        return fig

    df = pd.DataFrame(rows).sort_values("n_diseases", ascending=True)
    fig = go.Figure(go.Bar(
        y=df["gene"], x=df["n_diseases"],
        orientation="h",
        marker_color="#F97316",
        text=df["n_diseases"],
        textposition="outside",
        customdata=df[["diseases"]].values,
        hovertemplate="<b>%{y}</b><br>Bridges %{x} diseases:<br>%{customdata[0]}<extra></extra>",
    ))
    fig.update_layout(
        title=dict(text="Pleiotropic Genes — Bridging Multiple Diseases",
                   font=dict(size=15), x=0.5),
        xaxis=dict(title="Number of diseases", gridcolor="#f1f5f9"),
        yaxis=dict(title=""),
        plot_bgcolor="white", paper_bgcolor="white",
        height=max(300, len(df) * 30 + 100),
        margin=dict(l=90, r=60, t=70, b=50),
    )
    return fig


def build_category_pie(G: nx.Graph) -> go.Figure:
    """Pie chart of diseases by category."""
    from data.hpo_data import HPO_DATA, CATEGORY_COLORS
    from collections import Counter
    cats = Counter(v.get("category","Other") for v in HPO_DATA.values())
    labels = list(cats.keys())
    values = list(cats.values())
    colors = [CATEGORY_COLORS.get(l,"#94A3B8") for l in labels]
    fig = go.Figure(go.Pie(
        labels=labels, values=values,
        marker=dict(colors=colors, line=dict(color="#fff", width=2)),
        textinfo="label+value",
        hovertemplate="<b>%{label}</b><br>%{value} diseases (%{percent})<extra></extra>",
        hole=0.4,
    ))
    fig.update_layout(
        title=dict(text="Disease Distribution by Category<br><sup>32 conditions across 10 categories</sup>",
                   font=dict(size=15), x=0.5),
        paper_bgcolor="white",
        height=420,
        margin=dict(l=20, r=20, t=80, b=20),
        legend=dict(font=dict(size=11)),
    )
    return fig
