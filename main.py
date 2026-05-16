#!/usr/bin/env python3
"""
main.py  —  HPO Gene–Phenotype Knowledge Graph Builder
=======================================================
Builds the NetworkX knowledge graph from 32 HPO-annotated conditions
and generates all interactive HTML outputs.

Usage
-----
    # Full graph (all 32 conditions, genes only — no phenotype nodes)
    python main.py

    # With HPO phenotype nodes visible
    python main.py --phenotypes

    # Filter to specific conditions
    python main.py --conditions "Type 2 Diabetes" "Alzheimer Disease"

    # Analyse a specific gene
    python main.py --gene PTPN22

    # Find path between two genes
    python main.py --path PTPN22 STAT4
"""

import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from src.graph_builder import (
    build_graph, get_graph_stats, get_gene_summary,
    get_disease_summary, shortest_path_between_genes, shared_diseases,
)
from src.visualiser import (
    render_pyvis,
    build_degree_distribution,
    build_disease_gene_heatmap,
    build_pleiotropy_chart,
    build_category_pie,
)
from data.hpo_data import HPO_DATA


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hpo_explorer",
        description="HPO Gene–Disease Knowledge Graph · 32 conditions · 96 genes",
    )
    p.add_argument("--phenotypes", action="store_true",
                   help="Include HPO phenotype term nodes in the graph")
    p.add_argument("--conditions", nargs="+", metavar="CONDITION",
                   help="Subset of conditions to include")
    p.add_argument("--no-physics", action="store_true",
                   help="Disable force-directed layout (faster render)")
    p.add_argument("--gene", metavar="GENE",
                   help="Print summary for a specific gene")
    p.add_argument("--path", nargs=2, metavar=("GENE_A","GENE_B"),
                   help="Find shortest path between two genes in the graph")
    p.add_argument("--output-dir", default="output",
                   help="Output directory (default: output/)")
    p.add_argument("--list-conditions", action="store_true",
                   help="List all 32 available conditions and exit")
    return p


def main() -> int:
    args   = build_parser().parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── List conditions ────────────────────────────────────────────────────
    if args.list_conditions:
        print("\n32 HPO-annotated conditions in this dataset:\n")
        from collections import defaultdict
        by_cat: dict = defaultdict(list)
        for name, info in HPO_DATA.items():
            by_cat[info.get("category","Other")].append(name)
        for cat in sorted(by_cat):
            print(f"  [{cat}]")
            for d in sorted(by_cat[cat]):
                ng = len(HPO_DATA[d]["genes"])
                print(f"    • {d}  ({ng} gene{'s' if ng!=1 else ''})")
        print()
        return 0

    # ── Build graph ────────────────────────────────────────────────────────
    print(f"\n[1/5] Building knowledge graph …")
    G = build_graph(
        include_phenotype_nodes = args.phenotypes,
        conditions              = args.conditions,
    )
    stats = get_graph_stats(G)
    print(f"      Diseases  : {stats['diseases']}")
    print(f"      Genes     : {stats['genes']}")
    print(f"      Phenotypes: {stats['phenotypes']}")
    print(f"      Edges     : {stats['total_edges']}")
    print(f"      Pleiotropic genes: {stats['pleiotropic_genes']}")
    print(f"      → {stats['pleiotropic_list']}")

    # ── Gene query ─────────────────────────────────────────────────────────
    if args.gene:
        g = args.gene.upper()
        if g not in G.nodes:
            print(f"\nGene '{g}' not found in graph.")
            return 1
        d = G.nodes[g]
        print(f"\n── Gene: {g} ──────────────────────────────")
        print(f"  Diseases   : {', '.join(d.get('diseases',[]))}")
        print(f"  Degree     : {G.degree(g)}")
        print(f"  Pleiotropic: {len(d.get('diseases',[])) > 1}")
        return 0

    # ── Path query ─────────────────────────────────────────────────────────
    if args.path:
        a, b = args.path[0].upper(), args.path[1].upper()
        path = shortest_path_between_genes(G, a, b)
        shared = shared_diseases(G, a, b)
        print(f"\n── Path: {a} → {b} ─────────────────────────")
        if path:
            print(f"  Shortest path ({len(path)-1} hops): {' → '.join(path)}")
        else:
            print("  No path found.")
        print(f"  Shared diseases: {shared or 'none'}")
        return 0

    # ── Pyvis interactive graph ────────────────────────────────────────────
    print(f"\n[2/5] Rendering interactive network graph …")
    graph_path = outdir / "hpo_knowledge_graph.html"
    render_pyvis(
        G,
        output_path     = graph_path,
        show_phenotypes = args.phenotypes,
        physics         = not args.no_physics,
    )
    print(f"      → {graph_path}")

    # ── Plotly charts → HTML ───────────────────────────────────────────────
    print(f"[3/5] Building analytics charts …")
    charts = {
        "degree_distribution.html": build_degree_distribution(G),
        "genes_per_category.html":  build_disease_gene_heatmap(G),
        "pleiotropy.html":          build_pleiotropy_chart(G),
        "category_pie.html":        build_category_pie(G),
    }
    for fname, fig in charts.items():
        p = outdir / fname
        fig.write_html(str(p), include_plotlyjs="cdn", full_html=True)
        print(f"      → {p}")

    # ── Summary tables ─────────────────────────────────────────────────────
    print(f"[4/5] Writing summary tables …")
    import csv
    gene_rows = get_gene_summary(G)
    with open(outdir / "gene_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=gene_rows[0].keys())
        w.writeheader(); w.writerows(gene_rows)
    disease_rows = get_disease_summary(G)
    with open(outdir / "disease_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=disease_rows[0].keys())
        w.writeheader(); w.writerows(disease_rows)
    print(f"      → {outdir}/gene_summary.csv  ({len(gene_rows)} genes)")
    print(f"      → {outdir}/disease_summary.csv  ({len(disease_rows)} diseases)")

    # ── Dashboard index page ───────────────────────────────────────────────
    print(f"[5/5] Building dashboard index …")
    _build_index(outdir, stats)
    print(f"      → {outdir}/index.html")

    print(f"\n✅ All outputs written to: {outdir.resolve()}/")
    print(f"   Open {outdir}/index.html in your browser to explore.\n")
    return 0


def _build_index(outdir: Path, stats: dict):
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>HPO Knowledge Graph — Dashboard</title>
<style>
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'Segoe UI',system-ui,sans-serif;background:#0F172A;color:#E2E8F0;line-height:1.6}}
  .header{{background:linear-gradient(135deg,#1E3A5F,#1a5276);padding:32px 40px;border-bottom:1px solid #334155}}
  .header h1{{font-size:24px;font-weight:700;color:#F1F5F9}}
  .header p{{color:#94A3B8;font-size:14px;margin-top:6px}}
  .container{{max-width:1100px;margin:0 auto;padding:32px 24px}}
  .stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:16px;margin-bottom:36px}}
  .stat{{background:#1E293B;border-radius:12px;padding:18px;text-align:center;border:1px solid #334155}}
  .stat .num{{font-size:30px;font-weight:700;color:#60A5FA}}
  .stat .lbl{{font-size:11px;color:#94A3B8;text-transform:uppercase;letter-spacing:.5px;margin-top:4px}}
  .grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:20px}}
  .card{{background:#1E293B;border-radius:14px;padding:24px;border:1px solid #334155;
         text-decoration:none;color:inherit;display:block;transition:border-color .2s,transform .15s}}
  .card:hover{{border-color:#60A5FA;transform:translateY(-2px)}}
  .card-icon{{font-size:32px;margin-bottom:12px}}
  .card h3{{font-size:16px;font-weight:600;color:#F1F5F9;margin-bottom:6px}}
  .card p{{font-size:13px;color:#94A3B8}}
  .footer{{text-align:center;color:#475569;font-size:12px;margin-top:48px;padding-bottom:24px}}
  .badge{{display:inline-block;background:#1E3A5F;color:#60A5FA;border-radius:99px;
          padding:2px 10px;font-size:11px;margin:2px}}
</style>
</head>
<body>
<div class="header">
  <h1>🧬 HPO Gene–Phenotype Knowledge Graph</h1>
  <p>Interactive exploration of gene–disease associations across 32 medical conditions · HPO · OMIM · ClinVar · GWAS Catalog</p>
  <div style="margin-top:12px">
    <span class="badge">32 Conditions</span>
    <span class="badge">96 Unique Genes</span>
    <span class="badge">108 Gene–Disease Links</span>
    <span class="badge">11 Pleiotropic Genes</span>
    <span class="badge">10 Disease Categories</span>
  </div>
</div>
<div class="container">
  <div class="stats">
    <div class="stat"><div class="num">{stats['diseases']}</div><div class="lbl">Diseases</div></div>
    <div class="stat"><div class="num">{stats['genes']}</div><div class="lbl">Genes</div></div>
    <div class="stat"><div class="num">{stats['gene_disease_links']}</div><div class="lbl">Gene–Disease Links</div></div>
    <div class="stat"><div class="num">{stats['pleiotropic_genes']}</div><div class="lbl">Pleiotropic Genes</div></div>
    <div class="stat"><div class="num">{stats['total_edges']}</div><div class="lbl">Graph Edges</div></div>
    <div class="stat"><div class="num">{stats['density']}</div><div class="lbl">Graph Density</div></div>
  </div>
  <div class="grid">
    <a class="card" href="hpo_knowledge_graph.html" target="_blank">
      <div class="card-icon">🕸️</div>
      <h3>Interactive Knowledge Graph</h3>
      <p>Force-directed Pyvis network · Drag, zoom, hover for details · Gene–disease–phenotype connections</p>
    </a>
    <a class="card" href="degree_distribution.html" target="_blank">
      <div class="card-icon">📊</div>
      <h3>Gene Connectivity Chart</h3>
      <p>Top 25 genes by graph degree · Orange = pleiotropic (bridges multiple diseases)</p>
    </a>
    <a class="card" href="genes_per_category.html" target="_blank">
      <div class="card-icon">🗂️</div>
      <h3>Genes per Disease Category</h3>
      <p>How many unique genes are linked to each disease category (Metabolic, Neurological, etc.)</p>
    </a>
    <a class="card" href="pleiotropy.html" target="_blank">
      <div class="card-icon">🔀</div>
      <h3>Pleiotropic Gene Analysis</h3>
      <p>Genes that connect multiple diseases — key drug targets and shared pathways</p>
    </a>
    <a class="card" href="category_pie.html" target="_blank">
      <div class="card-icon">🥧</div>
      <h3>Disease Category Breakdown</h3>
      <p>32 conditions distributed across 10 clinical categories</p>
    </a>
    <a class="card" href="gene_summary.csv" target="_blank">
      <div class="card-icon">📋</div>
      <h3>Gene Summary CSV</h3>
      <p>Per-gene table: diseases, degree, pleiotropic flag, PubMed reference</p>
    </a>
  </div>
  <div style="margin-top:36px;background:#1E293B;border-radius:12px;padding:20px;border:1px solid #334155">
    <h3 style="color:#F1F5F9;margin-bottom:10px">Pleiotropic Genes (bridging ≥2 diseases)</h3>
    <p style="color:#94A3B8;font-size:13px;line-height:1.8">
      {" ".join(f'<span class="badge">{g}</span>' for g in stats.get("pleiotropic_list", []))}
    </p>
  </div>
</div>
<div class="footer">
  Data: HPO (hpo.jax.org) · OMIM (omim.org) · ClinVar · GWAS Catalog &nbsp;|&nbsp;
  Built by <strong>Rayala Madhu Bhanu Varma</strong> · Lead Data Scientist · GenepoweRx
</div>
</body>
</html>"""
    (outdir / "index.html").write_text(html, encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
