"""
src/graph_builder.py
--------------------
Builds a multi-type NetworkX knowledge graph from the HPO dataset.

Node types
----------
  disease  — one per condition (32 total)
  gene     — one per unique gene (96 total)
  phenotype — HPO terms (shared across diseases)

Edge types
----------
  gene_disease   — gene causes / is risk factor for disease
  gene_phenotype — gene is annotated with HPO phenotype term
  disease_phenotype — disease is characterised by HPO term
"""

from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import networkx as nx
from collections import defaultdict
from data.hpo_data import HPO_DATA, HPO_TERM_LABELS, CATEGORY_COLORS, EFFECT_COLORS


# ── Build graph ───────────────────────────────────────────────────────────────

def build_graph(
    include_phenotype_nodes: bool = True,
    conditions: list[str] | None = None,
) -> nx.Graph:
    """
    Build and return the full HPO knowledge graph.

    Args:
        include_phenotype_nodes: If True, add HPO phenotype term nodes.
        conditions: Subset of conditions to include (None = all 32).

    Returns:
        nx.Graph with node/edge attributes for visualisation.
    """
    G = nx.Graph()
    source = HPO_DATA
    if conditions:
        source = {k: v for k, v in HPO_DATA.items() if k in conditions}

    # ── Disease nodes ─────────────────────────────────────────────────────────
    for disease, info in source.items():
        cat   = info.get("category", "Other")
        color = CATEGORY_COLORS.get(cat, "#94A3B8")
        G.add_node(
            disease,
            node_type  = "disease",
            category   = cat,
            hpo_id     = info.get("hpo_id", ""),
            omim       = info.get("omim", ""),
            description= info.get("description", ""),
            color      = color,
            size       = 30,
            shape      = "diamond",
            label      = disease,
        )

    # ── Gene nodes + edges ────────────────────────────────────────────────────
    gene_disease_count: dict[str, int] = defaultdict(int)

    for disease, info in source.items():
        for gene, gdata in info["genes"].items():
            gene_disease_count[gene] += 1
            effect = gdata.get("effect", "risk")

            # add gene node (update size if already exists from another disease)
            if gene not in G.nodes:
                G.add_node(
                    gene,
                    node_type = "gene",
                    color     = "#64748B",
                    size      = 18,
                    shape     = "dot",
                    label     = gene,
                    diseases  = [],
                    pubmeds   = [],
                )
            G.nodes[gene]["diseases"] = G.nodes[gene].get("diseases", []) + [disease]
            if gdata.get("pubmed"):
                G.nodes[gene]["pubmeds"] = G.nodes[gene].get("pubmeds", []) + [gdata["pubmed"]]

            # gene–disease edge
            G.add_edge(
                gene, disease,
                edge_type = "gene_disease",
                effect    = effect,
                pubmed    = gdata.get("pubmed"),
                evidence  = gdata.get("evidence", ""),
                color     = EFFECT_COLORS.get(effect, "#94A3B8"),
                width     = 3 if effect == "causal" else 1.5,
            )

            # phenotype nodes + edges
            if include_phenotype_nodes:
                for hpo_id in gdata.get("hpo_terms", []):
                    label = HPO_TERM_LABELS.get(hpo_id, hpo_id)
                    if hpo_id not in G.nodes:
                        G.add_node(
                            hpo_id,
                            node_type = "phenotype",
                            label     = label,
                            color     = "#A78BFA",
                            size      = 10,
                            shape     = "dot",
                        )
                    G.add_edge(
                        gene, hpo_id,
                        edge_type = "gene_phenotype",
                        color     = "#C4B5FD",
                        width     = 0.8,
                    )
                    if disease in G.nodes and not G.has_edge(disease, hpo_id):
                        G.add_edge(
                            disease, hpo_id,
                            edge_type = "disease_phenotype",
                            color     = "#DDD6FE",
                            width     = 0.5,
                        )

    # Update gene node sizes by disease-degree (pleiotropic genes appear bigger)
    for gene in G.nodes:
        if G.nodes[gene].get("node_type") == "gene":
            n = gene_disease_count.get(gene, 1)
            G.nodes[gene]["size"] = 14 + n * 8   # bigger = more diseases
            if n > 1:
                G.nodes[gene]["color"] = "#F97316"  # orange = pleiotropic

    return G


# ── Analytics helpers ─────────────────────────────────────────────────────────

def get_graph_stats(G: nx.Graph) -> dict:
    disease_nodes  = [n for n, d in G.nodes(data=True) if d.get("node_type") == "disease"]
    gene_nodes     = [n for n, d in G.nodes(data=True) if d.get("node_type") == "gene"]
    pheno_nodes    = [n for n, d in G.nodes(data=True) if d.get("node_type") == "phenotype"]
    gene_dis_edges = [(u,v) for u,v,d in G.edges(data=True) if d.get("edge_type") == "gene_disease"]
    pleiotropic    = [g for g in gene_nodes if G.nodes[g].get("color") == "#F97316"]
    return {
        "diseases":         len(disease_nodes),
        "genes":            len(gene_nodes),
        "phenotypes":       len(pheno_nodes),
        "gene_disease_links": len(gene_dis_edges),
        "pleiotropic_genes": len(pleiotropic),
        "total_edges":      G.number_of_edges(),
        "density":          round(nx.density(G), 4),
        "components":       nx.number_connected_components(G),
        "pleiotropic_list": sorted(pleiotropic),
    }


def get_gene_summary(G: nx.Graph) -> list[dict]:
    """Return per-gene summary rows for the data table."""
    rows = []
    for n, d in G.nodes(data=True):
        if d.get("node_type") != "gene":
            continue
        diseases  = list(set(d.get("diseases", [])))
        pubmeds   = list(set(d.get("pubmeds", [])))
        degree    = G.degree(n)
        rows.append({
            "gene":       n,
            "diseases":   ", ".join(diseases),
            "n_diseases": len(diseases),
            "degree":     degree,
            "pleiotropic": len(diseases) > 1,
            "pubmed":     pubmeds[0] if pubmeds else None,
        })
    return sorted(rows, key=lambda x: -x["n_diseases"])


def get_disease_summary(G: nx.Graph) -> list[dict]:
    """Return per-disease summary rows."""
    rows = []
    for n, d in G.nodes(data=True):
        if d.get("node_type") != "disease":
            continue
        neighbors = list(G.neighbors(n))
        genes     = [x for x in neighbors if G.nodes[x].get("node_type") == "gene"]
        rows.append({
            "disease":   n,
            "category":  d.get("category", ""),
            "n_genes":   len(genes),
            "hpo_id":    d.get("hpo_id", ""),
            "omim":      d.get("omim", ""),
            "description": d.get("description", ""),
        })
    return sorted(rows, key=lambda x: -x["n_genes"])


def shortest_path_between_genes(G: nx.Graph, gene_a: str, gene_b: str) -> list | None:
    """Find shortest path between two genes through disease / phenotype nodes."""
    try:
        return nx.shortest_path(G, source=gene_a, target=gene_b)
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return None


def shared_diseases(G: nx.Graph, gene_a: str, gene_b: str) -> list[str]:
    """Return list of diseases shared by both genes."""
    da = set(G.nodes[gene_a].get("diseases", [])) if gene_a in G.nodes else set()
    db = set(G.nodes[gene_b].get("diseases", [])) if gene_b in G.nodes else set()
    return sorted(da & db)
