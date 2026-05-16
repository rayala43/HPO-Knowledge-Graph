# 🧬 HPO Gene–Phenotype Knowledge Graph Explorer

A NetworkX knowledge graph linking **96 genes → 108 associations → 32 diseases** using curated HPO, OMIM, ClinVar, and GWAS Catalog data. Outputs fully interactive HTML dashboards — no server required.

Built from real-world clinical genomics experience at GenepoweRx.

---

## Quick Start

```bash
git clone https://github.com/YOUR_USERNAME/hpo-knowledge-graph
cd hpo-knowledge-graph
pip install -r requirements.txt

# Build all outputs (opens output/index.html)
python main.py

# Include HPO phenotype term nodes
python main.py --phenotypes

# Filter to specific conditions
python main.py --conditions "Type 2 Diabetes" "Alzheimer Disease" "Parkinson Disease"

# Query a gene
python main.py --gene PTPN22

# Find shortest path between two genes
python main.py --path PTPN22 STAT4

# List all 32 conditions
python main.py --list-conditions
```

Open `output/index.html` in any browser — no server needed.

---

## Outputs

| File | Description |
|---|---|
| `output/index.html` | Dashboard with KPI cards + links to all visualisations |
| `output/hpo_knowledge_graph.html` | **Interactive Pyvis network** — drag, zoom, hover |
| `output/degree_distribution.html` | Top 25 genes by connectivity (pleiotropic genes highlighted) |
| `output/genes_per_category.html` | Unique genes per disease category |
| `output/pleiotropy.html` | Genes bridging multiple diseases |
| `output/category_pie.html` | Disease distribution across 10 categories |
| `output/gene_summary.csv` | Per-gene table with disease links and PubMed IDs |
| `output/disease_summary.csv` | Per-disease table with gene count and OMIM ID |

---

## Dataset

`data/hpo_data.py` — **32 conditions · 96 genes · 108 gene–disease links**

All annotations are peer-reviewed and sourced from:
- **HPO** (hpo.jax.org) — phenotype term IDs and gene–phenotype links
- **OMIM** (omim.org) — disease identifiers
- **ClinVar** — clinical significance classifications
- **GWAS Catalog** (ebi.ac.uk/gwas) — GWAS-derived associations

### Disease Categories

| Category | Conditions |
|---|---|
| Metabolic | Type 2 Diabetes, Type 1 Diabetes, Familial Hypercholesterolemia, PKU, Wilson, Gaucher, Fabry |
| Neurological | Alzheimer, Parkinson, Huntington, ALS, Autism, Schizophrenia, NF1, Tuberous Sclerosis, Fragile X |
| Cardiovascular | Coronary Artery Disease, Hypertension |
| Oncology | Breast Cancer, Colorectal Cancer |
| Autoimmune | Rheumatoid Arthritis, Lupus (SLE), Celiac Disease |
| Hematological | Sickle Cell Disease, Hemophilia A, Hemophilia B |
| Chromosomal | Down Syndrome, Turner Syndrome, Klinefelter Syndrome |
| Connective Tissue | Marfan Syndrome |
| Muscular | Duchenne Muscular Dystrophy |
| Respiratory | Cystic Fibrosis |

### Pleiotropic Genes (bridge ≥2 diseases)

| Gene | Diseases |
|---|---|
| PTPN22 | Type 1 Diabetes · Rheumatoid Arthritis · Lupus (SLE) |
| APOE | Coronary Artery Disease · Alzheimer Disease |
| PCSK9 | Coronary Artery Disease · Familial Hypercholesterolemia |
| LDLR | Coronary Artery Disease · Familial Hypercholesterolemia |
| APOB | Coronary Artery Disease · Familial Hypercholesterolemia |
| STAT4 | Rheumatoid Arthritis · Lupus (SLE) |
| GBA | Parkinson Disease · Gaucher Disease |
| HLA-DRB1 | Type 1 Diabetes · Rheumatoid Arthritis |
| FMR1 | Autism Spectrum Disorder · Fragile X Syndrome |
| APP | Alzheimer Disease · Down Syndrome |
| SOD1 | ALS · Down Syndrome |

---

## Project Structure

```
hpo_explorer/
├── main.py                    # CLI entry point
├── src/
│   ├── graph_builder.py       # NetworkX graph construction + analytics
│   └── visualiser.py          # Pyvis HTML + Plotly chart builders
├── data/
│   └── hpo_data.py            # Curated HPO dataset (32 conditions)
├── output/                    # Generated HTML/CSV files land here
├── requirements.txt
└── README.md
```

---

## Graph Schema

```
Node types
──────────
  disease   (32)  — HPO/OMIM condition  [diamond shape]
  gene      (96)  — HGNC symbol         [circle; orange = pleiotropic]
  phenotype (N)   — HPO term (optional) [small dot; --phenotypes flag]

Edge types
──────────
  gene_disease    — gene linked to disease (effect: causal / risk / protective)
  gene_phenotype  — gene annotated with HPO term
  disease_phenotype — disease characterised by HPO term
```

---

## CLI Reference

```
python main.py [OPTIONS]

Options:
  --phenotypes           Include HPO phenotype term nodes
  --conditions NAME ...  Filter to specific conditions
  --no-physics           Disable force-directed layout (faster)
  --gene GENE            Query a specific gene
  --path GENE_A GENE_B   Find shortest path between two genes
  --output-dir DIR       Output directory (default: output/)
  --list-conditions      List all 32 conditions with gene counts
```

---

## Tech Stack

| Layer | Library |
|---|---|
| Graph engine | NetworkX |
| Interactive network | Pyvis |
| Charts | Plotly |
| Data | pandas, numpy |

---

## License
MIT
