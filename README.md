# Proteomics Analysis of the 5xFAD Alzheimer's Mouse Model

## Publication

> **Transcending the amyloid-beta dominance paradigm in Alzheimer's disease: An exploration of behavioural, metabolic, and gut microbiota phenotypes in 5xFAD mice**
>
> Medina-Vera D., Zambrana-Infantes E.N., López-Gambero A.J., Verheul-Campos J., Santín L.J., Baixeras E., Suarez J., Pavon F.J., Rosell-Valle C., Rodríguez de Fonseca F.

---

## Project Overview

This repository contains the R scripts used for the **shotgun proteomics analysis** of the 5xFAD mouse model of Alzheimer's disease, comparing transgenic animals against wild-type (WT) controls.

The **5xFAD** mouse carries five familial AD mutations (*APP Swedish K670N/M671L*, *APP Florida I716V*, *APP London V717I*, *PSEN1 M146L*, *PSEN1 L286V*) and is one of the most aggressive amyloid models available, developing intraneuronal Aβ42 accumulation and plaque deposition as early as 2 months of age.

The pipeline covers the full workflow from raw LFQ (Label-Free Quantification) abundance data exported from Proteome Discoverer to statistical analysis, visualization, and protein network/enrichment analysis.

---

## Key Differences vs 3×Tg-AD Pipeline

| Parameter | 3×Tg-AD | 5xFAD |
|-----------|---------|-------|
| p-value threshold | 0.05 | **0.01** |
| FC threshold | 25% (log₂ ±0.32/−0.41) | **50% (log₂ ±0.58/−1.00)** |
| Coverage filter (5xFAD group) | ≥3 valid values | **≥4 males / ≥2 females** |
| Input file | 230414-LFQ_Filtradas.xlsx | **230223-LFQ_Incl-Isoforms_FILTRADAS.xlsx** |

---

## Results Summary

| Comparison | Proteins detected | Up-regulated | Down-regulated |
|------------|:-----------------:|:------------:|:--------------:|
| 5xFAD vs WT | 2,784 | 417 | 111 |

---

## Repository Structure

```
proteomics-5xFAD/
│
├── R/
│   ├── preprocessing/
│   │   └── 01_preprocessing.R            # Data loading, cleaning, imputation & normalization
│   ├── analysis/
│   │   └── 02_differential_expression.R  # Statistical tests, log2FC, volcano table
│   └── visualization/
│       └── 03_visualization.R            # Volcano, heatmap, PCA, network, enrichment plots
│
├── data/
│   ├── raw/                              # Input files (not tracked — see below)
│   └── processed/                        # Intermediate outputs
│
├── results/
│   ├── figures/                          # All generated figures (tif, svg, pdf, png)
│   └── tables/                           # Differentially expressed proteins, full stats table
│
├── docs/
│   └── pipeline_overview.md              # Step-by-step methodology description
│
├── .gitignore
└── README.md
```

---

## Analysis Pipeline

### 1. Preprocessing (`01_preprocessing.R`)
- Loads raw LFQ abundance data from Proteome Discoverer output (includes isoforms)
- Fills missing gene symbols using a manually curated dictionary (21 entries)
- Separates samples by genotype (WT / 5xFAD) and sex (male / female)
- **Outlier removal** per protein per group using a 3×IQR rule
- Coverage filter: ≥3 valid values (WT groups), ≥4 (5xFAD males), ≥2 (5xFAD females)
- Retains only proteins detected across all four sex × genotype groups
- **Median imputation** within groups
- **Normalization** by a down-scaling correction factor (column mean / minimum column mean)

### 2. Differential Expression (`02_differential_expression.R`)
- Compares 5xFAD vs WT across all animals
- Normality assessed per protein with **Shapiro-Wilk test**:
  - Both groups normal → **Welch's t-test**
  - At least one non-normal → **Wilcoxon rank-sum test**
- Computes log₂(Fold Change) and −log₁₀(p-value)
- Classifies proteins as **up-regulated** or **down-regulated**
  - Thresholds: **p < 0.01** and **|log₂FC| > log₂(1.5)** [±0.58]

### 3. Visualization (`03_visualization.R`)
- **PCA** of normalized abundances (colored by genotype, shaped by sex)
- **Volcano plot** with labeled top candidates
- **Heatmap** of differentially expressed proteins (log₁₀ abundances, hierarchical clustering)
- **Protein–protein interaction network** (STRING-derived edges, ggraph / Fruchterman–Reingold layout)
- **GO and KEGG enrichment** bar plots (top 20 terms per category, FDR-colored)

All figures saved in TIFF (600 dpi), SVG, PDF, and PNG formats.

---

## Input Data

| File | Location | Description |
|------|----------|-------------|
| `230223-LFQ_Incl-Isoforms_FILTRADAS.xlsx` | `data/raw/` | LFQ protein abundances from Proteome Discoverer (includes isoforms) |
| `sample-metadata.xlsx` | `data/raw/` | Sample annotations: ID, Genotype, sex |
| `Enrichment_Results_5xfad_005.xlsx` | `data/processed/` | GO/KEGG enrichment results from STRING |
| `RED_5xfad_005.xlsx` | `data/processed/` | Protein–protein interaction edges from STRING |

> ⚠️ Raw data files are not tracked in this repository. Available upon reasonable request from the corresponding author: [fernando.rodriguez@ibima.eu](mailto:fernando.rodriguez@ibima.eu)

---

## Dependencies

```r
install.packages(c("readxl", "openxlsx", "dplyr", "tidyverse", "data.table",
                   "effsize", "car", "ggplot2", "ggrepel", "pheatmap",
                   "RColorBrewer", "gplots", "FactoMineR", "factoextra",
                   "missMDA", "VIM", "igraph", "ggraph", "edgeR",
                   "tidytext", "tm", "formattable"))
```

---

## Usage

```r
setwd("path/to/proteomics-5xFAD")

source("R/preprocessing/01_preprocessing.R")
source("R/analysis/02_differential_expression.R")
source("R/visualization/03_visualization.R")
```

---

## Authors

| Role | Name | Affiliation |
|------|------|-------------|
| Proteomics analysis | **Julia Verheul-Campos** | IBIMA Plataforma BIONAND, Málaga |
| Study design & coordination | D. Medina-Vera, F.J. Pavon, F. Rodríguez de Fonseca | IBIMA, Málaga |
| Animal experiments | E.N. Zambrana-Infantes, A.J. López-Gambero, C. Rosell-Valle | IBIMA, Málaga |
| Collaborators | L.J. Santín, E. Baixeras, J. Suarez | Universidad de Málaga |

**Correspondence:** [fernando.rodriguez@ibima.eu](mailto:fernando.rodriguez@ibima.eu)

---

## Citation

```bibtex
@article{MedinaVera2025_5xFAD,
  author    = {Medina-Vera, D. and Zambrana-Infantes, E.N. and López-Gambero, A.J.
               and Verheul-Campos, J. and Santín, L.J. and Baixeras, E.
               and Suarez, J. and Pavon, F.J. and Rosell-Valle, C.
               and Rodríguez de Fonseca, F.},
  title     = {Transcending the amyloid-beta dominance paradigm in {Alzheimer's} disease:
               An exploration of behavioural, metabolic, and gut microbiota phenotypes
               in {5xFAD} mice},
  journal   = {(journal pending)},
  year      = {2025}
}
```

---

## License

MIT License — see `LICENSE` for details.
