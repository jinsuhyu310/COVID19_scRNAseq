# COVID19_scRNAseq

Single-cell RNA sequencing analysis of peripheral blood mononuclear cells (PBMCs) from COVID-19 patients receiving TDR (Tocilizumab + Dexamethasone + Remdesivir) treatment, comparing immune responses between patients with poor and good prognosis.

---

## Study Design

| Group | Patients | Timepoints |
|-------|----------|------------|
| Poor prognosis | P410, P413, P466 | Day 1 (pre-treatment), Day 7 (post-treatment) |
| Good prognosis | P447, P448, P468 | Day 1 (pre-treatment), Day 7 (post-treatment) |

12 samples total (6 patients × 2 timepoints), processed with 10x Genomics Chromium.

---

## Repository Structure

```
COVID19_scRNAseq/
└── scripts/
    ├── QC&Preprocessing.R
    ├── Annotation.R
    ├── DEG_GO_Analysis.R
    ├── CD14_Monocyte_Analysis.R
    ├── Monocyte_Analysis.R
    ├── B_Analysis.R
    └── T_Analysis.R
```

---

## Scripts

Raw 10x `.h5` files are loaded, QC-filtered (`nFeature_RNA 200–6000`, `percent.mt < 30`), and integrated across 12 samples using CCA (`QC&Preprocessing.R`). Cell types are annotated with SingleR (Monaco Immune reference) and Azimuth (PBMC multimodal reference), followed by manual refinement to produce major and fine-grained labels (`Annotation.R`). Whole-PBMC DEG analysis compares pre- and post-treatment within each prognosis group using Wilcoxon tests, with results visualized as volcano plots, a quadrant DEG plot (log2FC Poor vs Good), GO enrichment barplots, and a cross-cell-type heatmap (`DEG_GO_Analysis.R`). CD14 Monocytes are subclustered and Cluster 4 (Poor_d7) is characterized by DEG, GO, and GSEA KEGG analysis alongside S100A8/A9/A12 expression (`CD14_Monocyte_Analysis.R`). CD14 and CD16 Monocytes are re-integrated with SCT and analyzed for IL6+ cell distribution, Hallmark/KEGG module scores (IL6-JAK-STAT3, TNFα-NFκB, Glycolysis, Hypoxia), and temporal changes (Δ d7−d1) per cluster, with detailed characterization of Cluster 7 (`Monocyte_Analysis.R`). B cell subsets (naive B, memory B, Plasmablast) are examined for proportional changes, DEG (including memory B vs all), GO enrichment, and TNF/IL6 signaling module scores (`B_Analysis.R`). CD8 T cells are annotated by module score (Naive-like, Early Activated, Effector Memory, Exhausted), trajectory-mapped with Monocle3, and analyzed for cluster proportions, DEG, TEM gene expression patterns, ROS-related genes, and cytotoxicity/exhaustion scores across CD8 TCM/TEM and CD4 CTL (`T_Analysis.R`).

---

## Dependencies

```r
# Core
Seurat (v5), SeuratDisk, SeuratExtend
SingleR, celldex
monocle3

# Enrichment
clusterProfiler, org.Hs.eg.db, enrichplot
msigdbr, KEGGREST, AnnotationDbi

# Visualization
ggplot2, ggrepel, patchwork, EnhancedVolcano
ggpubr, rstatix, pheatmap, RColorBrewer
scales, viridis

# Utilities
dplyr, tidyr, tibble, stringr
Matrix, future
```

---

## Data Availability

Raw data will be deposited in GEO upon publication.

---

## Contact

Jinsuh Yu — Pusan National University (jinsuhyu310@gmail.com)
