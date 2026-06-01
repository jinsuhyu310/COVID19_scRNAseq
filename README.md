COVID19_scRNAseq
Single-cell RNA sequencing analysis of peripheral blood mononuclear cells (PBMCs) from COVID-19 patients receiving TDR (Tocilizumab + Dexamethasone + Remdesivir) treatment, comparing immune responses between patients with poor and good prognosis.

Study Design
GroupPatientsTimepointsPoor prognosisP410, P413, P466Day 1 (pre-treatment), Day 7 (post-treatment)Good prognosisP447, P448, P468Day 1 (pre-treatment), Day 7 (post-treatment)
12 samples total (6 patients × 2 timepoints), processed with 10x Genomics Chromium.

Repository Structure
  COVID19_scRNAseq/
  └── scripts/
      ├── QC&Preprocessing.R
      ├── Annotation.R
      ├── DEG_GO_Analysis.R
      ├── CD14_Monocyte_Analysis.R
      ├── Monocyte_Analysis.R
      ├── B_Analysis.R
      └── T_Analysis.R


Dependencies
r# Core
Seurat (v5), SeuratDisk, SeuratExtend
SingleR, celldex
monocle3

Enrichment
clusterProfiler, org.Hs.eg.db, enrichplot
msigdbr, KEGGREST, AnnotationDbi

Visualization
ggplot2, ggrepel, patchwork, EnhancedVolcano
ggpubr, rstatix, pheatmap, RColorBrewer
scales, viridis

Utilities
dplyr, tidyr, tibble, stringr
Matrix, future


Data Availability
Raw data will be deposited in GEO upon publication.

Contact
Jinsuh Yu — Pusan National University (jinsuhyu310@gmail.com)
