##################################################################### SingleR annotation ####################################

setwd("~/covid_scRNA/covid_R/NEW/Annotation")
merged_int <-readRDS("~/covid_scRNA/covid_R/NEW/merged_integration_all.rds")
DefaultAssay(merged_int) <- "RNA"
merged_int <- NormalizeData(merged_int, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
merged_int <- JoinLayers(merged_int)

library(SingleR)
library(celldex)

ref <- celldex::MonacoImmuneData()
test_mat <- GetAssayData(merged_int, assay = "RNA", layer = "data")


pred_cells <- SingleR(
  test = as.matrix(test_mat),
  ref = ref,
  labels = ref$label.fine
)

merged_int$SingleR_label <- pred_cells$labels

Idents(merged_int) <- "seurat_clusters"
pred_clusters <- SingleR(
  test = as.matrix(test_mat),
  ref = ref,
  labels = ref$label.fine,
  clusters = Idents(merged_int)
)

cluster_map <- setNames(pred_clusters$labels, rownames(pred_clusters))
cluster_ids <- as.character(Idents(merged_int))
singler_cluster_labels <- cluster_map[cluster_ids]
names(singler_cluster_labels) <- colnames(merged_int)

merged_int <- AddMetaData(merged_int, metadata = singler_cluster_labels, col.name = "SingleR_cluster")

table(merged_int$SingleR_label)
sort(table(merged_int$SingleR_label), decreasing = TRUE) %>% head(30)

merged_int@meta.data <- merged_int@meta.data %>%
  mutate(
    SingleR_label_simple = case_when(
      SingleR_label %in% c(
        "Classical monocytes" 
      ) ~ "CD14 Mono",
      SingleR_label %in% c(
        "Non classical monocytes"
      ) ~ "CD16 Mono",
      SingleR_label %in% c(
        "Intermediate monocytes"
      ) ~ "Intermediate Mono",
      SingleR_label %in% c(
        "Naive CD4 T cells", "Terminal effector CD4 T cells", "Th1 cells", "Th2 cells", 
        "Th17 cells", "Th1/Th17 cells", "Follicular helper T cells", "T regulatory cells"
      ) ~ "CD4 T",
      
      SingleR_label %in% c(
        "Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", 
        "Terminal effector CD8 T cells"
      ) ~ "CD8 T",
      
      SingleR_label %in% c("MAIT cells", "Vd2 gd T cells", "Non-Vd2 gd T cells") ~ "gd T cells",
      
      SingleR_label %in% c(
        "Naive B cells", "Switched memory B cells", "Non-switched memory B cells",
        "Exhausted B cells", "Plasmablasts"
      ) ~ "B",
      
      SingleR_label %in% c("Natural killer cells") ~ "NK",
      
      SingleR_label %in% c("Myeloid dendritic cells", "Plasmacytoid dendritic cells") ~ "DC",
      
      SingleR_label %in% c("Low-density basophils", "Low-density neutrophils") ~ "Granulocytes",
      
      SingleR_label %in% c("Progenitor cells") ~ "Progenitors",
      TRUE ~ "Other"
    )
  )

table(merged_int$SingleR_label_simple)

my_color <- c(
  "B"="#33A02C","CD14 Mono"="#FD7F27","CD16 Mono"="#FDBF6F",
  "CD4 T"="#1F78B4","CD8 T"="#E31A1C","NK"="#4ECCCB",
  "DC"="#B15928","gd T cells"="#999999","Granulocytes"="#FF7F00",
  "Progenitors"="#CAB2D6","Intermediate Mono"="#6A3D9A","Other"="black"
)

p_singler <- DimPlot(
  merged_int,
  reduction = "umap",        
  group.by = "SingleR_label_simple",
  label = TRUE,
  repel = TRUE, label.size = 4
) + 
  scale_color_manual(values = my_color) + 
  ggtitle("SingleR Annotation") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

dir.create("./png/", showWarnings = FALSE, recursive = TRUE)
ggsave("./png/UMAP_SingleR_Annotation.png", p_singler, width = 8, height = 6, dpi = 400)

sample_order <- c("P410-d1","P410-d7","P413-d1","P413-d7",
                  "P466-d1","P466-d7","P447-d1","P447-d7",
                  "P448-d1","P448-d7","P468-d1","P468-d7")
merged_int$orig.ident <- factor(merged_int$orig.ident, levels = sample_order)

p1 <- DimPlot(
  merged_int, reduction = "umap",
  group.by = "SingleR_label_simple",
  split.by = "orig.ident",
  label = TRUE, repel = TRUE, label.size = 3,
  ncol = 6
) + scale_color_manual(values = my_color) +
  ggtitle("SingleR Annotation")

ggsave("./png/UMAP_SingleR_Annotation_split.png", p1, width = 18, height = 8, dpi = 300)


##################################################################### Azimuth annotation #######################################

reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")

DefaultAssay(reference) <- "SCT"
DefaultAssay(merged_int) <- "RNA"
merged_int <- SCTransform(merged_int, verbose = FALSE)
DefaultAssay(merged_int) <- "SCT"

transfer_features <- intersect(VariableFeatures(reference), VariableFeatures(merged_int))

anchors <- FindTransferAnchors(
  reference = reference,
  query = merged_int,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:20,
  recompute.residuals = FALSE
)

merged_int <- MapQuery(
  anchorset = anchors,
  query = merged_int,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap")

View(merged_int@meta.data)
merged_int$predicted.celltype.major <- merged_int$predicted.celltype.l2
merged_int$predicted.celltype.major[merged_int$predicted.celltype.l2 %in%
                                      c("B intermediate", "B memory", "B naive", "Plasmablast")] <- "B"

merged_int$predicted.celltype.major[merged_int$predicted.celltype.l2 %in% c("CD14 Mono")] <- "CD14 Mono"
merged_int$predicted.celltype.major[merged_int$predicted.celltype.l2 %in% c("CD16 Mono")] <- "CD16 Mono"
merged_int$predicted.celltype.major[merged_int$predicted.celltype.l2 %in% 
                                      c("CD4 Naive","CD4 TEM","CD4 TCM","CD4 CTL","CD4 Proliferating","Treg")] <- "CD4 T"
merged_int$predicted.celltype.major[merged_int$predicted.celltype.l2 %in% 
                                      c("CD8 Naive","CD8 TEM","CD8 TCM","CD8 Proliferating")] <- "CD8 T"
merged_int$predicted.celltype.major[merged_int$predicted.celltype.l2 %in%
                                      c("NK","NK_CD56bright","NK Proliferating")] <- "NK"
merged_int$predicted.celltype.major[merged_int$predicted.celltype.l2 %in%
                                      c("cDC2","pDC")] <- "DC"
merged_int$predicted.celltype.major[merged_int$predicted.celltype.l2 %in%
                                      c("dnT","gdT","MAIT","ILC")] <- "other T"
merged_int$predicted.celltype.major[merged_int$predicted.celltype.l2 %in% c("HSPC")] <- "HSPC"

my_colorss <- c(
  "B"="#33A02C","Mono"="#FD7F27",
  "CD4 T"="#1F78B4","CD8 T"="#E31A1C","NK"="#4ECCCB",
  "DC"="#B15928","other"="#999999",
  "other T"="#CAB2D6"
)

p_singler <- DimPlot(
  merged_int,
  reduction = "umap",        
  group.by = "predicted.celltype.l1",
  label = TRUE,
  repel = TRUE, label.size = 4
) + 
  scale_color_manual(values = my_colorss) + 
  ggtitle("Azimuth Annotation") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

dir.create("./png/", showWarnings = FALSE, recursive = TRUE)
ggsave("./png/UMAP_Azimuth_Annotation.png", p_singler, width = 7, height = 6, dpi = 400)


############################################################ Azimuth + SingleR overwrite annotation ##################################

stopifnot("predicted.celltype.major" %in% colnames(merged_int@meta.data))
stopifnot("SingleR_label_simple"   %in% colnames(merged_int@meta.data))

## 1) Doublet → SingleR_label_simple
doublet_idx <- which(merged_int$predicted.celltype.major == "Doublet")
valid_doublet_idx <- doublet_idx[!is.na(merged_int$SingleR_label_simple[doublet_idx])]

n_overwrite1 <- sum(
  merged_int$predicted.celltype.major[valid_doublet_idx] != 
    merged_int$SingleR_label_simple[valid_doublet_idx],
  na.rm = TRUE
)
merged_int$predicted.celltype.major[valid_doublet_idx] <- 
  merged_int$SingleR_label_simple[valid_doublet_idx]

cat("1) Number of cells overwritten from Doublet to SingleR_label_simple:", n_overwrite1, "\n")

## 2) Overwrite with SingleR monocyte lineage
mono_like <- c("CD14 Mono", "CD16 Mono", "Intermediate Mono")

mono_idx <- which(merged_int$SingleR_label_simple %in% mono_like)

n_overwrite2 <- sum(
  merged_int$predicted.celltype.major[mono_idx] != 
    merged_int$SingleR_label_simple[mono_idx],
  na.rm = TRUE
)

merged_int$predicted.celltype.major[mono_idx] <- 
  merged_int$SingleR_label_simple[mono_idx]

cat("2) Number of cells overwritten with SingleR monocyte lineage:", n_overwrite2, "\n")

## 3) Remove Platelet / Eryth / Granulocytes / Progenitors
remove_targets <- c("Platelet","Eryth","Granulocytes","Progenitors")

keep_cells <- rownames(merged_int@meta.data)[
  !(merged_int$predicted.celltype.major %in% remove_targets)
]
n_removed <- ncol(merged_int) - length(keep_cells)

merged_int <- subset(merged_int, cells = keep_cells)

cat("3) Total number of removed cells (Platelet/Eryth/Granulocytes/Progenitors):", n_removed, "\n")

## 4) Reassign Intermediate Mono (all to CD16 Mono)

cd14_markers <- c("CD14", "LYZ", "S100A8", "S100A9", "VCAN", "FCN1")
p_cd14 <- FeaturePlot(
  merged_int,
  features  = cd14_markers,
  reduction = "umap",
  cols      = c("lightgrey", "red"), ncol=3,
  order     = TRUE
)

p_cd14 <- p_cd14 & theme(
  plot.title = element_text(size = 9, face = "bold"),  
  axis.title = element_text(size = 8),
  axis.text  = element_text(size = 7)
)

p_cd14 <- p_cd14 + plot_annotation(
  title = "CD14 marker expression",
  theme = theme(plot.title = element_text(size = 10, face = "bold"))
)

ggsave(
  filename = file.path(outdir, "FeaturePlot_CD14_markers.png"),
  plot     = p_cd14,
  width    = 9,
  height   = 5,
  dpi      = 300
)

cd16_markers <- c("FCGR3A", "MS4A7", "LST1", "CX3CR1", "MSR1", "IFITM3")
p_cd16 <- FeaturePlot(
  merged_int,
  features  = cd16_markers,
  reduction = "umap",
  cols      = c("lightgrey", "red"), ncol=3,
  order     = TRUE
)
p_cd16 <- p_cd16 & theme(
  plot.title = element_text(size = 9, face = "bold"),
  axis.title = element_text(size = 8),
  axis.text  = element_text(size = 7)
)

p_cd16 <- p_cd16 + plot_annotation(
  title = "CD16 marker expression",
  theme = theme(plot.title = element_text(size = 10, face = "bold"))
)

ggsave(
  filename = file.path(outdir, "FeaturePlot_CD16_markers.png"),
  plot     = p_cd16,
  width    = 9,
  height   = 5,
  dpi      = 300
)


mono_levels <- c("CD14 Mono", "CD16 Mono", "Intermediate Mono")

merged_int$Mono_group <- "Other"
merged_int$Mono_group[merged_int$predicted.celltype.major %in% mono_levels] <-
  merged_int$predicted.celltype.major[merged_int$predicted.celltype.major %in% mono_levels]

merged_int$Mono_group <- factor(merged_int$Mono_group,
                                levels = c(mono_levels, "Other"))
p_mono_map <- DimPlot(
  merged_int,
  group.by = "Mono_group",
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) + ggtitle("Monocyte UMAP")

dir.create("~/covid_scRNA/covid_R/NEW/Annotation/png", showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = "~/covid_scRNA/covid_R/NEW/Annotation/png/UMAP_Mono_group.png",
  plot     = p_mono_map,
  width    = 7,
  height   = 6,
  dpi      = 300
)


intermediate_cells <- colnames(merged_int)[
  merged_int$predicted.celltype.major == "Intermediate Mono"
]

cat("4) Number of cells relabeled from Intermediate Mono to CD16 Mono:", length(intermediate_cells), "\n")

merged_int$predicted.celltype.major[intermediate_cells] <- "CD16 Mono"

tab_after <- table(merged_int$predicted.celltype.major, useNA = "ifany")
cat("\n[After all reassign + filtering]\n")
print(tab_after)

DimPlot(
  merged_int,
  group.by  = "predicted.celltype.major",
  reduction = "umap",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("Final predicted.celltype.major")

########################################################################## re-clustering ############################################## 

DefaultAssay(merged_int) <- "integrated"   

set.seed(1234)

merged_int <- FindNeighbors(merged_int, dims = 1:15)
merged_int <- FindClusters(merged_int, resolution = 0.6)
merged_int <- RunUMAP(merged_int, dims = 1:15)

DimPlot(
  merged_int,
  group.by  = "predicted.celltype.major",  
  reduction = "umap",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("Final predicted.celltype.major (reclustered)")

saveRDS(merged_int, "~/covid_scRNA/covid_R/NEW/merged_int_annotated_final.rds")


########################################################################### UMAP ###############################################

my_colors <-c("B"="#33A02C", "CD14 Mono"="#FD7F27", "CD16 Mono"="#DAA520",
              "CD4 T"="#1F78B4", "CD8 T"="#E31A1C", "NK"="#4ECCCB", "HSPC"="#984EA3",
              "DC"="#B15928", "other T"="#CAB2D6")

DefaultAssay(merged_int)

Markers <- FeaturePlot(
  object = merged_int,
  features = c("CD14","FCGR3A","MS4A7","LYZ","CD3E","CD3G","NKG7","CD79A","CD4","CD8A"),
  reduction = "umap",
  pt.size = 0.3, min.cutoff = "q10",
  max.cutoff = "q90", ncol=5,  cols= c("lightgrey", "blue"),
)
png("./png/DEG_Markers.png", width = 5500, height = 2000, res = 300)
print(Markers)
dev.off()

###################################

predicted.celltype.major <- DimPlot(
  merged_int,
  reduction = "umap",
  group.by = "predicted.celltype.major",
  label = TRUE,
  label.size = 4,
  repel = TRUE,
  raster = FALSE,
  cols = my_colors
) + ggtitle ("") +theme(plot.title = element_text(hjust =0.5, face= "bold"))
png("./png/UMAP.png", width=4300, height=3500, res=600)
print(predicted.celltype.major)
dev.off()

merged_int$Prognosis_Group <- factor(
  merged_int$Prognosis_Group,
  levels = c("Poor_d1", "Poor_d7", "Good_d1", "Good_d7")
)

predicted.celltype.major <- DimPlot(
  merged_int,
  reduction = "umap",
  group.by = "predicted.celltype.major",
  split.by = "Prognosis_Group",
  label = TRUE,
  label.size = 3,
  repel = TRUE,
  raster = FALSE,
  cols = my_colors,
  ncol = 2
)+ theme(
  plot.title   = element_blank())
png("./png/UMAP_prognosis_group.png", width = 3000, height = 2900, res = 300)
print(predicted.celltype.major)
dev.off()


predicted.celltype.major <- DimPlot(
  merged_int,
  reduction = "umap",
  group.by = "predicted.celltype.major",
  split.by = "Group",
  label = TRUE,
  label.size = 3,
  repel = TRUE,
  raster = FALSE,
  cols = my_colors, ncol=2
)
png("./png/UMAP_group.png", width=4500, height=3000, res=600)
print(predicted.celltype.major)
dev.off()

##################################################### cell annotation fine tune

my_colors_fine <- c(
  # Existing major lineage
  "B"          = "#33A02C",  # As is
  "CD14 Mono"  = "#FD7F27",  # As is
  "CD16 Mono"  = "#DAA520",  # Similar gold
  "DC"         = "#B15928",  # As is
  "HSPC"       = "#984EA3",  # As is
  "NK"         = "#4ECCCB",  # As is
  "other T"    = "#CAB2D6",  # As is
  
  # CD4 T subtypes – Blue lineage
  "CD4 Naive"        = "#6BAED6",
  "CD4 TCM"          = "#3182BD",
  "CD4 TEM"          = "#08519C",
  "CD4 CTL"          = "#4292C6",
  "CD4 Proliferating"= "#08306B",
  
  # CD8 T subtypes – Red lineage
  "CD8 Naive"        = "#FCAE91",
  "CD8 TCM"          = "#FB6A4A",
  "CD8 TEM"          = "#E31A1C",  # Original CD8 T color
  "CD8 Proliferating"= "#A50F15",
  
  # Treg – Distinctly visible
  "Treg"             = "#FFD92F"
)


merged_int$predicted.celltype.major_fine <- merged_int$predicted.celltype.major

sel_cd4_cd8 <- merged_int$predicted.celltype.major %in% c("CD4 T", "CD8 T")

merged_int$predicted.celltype.major_fine[sel_cd4_cd8] <- 
  as.character(merged_int$predicted.celltype.l2[sel_cd4_cd8])
merged_int$predicted.celltype.major_fine <- factor(merged_int$predicted.celltype.major_fine)

merged_int$predicted.celltype.major_fine <- as.character(merged_int$predicted.celltype.major_fine)

idx_B <- merged_int$predicted.celltype.major_fine == "B"

merged_int$predicted.celltype.major_fine[idx_B] <- as.character(
  merged_int$predicted.celltype.l2[idx_B]
)

merged_int$predicted.celltype.major_fine <- factor(merged_int$predicted.celltype.major_fine)

table(merged_int$predicted.celltype.major_fine)

p_umap <- DimPlot(
  merged_int,
  reduction = "umap",
  group.by  = "predicted.celltype.major_fine",
  label     = FALSE,          # ← Turn off the original label=TRUE
  repel     = FALSE,
  raster    = FALSE,
  cols      = my_colors_fine
) +
  ggtitle("") +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"   # ← Remove legend
  )

# UMAP coordinates + celltype information
umap_df <- as.data.frame(Embeddings(merged_int, "umap"))
umap_df$predicted.celltype.major_fine <- merged_int$predicted.celltype.major_fine

# Auto-check axis names (usually "umap_1", "umap_2")
xcol <- colnames(umap_df)[1]
ycol <- colnames(umap_df)[2]

# Center coordinates of each celltype (using median)
label_df <- umap_df %>%
  dplyr::group_by(predicted.celltype.major_fine) %>%
  dplyr::summarise(
    x = median(.data[[xcol]], na.rm = TRUE),
    y = median(.data[[ycol]], na.rm = TRUE),
    .groups = "drop"
  )

# Add box labels
p_umap_labeled <- p_umap +
  geom_label_repel(
    data        = label_df,
    aes(
      x     = x,
      y     = y,
      label = predicted.celltype.major_fine,
      fill  = predicted.celltype.major_fine 
    ),
    size          = 3.3,
    label.size    = 0.25,                  
    label.padding = unit(0.2, "lines"),
    color         = "black",                
    fontface      = "bold",
    max.overlaps  = Inf,
    inherit.aes   = FALSE
  ) +
  scale_fill_manual(values = my_colors_fine)    # ← Match with UMAP colors

png("./png/UMAP_fine.png", width = 3700, height = 3650, res = 600)
print(p_umap_labeled)
dev.off()

merged_int$Prognosis_Group <- factor(
  merged_int$Prognosis_Group,
  levels = c("Poor_d1", "Poor_d7", "Good_d1", "Good_d7")
)

# 2) Plot UMAP
predicted.celltype.major <- DimPlot(
  merged_int,
  reduction = "umap",
  group.by  = "predicted.celltype.major_fine",
  split.by  = "Prognosis_Group",
  label     = FALSE,
  label.size = 3,
  repel     = TRUE,
  raster    = FALSE,
  cols      = my_colors_fine,
  ncol      = 2        # 2 x 2 array in order of Poor_d1, Poor_d7, Good_d1, Good_d7
) + 
  theme(
    plot.title = element_blank()
  )

png("./png/UMAP_prognosis_group_fine.png",
    width = 3000, height = 2700, res = 300)
print(predicted.celltype.major)
dev.off()

merged_int$Prognosis_Group <- factor(
  merged_int$Prognosis_Group,
  levels = c("Poor_d1", "Poor_d7", "Good_d1", "Good_d7")
)

## 0) Safety check (optional)
stopifnot(
  "predicted.celltype.major"      %in% colnames(merged_int@meta.data),
  "predicted.celltype.major_fine" %in% colnames(merged_int@meta.data)
)

## 1) Default everything to "fine" version
merged_int$predicted.celltype.major_nofine <- as.character(
  merged_int$predicted.celltype.major_fine
)

## 2) Overwrite with coarse "B" only for B lineage cells
b_mask <- merged_int$predicted.celltype.major == "B"

merged_int$predicted.celltype.major_nofine[b_mask] <- "B"

## 3) Organize factor levels (in desired order; can be omitted if not needed)
merged_int$predicted.celltype.major_nofine <- factor(
  merged_int$predicted.celltype.major_nofine,
  levels = c(
    "B",                     # B is a coarse label
    "CD14 Mono",
    "CD16 Mono",
    "CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM",
    "CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM",
    "DC", "HSPC", "NK",
    "other T",
    "Plasmablast", "B intermediate", "B memory", "B naive", # Add fine B at the end if necessary
    "Treg"
  )
)

## 4) Verification
table(merged_int$predicted.celltype.major)
table(merged_int$predicted.celltype.major_fine)
table(merged_int$predicted.celltype.major_nofine)

# 2) Plot UMAP
predicted.celltype.major <- DimPlot(
  merged_int,
  reduction = "ref.umap",
  group.by  = "predicted.celltype.major_nofine",
  split.by  = "Prognosis_Group",
  label     = FALSE,
  label.size = 5,
  repel     = TRUE,
  raster    = FALSE,
  cols      = my_colors_fine,
  ncol      = 4     
) + 
  theme(
    plot.title = element_blank(),
    strip.text.x = element_text(size = 13.5, face = "bold")
  )

png("./png/refUMAP_prognosis_group_fine.png",
    width = 9000, height = 2500, res = 600)
print(predicted.celltype.major)
dev.off()

predicted.celltype.major <- DimPlot(
  merged_int,
  reduction = "umap",
  group.by = "predicted.celltype.major_fine",
  split.by = "orig.ident",
  label = FALSE,
  label.size = 3,
  repel = TRUE,
  raster = FALSE,
  cols = my_colors_fine, ncol=6
)+ggtitle("")
png("./png/UMAP_indivudual_fine.png", width=9000, height=4500, res=600)
print(predicted.celltype.major)
dev.off()

predicted.celltype.major <- DimPlot(
  merged_int,
  reduction = "ref.umap",
  group.by = "predicted.celltype.major_fine",
  split.by = "orig.ident",
  label = FALSE,
  label.size = 3,
  repel = TRUE,
  raster = FALSE,
  cols = my_colors_fine, ncol=6
)+ggtitle("")
png("./png/refUMAP_indivudual_fine.png", width=10000, height=4000, res=600)
print(predicted.celltype.major)
dev.off()

############################################ ref.umap (fine version)

library(dplyr)
library(ggrepel)

## 1) Basic UMAP (ref.umap) plot
p_umap <- DimPlot(
  merged_int,
  reduction = "ref.umap",
  group.by  = "predicted.celltype.major_fine",
  label     = FALSE,
  repel     = FALSE,
  raster    = FALSE,
  cols      = my_colors_fine
) +
  ggtitle("") +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

## 2) ref.umap coordinates + celltype information
ref_df <- as.data.frame(Embeddings(merged_int, "ref.umap"))
ref_df$predicted.celltype.major_fine <- merged_int$predicted.celltype.major_fine

# Axis names (e.g., "refUMAP_1", "refUMAP_2")
xcol <- colnames(ref_df)[1]
ycol <- colnames(ref_df)[2]

## 3) Center coordinates of each celltype (median) - without across
label_df <- ref_df %>%
  dplyr::group_by(predicted.celltype.major_fine) %>%
  dplyr::summarise(
    x = median(.data[[xcol]], na.rm = TRUE),
    y = median(.data[[ycol]], na.rm = TRUE),
    .groups = "drop"
  )

## 4) Add label boxes (with celltype color as background)
p_umap_labeled <- p_umap +
  geom_label_repel(
    data = label_df,
    aes(
      x     = x,
      y     = y,
      label = predicted.celltype.major_fine,
      fill  = predicted.celltype.major_fine
    ),
    size          = 3.3,
    label.size    = 0.25,
    label.padding = unit(0.2, "lines"),
    color         = "black",
    fontface      = "bold",
    max.overlaps  = Inf,
    inherit.aes   = FALSE
  ) +
  scale_fill_manual(values = my_colors_fine)

png("./png/refUMAP_fine.png", width = 3700, height = 3650, res = 600)
print(p_umap_labeled)
dev.off()




###################### umap_individual
sample_order <- c(
  "P410-d1","P410-d7","P413-d1","P413-d7",
  "P466-d1","P466-d7","P447-d1","P447-d7",
  "P448-d1","P448-d7","P468-d1","P468-d7"
)

merged_int$orig.ident <- factor(merged_int$orig.ident, levels = sample_order)
celltypes_in_data <- unique(merged_int@meta.data$predicted.celltype.major)
colors_use <- my_colors[names(my_colors) %in% celltypes_in_data]

p_azimuth <- DimPlot(
  merged_int,
  reduction = "umap",
  group.by = "predicted.celltype.major",
  split.by = "orig.ident",
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  raster = FALSE,
  cols = colors_use
)

p_azimuth <- p_azimuth +
  ggplot2::facet_wrap(~ orig.ident, ncol = 6, scales = "fixed") +
  ggplot2::ggtitle("") +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    legend.position = "right"
  )

dir.create("./png", showWarnings = FALSE, recursive = TRUE)

png("./png/UMAP_split.png", 
    width = 12000, height = 5000, res = 600)
print(p_azimuth)
dev.off()


########### ref_umap_individual
p_azimuth <- DimPlot(
  merged_int,
  reduction = "ref.umap",
  group.by = "predicted.celltype.major",
  split.by = "orig.ident",
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  raster = FALSE,
  cols = colors_use
)

p_azimuth <- p_azimuth +
  ggplot2::facet_wrap(~ orig.ident, ncol = 6, scales = "fixed") +
  ggplot2::ggtitle("") +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    legend.position = "right"
  )

dir.create("./png", showWarnings = FALSE, recursive = TRUE)

png("./png/refUMAP_split.png", 
    width = 12000, height = 5000, res = 600)
print(p_azimuth)
dev.off()