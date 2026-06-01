############################################################################# QC #########################################
setwd("~/covid_scRNA/covid_R/NEW/Merge")

samples <- c("P410-d1","P410-d7","P413-d1","P413-d7",
             "P447-d1","P447-d7","P448-d1","P448-d7",
             "P466-d1","P466-d7","P468-d1","P468-d7")

prognosis_info <- c(
  "P410" = "Poor", "P413" = "Poor", "P466" = "Poor",
  "P447" = "Good", "P448" = "Good", "P468" = "Good"
)

base_dir <- "~/covid_scRNA/A_TBD230150_17225_20230317/02_cellranger_file"

seurat_list <- list()
for (sample in samples) {
  h5_path <- file.path(base_dir, sample, "outs", "raw_feature_bc_matrix.h5")
  seurat_data <- Read10X_h5(filename = h5_path)
  obj <- CreateSeuratObject(counts = seurat_data, project = sample, min.cells = 3)
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 30)
  
  group <- sub(".*-(d[17])$", "\\1", sample, ignore.case = TRUE)
  sample_name <- sub("-d[17]$", "", sample, ignore.case = TRUE)
  obj$Group <- group
  obj$Sample <- sample_name
  if (! sample_name %in% names(prognosis_info)) stop(paste("Sample", sample_name, "not found in prognosis_info"))
  obj$Prognosis <- unname(prognosis_info[sample_name])
  obj$Prognosis_Group <- paste0(obj$Prognosis, "_", obj$Group)
  
  seurat_list[[sample]] <- obj
  saveRDS(obj, file = paste0("~/covid_scRNA/covid_R/NEW/QC_before_merge/", gsub("-", "", sample), ".rds"))
}

merged_raw <- merge(seurat_list[[1]], y = seurat_list[-1], project = "COVID_merge")
saveRDS(merged_raw, "~/covid_scRNA/covid_R/NEW/merged_raw.rds")


############################################################ Integration & Clustering ################################################

setwd("~/covid_scRNA/covid_R/NEW/Merge")
obj_list <- SplitObject(merged_raw, split.by = "orig.ident")

obj_list <- lapply(obj_list, function(x){
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  x
})

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 2000)
anchors  <- FindIntegrationAnchors(
  object.list     = obj_list,
  anchor.features = features,
  dims            = 1:20,  
  k.anchor        = 5      
)
merged_int <- IntegrateData(anchorset = anchors)  
saveRDS(merged_int, "~/covid_scRNA/covid_R/NEW/merged_integration_all.rds")

DefaultAssay(merged_int) <- "integrated"
merged_int <- ScaleData(merged_int, verbose = FALSE)
merged_int <- RunPCA(merged_int, npcs = 30, verbose = FALSE)
merged_int <- RunUMAP(merged_int, dims = 1:15)
merged_int <- FindNeighbors(merged_int, dims = 1:15)
merged_int <- FindClusters(merged_int, resolution = 0.6)

saveRDS(merged_int, "~/covid_scRNA/covid_R/NEW/merged_integration_all.rds")

my_pub_colors <- c(
  "Good_d1"="#80B1D3","Good_d7"="#483D8B",
  "Poor_d1"="#FDB462","Poor_d7"="#D95F02"
)

my_base_cols <- c(
  "B"               = "#33A02C",
  "CD14 Mono"       = "#FD7F27",
  "CD16 Mono"       = "#FDBF6F",
  "CD4 T"           = "#1F78B4",
  "CD8 T"           = "#E31A1C",
  "NK"              = "#4ECCCB",
  "DC"              = "#B15928",
  "gd T cells"      = "#999999",
  "Granulocytes"    = "#FF7F00",
  "Progenitors"     = "#CAB2D6",
  "Intermediate Mono" = "#6A3D9A",
  "Other"           = "black"
)

## 2) Check seurat_clusters levels & calculate the number of clusters
merged_int$seurat_clusters <- factor(merged_int$seurat_clusters)
cluster_ids <- levels(merged_int$seurat_clusters)
n_clusters  <- length(cluster_ids)

## 3) Create a palette for the number of clusters by mixing base tones
cluster_cols_vec <- colorRampPalette(unname(my_base_cols))(n_clusters)

names(cluster_cols_vec) <- cluster_ids

p1 <- DimPlot(
  merged_int,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("") +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = cluster_cols_vec)

p1

dir.create("./png", showWarnings = FALSE)
ggsave("./png/UMAP_Integrated_Clusters.png", p1, width = 6, height = 5, dpi = 600)

png("./png/UMAP.png", width=4300, height=3500, res=600)

p1 <- DimPlot(
  merged_int,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "Prognosis_Group",
  label = TRUE,
  repel = TRUE, ncol=2
) +
  ggtitle("UMAP (Integrated) by Prognosis_Group") 

dir.create("./png", showWarnings = FALSE)
ggsave("./png/UMAP_Integrated_Clusters_Split.png", p1, width = 8, height = 8, dpi = 300)

cols <- c(
  "Good" = "#1F78B4",  
  "Poor" = "#E41A1C")

sample_order <- c(
  "P410-d1","P410-d7","P413-d1","P413-d7",
  "P466-d1","P466-d7","P447-d1","P447-d7",
  "P448-d1","P448-d7","P468-d1","P468-d7"
)

poor_samples <- c("P410-d1","P410-d7","P413-d1","P413-d7","P466-d1","P466-d7")

merged_int$orig.ident <- factor(merged_int$orig.ident, levels = sample_order)

total_counts <- merged_int@meta.data %>%
  as.data.frame() %>%
  mutate(orig.ident = factor(orig.ident, levels = sample_order),
         Prognosis = ifelse(orig.ident %in% poor_samples, "Poor", "Good")) %>%
  group_by(orig.ident, Prognosis) %>%
  summarise(total_cells = n(), .groups = "drop")

p_total_cells <- ggplot(total_counts, aes(x = orig.ident, y = total_cells, fill = Prognosis)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = total_cells), vjust = -0.3, size = 4) +
  scale_fill_manual(values = cols) +
  labs(x = "Sample", y = "Total Number of Cells", fill = "Prognosis") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, face = "bold"))
ggsave("./png/Seurat_clusters_total_cell_counts_integration.png",
       plot = p_total_cells, width = 12, height = 7, dpi = 300)