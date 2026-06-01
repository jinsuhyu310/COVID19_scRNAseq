# =============================================================================
# 06_T_Cell_Analysis.R
# CD8 T Cell Analysis:
#   UMAP visualization → Cluster highlighting → Monocle3 trajectory →
#   Marker dotplot → Cell proportion → DEG → Cytotoxicity/Exhaustion scores
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(patchwork)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(monocle3)
library(Matrix)

# -----------------------------------------------------------------------------
# 0. Shared settings
# -----------------------------------------------------------------------------
GROUP_LEVELS <- c("Poor_d1", "Poor_d7", "Good_d1", "Good_d7")
CLUSTER_COL  <- "seurat_clusters"

my_pub_colors <- c(
  "Good_d1" = "#80B1D3", "Good_d7" = "#483D8B",
  "Poor_d1" = "#FDB462", "Poor_d7" = "#D95F02"
)

# Marker gene sets
NAIVE_GENES <- c("CCR7","IL7R","TCF7","LEF1","SELL","LTB","MAL","TRAC")
EARLY_GENES <- c("FOS","JUN","JUNB","FOSB","NR4A1","NR4A2","NR4A3",
                 "IER2","EGR1","DUSP1","DUSP2")
EFF_GENES   <- c("NKG7","PRF1","GZMB","GNLY","CTSW","FGFBP2","IFNG","GZMH")
EXH_GENES   <- c("PDCD1","LAG3","TIGIT","HAVCR2","TOX","CXCL13","CTLA4","ENTPD1")
CYTO_GENES  <- c("NKG7","GNLY","PRF1","GZMB","GZMH","GZMK","IFNG")

# =============================================================================
# SECTION 1: Helper Functions
# =============================================================================

# Auto-detect UMAP name
get_umap_name <- function(obj) {
  if ("umap" %in% names(obj@reductions)) "umap"
  else if ("ref.umap" %in% names(obj@reductions)) "ref.umap"
  else stop("No UMAP reduction found (umap / ref.umap).")
}

# Safe layer join for Seurat v5
safe_join_layers <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
    suppressWarnings(tryCatch(JoinLayers(obj, assay = assay),
                              error = function(e) obj))
  } else obj
}

# Cluster color palette (auto-scaled)
make_cluster_colors <- function(cluster_levels) {
  n <- length(cluster_levels)
  cols <- if (n <= 12) {
    RColorBrewer::brewer.pal(max(n, 3), "Set3")[seq_len(n)]
  } else if (n <= 20) {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n)
  } else {
    grDevices::colorRampPalette(
      c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
        "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
        "#8DA0CB","#E78AC3")
    )(n)
  }
  setNames(cols, cluster_levels)
}

# Module score annotation per cluster (majority vote)
annotate_clusters_by_score <- function(obj, assay = "RNA",
                                       cluster_col = CLUSTER_COL) {
  obj <- safe_join_layers(obj, assay)
  DefaultAssay(obj) <- assay
  
  pick_present <- function(gs) gs[gs %in% rownames(obj)]
  ng <- pick_present(NAIVE_GENES)
  eg <- pick_present(EARLY_GENES)
  fg <- pick_present(EFF_GENES)
  xg <- pick_present(EXH_GENES)
  
  obj <- AddModuleScore(obj, list(ng), name = "S_Naive")
  obj <- AddModuleScore(obj, list(eg), name = "S_Early")
  obj <- AddModuleScore(obj, list(fg), name = "S_Eff")
  obj <- AddModuleScore(obj, list(xg), name = "S_Exh")
  
  md <- obj@meta.data %>% as.data.frame()
  cl <- md[[cluster_col]]
  if (is.data.frame(cl)) cl <- cl[, 1]
  md$cluster_chr <- trimws(as.character(cl))
  
  cl_mean <- md %>%
    dplyr::group_by(cluster_chr) %>%
    dplyr::summarise(
      n     = dplyr::n(),
      naive = mean(S_Naive1, na.rm = TRUE),
      early = mean(S_Early1, na.rm = TRUE),
      eff   = mean(S_Eff1,   na.rm = TRUE),
      exh   = mean(S_Exh1,   na.rm = TRUE),
      .groups = "drop"
    )
  
  thr <- list(
    naive = quantile(cl_mean$naive, 0.75, na.rm = TRUE),
    exh   = quantile(cl_mean$exh,   0.80, na.rm = TRUE),
    early = quantile(cl_mean$early, 0.85, na.rm = TRUE),
    eff   = quantile(cl_mean$eff,   0.65, na.rm = TRUE)
  )
  med_eff <- median(cl_mean$eff, na.rm = TRUE)
  
  cl_annot <- cl_mean %>%
    dplyr::mutate(
      state = dplyr::case_when(
        naive >= thr$naive & eff < med_eff  ~ "CD8 Naive-like",
        exh   >= thr$exh   & eff >= med_eff ~ "CD8 Exhausted",
        early >= thr$early                  ~ "CD8 Early Activated",
        eff   >= thr$eff                    ~ "CD8 Effector Memory",
        TRUE                                ~ "Unassigned"
      )
    )
  
  if (sum(cl_annot$state == "CD8 Exhausted") == 0) {
    top <- cl_annot$cluster_chr[which.max(cl_annot$exh)]
    cl_annot$state[cl_annot$cluster_chr == top] <- "CD8 Exhausted"
    message("[INFO] No Exhausted cluster -> forced: ", top)
  }
  
  map_state <- setNames(cl_annot$state, cl_annot$cluster_chr)
  state_vec  <- map_state[md$cluster_chr]
  names(state_vec) <- rownames(md)
  obj@meta.data[names(state_vec), "CD8_state4_conf"] <- unname(state_vec)
  
  message("[OK] CD8_state4_conf added.")
  print(table(obj$CD8_state4_conf, useNA = "ifany"))
  list(obj = obj, cl_annot = cl_annot)
}

# Monocle3 CDS builder (reuses Seurat UMAP)
build_cds_from_seurat <- function(obj, umap_name = NULL) {
  DefaultAssay(obj) <- "RNA"
  obj <- safe_join_layers(obj)
  if (is.null(umap_name)) umap_name <- get_umap_name(obj)
  
  counts    <- GetAssayData(obj, assay = "RNA", slot = "counts")
  cell_meta <- obj@meta.data
  gene_meta <- data.frame(gene_short_name = rownames(counts),
                          row.names = rownames(counts))
  
  cds <- new_cell_data_set(counts, cell_meta, gene_meta)
  um  <- Embeddings(obj, reduction = umap_name)[colnames(cds), 1:2]
  colnames(um) <- c("UMAP_1","UMAP_2")
  reducedDims(cds)$UMAP <- um
  
  colData(cds)$seurat_clusters <- factor(
    as.character(Idents(obj))[match(colnames(cds), colnames(obj))]
  )
  if ("CD8_state4_conf" %in% colnames(obj@meta.data)) {
    colData(cds)$CD8_state4_conf <- factor(
      obj@meta.data[colnames(cds), "CD8_state4_conf"]
    )
  }
  
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = TRUE)
  cds
}

# Root cells: core-naive high (TCF7/IL7R/CCR7 top p%)
select_root_cells <- function(cds, root_label = "CD8 Naive-like", p_root = 0.95) {
  core_genes <- intersect(c("TCF7","IL7R","CCR7"), rownames(cds))
  stopifnot(length(core_genes) >= 2)
  
  naive_cells <- rownames(colData(cds))[
    !is.na(colData(cds)$CD8_state4_conf) &
      colData(cds)$CD8_state4_conf == root_label
  ]
  stopifnot(length(naive_cells) >= 50)
  
  M_raw <- SingleCellExperiment::counts(cds)[core_genes, naive_cells, drop = FALSE]
  ok    <- Matrix::colSums(M_raw > 0) >= 2
  cand  <- naive_cells[ok]
  stopifnot(length(cand) > 30)
  
  core_score <- Matrix::colSums(log1p(M_raw[, cand, drop = FALSE]))
  thr        <- quantile(core_score, probs = p_root, na.rm = TRUE)
  root_cells <- names(core_score)[core_score >= thr]
  
  if (length(root_cells) < 30) {
    root_cells <- names(sort(core_score, decreasing = TRUE))[seq_len(min(100, length(core_score)))]
    message("[WARN] root fallback: top N=", length(root_cells))
  }
  message("[INFO] root_cells=", length(root_cells))
  
  cells_all  <- rownames(colData(cds))
  core_vec   <- setNames(rep(NA_real_, length(cells_all)), cells_all)
  core_vec[names(core_score)] <- core_score
  colData(cds)$core_naive_score <- core_vec
  list(cds = cds, root_cells = root_cells)
}

# Cytotoxicity/Exhaustion module score helper
add_module_scores_safe <- function(obj, assay = "RNA", prefix = "Tcell") {
  obj <- safe_join_layers(obj, assay)
  DefaultAssay(obj) <- assay
  present   <- rownames(obj[[assay]])
  cyto_use  <- intersect(CYTO_GENES,  present)
  exh_use   <- intersect(EXH_GENES,   present)
  stopifnot(length(cyto_use) >= 3, length(exh_use) >= 3)
  
  obj <- AddModuleScore(obj, list(cyto_use),
                        name = paste0(prefix, "_CYTO_"), assay = assay)
  obj <- AddModuleScore(obj, list(exh_use),
                        name = paste0(prefix, "_EXH_"),  assay = assay)
  list(obj      = obj,
       cyto_col = paste0(prefix, "_CYTO_1"),
       exh_col  = paste0(prefix, "_EXH_1"),
       cyto_use = cyto_use,
       exh_use  = exh_use)
}

make_dot_df <- function(obj, group_col, subtype_col, score_col) {
  obj@meta.data %>% as.data.frame() %>%
    dplyr::transmute(
      group   = as.character(.data[[group_col]]),
      subtype = as.character(.data[[subtype_col]]),
      score   = as.numeric(.data[[score_col]])
    ) %>%
    dplyr::filter(!is.na(group), !is.na(subtype), is.finite(score)) %>%
    dplyr::group_by(subtype, group) %>%
    dplyr::summarise(avg     = mean(score, na.rm = TRUE),
                     pct_pos = mean(score > 0, na.rm = TRUE),
                     .groups = "drop")
}

plot_dot_score <- function(df_sum, title = NULL) {
  ggplot(df_sum, aes(x = group, y = subtype)) +
    geom_point(aes(size = pct_pos, color = avg), alpha = 0.95) +
    scale_size_continuous(name   = "% score > 0",
                          labels = scales::percent_format(accuracy = 1),
                          limits = c(0, 1), range = c(1.5, 9)) +
    scale_color_viridis_c(name = "Mean module score") +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          plot.title   = element_text(face = "bold", hjust = 0.5))
}

plot_box_score <- function(df, title = NULL, facet_col = NULL) {
  p <- ggplot(df, aes(x = group, y = score)) +
    geom_boxplot(outlier.size = 0.2) +
    labs(title = title, x = NULL, y = "Module score") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold", hjust = 0.5))
  if (!is.null(facet_col)) p <- p + facet_wrap(as.formula(paste("~", facet_col)),
                                               scales = "free_y")
  p
}

save_fig <- function(p, filename, width, height, dpi = 300) {
  tryCatch({
    ggsave(filename, p, width = width, height = height, dpi = dpi,
           limitsize = FALSE, bg = "white")
    message("[Saved] ", filename)
  }, error = function(e) {
    fn2 <- sub("\\.png$", "_fallback.png", filename)
    ggsave(fn2, p, width = width * 0.9, height = height * 0.9,
           dpi = dpi, limitsize = FALSE, bg = "white")
    message("[Saved fallback] ", fn2)
  })
}

# =============================================================================
# SECTION 2: Load CD8 T Object & Preprocessing
# =============================================================================
message("\n>>> Preprocessing CD8 T cell object")

obj <- T_CD8
stopifnot(inherits(obj, "Seurat"))

DefaultAssay(obj) <- "RNA"
obj <- safe_join_layers(obj)
Idents(obj) <- CLUSTER_COL

obj$Prognosis_Group <- factor(
  trimws(as.character(obj$Prognosis_Group)),
  levels = GROUP_LEVELS
)

# Numeric cluster ordering
clu_chr    <- as.character(obj$seurat_clusters)
clu_num    <- suppressWarnings(as.integer(clu_chr))
clu_levels <- if (any(is.na(clu_num))) sort(unique(clu_chr)) else
  as.character(sort(unique(clu_num)))
obj$seurat_clusters <- factor(clu_chr, levels = clu_levels)
Idents(obj) <- CLUSTER_COL

cluster_levels <- levels(obj$seurat_clusters)
cluster_colors <- make_cluster_colors(cluster_levels)

umap_name <- get_umap_name(obj)
wd        <- getwd()
message("[INFO] Working directory: ", wd)

# =============================================================================
# SECTION 3: CD8 T UMAP Plots
# =============================================================================
message("\n>>> CD8 T UMAP plots")

# 3-A: Full 4-group split UMAP
p_split4 <- DimPlot(
  obj,
  reduction  = umap_name,
  group.by   = CLUSTER_COL,
  split.by   = "Prognosis_Group",
  label      = TRUE, label.size = 4,
  repel      = TRUE, raster = FALSE,
  cols       = cluster_colors, ncol = 2
) +
  ggtitle("CD8 T") +
  theme(plot.title        = element_text(hjust = 0.5, face = "bold"),
        strip.background  = element_blank(),
        strip.text        = element_text(size = 14, face = "bold"),
        legend.title      = element_blank())

save_fig(p_split4, "TCD8_UMAP_Prognosis_Group.png", 10, 9)

# 3-B: Pair UMAPs (Poor d1/d7, Good d1/d7, d1 cross, d7 cross)
save_pair_umap <- function(sub_obj, levels_use, filename, title_txt) {
  sub_obj$Prognosis_Group <- factor(
    as.character(sub_obj$Prognosis_Group), levels = levels_use
  )
  p <- DimPlot(sub_obj, reduction = umap_name, group.by = CLUSTER_COL,
               split.by = "Prognosis_Group", label = TRUE, label.size = 4,
               repel = TRUE, raster = FALSE, cols = cluster_colors, ncol = 2) +
    ggtitle(title_txt) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_fig(p, filename, 10, 4.5)
}

save_pair_umap(subset(obj, Prognosis_Group %in% c("Poor_d1","Poor_d7")),
               c("Poor_d1","Poor_d7"), "CD8_T_UMAP_Poor_d1_vs_d7.png",
               "CD8 T – Poor d1 vs d7")

save_pair_umap(subset(obj, Prognosis_Group %in% c("Good_d1","Good_d7")),
               c("Good_d1","Good_d7"), "CD8_T_UMAP_Good_d1_vs_d7.png",
               "CD8 T – Good d1 vs d7")

save_pair_umap(subset(obj, Prognosis_Group %in% c("Poor_d1","Good_d1")),
               c("Poor_d1","Good_d1"), "CD8_T_UMAP_Poor_d1_vs_Good_d1.png",
               "CD8 T – Poor vs Good (Day 1)")

save_pair_umap(subset(obj, Prognosis_Group %in% c("Poor_d7","Good_d7")),
               c("Poor_d7","Good_d7"), "CD8_T_UMAP_Poor_d7_vs_Good_d7.png",
               "CD8 T – Poor vs Good (Day 7)")

# =============================================================================
# SECTION 4: Cluster Highlight UMAP (Clusters 1 & 3)
# =============================================================================
message("\n>>> Highlight clusters 1 and 3")

FOCUS_CLUSTERS <- c("1","3")

obj$focus_1_3 <- factor(
  dplyr::case_when(
    as.character(obj[[CLUSTER_COL]][, 1]) == "1" ~ "Cluster 1",
    as.character(obj[[CLUSTER_COL]][, 1]) == "3" ~ "Cluster 3",
    TRUE ~ "Other clusters"
  ),
  levels = c("Cluster 1","Cluster 3","Other clusters")
)

pal_focus <- c("Cluster 1" = "#E64B35", "Cluster 3" = "#00A087",
               "Other clusters" = "#E0E0E0")

pC <- DimPlot(obj, reduction = umap_name, group.by = "focus_1_3",
              label = FALSE, pt.size = 0.35, cols = pal_focus, order = TRUE) +
  theme_classic(base_size = 10) +
  theme(plot.title   = element_blank(),
        axis.text    = element_text(size = 13.5),
        axis.title   = element_text(size = 15),
        axis.ticks   = element_blank(),
        axis.line    = element_line(linewidth = 0.6),
        legend.title = element_blank(),
        legend.text  = element_text(size = 12)) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

# Manual label positions
emb      <- Embeddings(obj, reduction = umap_name)
df_lab   <- data.frame(UMAP_1 = emb[,1], UMAP_2 = emb[,2],
                       label_g = obj$focus_1_3)
centers  <- df_lab %>%
  dplyr::filter(label_g != "Other clusters") %>%
  dplyr::group_by(label_g) %>%
  dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),
                   .groups = "drop") %>%
  dplyr::mutate(
    UMAP_1 = dplyr::case_when(label_g == "Cluster 1" ~ UMAP_1 - 3.2,
                              label_g == "Cluster 3" ~ UMAP_1 - 3.0,
                              TRUE ~ UMAP_1),
    UMAP_2 = dplyr::case_when(label_g == "Cluster 1" ~ UMAP_2 + 0.2,
                              label_g == "Cluster 3" ~ UMAP_2 - 0.1,
                              TRUE ~ UMAP_2)
  )

pC <- pC +
  geom_text(data = centers, aes(x = UMAP_1, y = UMAP_2, label = label_g),
            fontface = "bold", size = 4, color = "black")

save_fig(pC, "TCD8_UMAP_Highlight_Clusters_1_3.png", 5.5, 4)

# =============================================================================
# SECTION 5: CD8 T Cell State Annotation & Marker Dotplot
# =============================================================================
message("\n>>> CD8 state annotation by module score")

if (!"CD8_state4_conf" %in% colnames(obj@meta.data)) {
  res_annot <- annotate_clusters_by_score(obj)
  obj       <- res_annot$obj
  cl_annot  <- res_annot$cl_annot
  print(cl_annot %>% dplyr::arrange(state, dplyr::desc(n)))
} else {
  message("[INFO] CD8_state4_conf already exists.")
}

# Marker dotplot
pick_present <- function(gs) gs[gs %in% rownames(obj)]
gene_sets <- list(
  "Naive/Memory"       = pick_present(NAIVE_GENES),
  "Early activation"   = pick_present(EARLY_GENES),
  "Effector/Cytotoxic" = pick_present(EFF_GENES),
  "Exhaustion"         = pick_present(EXH_GENES)
)

p_dot <- DotPlot(obj, features = gene_sets, assay = "RNA") +
  RotatedAxis() +
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 9),
        strip.text  = element_text(size = 11, face = "bold"),
        plot.title  = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "CD8 T: Markers", x = NULL, y = "Cluster")

save_fig(p_dot, "CD8_T_DotPlot_Markers.png", 11, 5.5)

# =============================================================================
# SECTION 6: CD8 T Subtypes UMAP (Azimuth predicted.celltype.major_fine)
# =============================================================================
message("\n>>> CD8 T subtypes UMAP")

stopifnot("predicted.celltype.major_fine" %in% colnames(obj@meta.data))

cd8_subtype_levels <- c("CD8 Naive","CD8 TCM","CD8 TEM","CD8 Proliferating")
obj$predicted.celltype.major_fine <- factor(
  as.character(obj$predicted.celltype.major_fine),
  levels = cd8_subtype_levels
)

cd8_colors <- c("CD8 Naive" = "#1F78B4", "CD8 TCM" = "#33A02C",
                "CD8 TEM" = "#E31A1C", "CD8 Proliferating" = "#FF7F00")

p_subtypes <- DimPlot(obj, reduction = umap_name,
                      group.by = "predicted.celltype.major_fine",
                      cols = cd8_colors, label = FALSE, pt.size = 0.25) +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10),
        plot.title   = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "CD8 T Subtypes")

# Add centroid labels
emb_df  <- as.data.frame(Embeddings(obj, umap_name))
colnames(emb_df)[1:2] <- c("UMAP1","UMAP2")
df_all  <- cbind(obj@meta.data, emb_df)

centroids <- df_all %>%
  dplyr::filter(!is.na(predicted.celltype.major_fine)) %>%
  dplyr::group_by(predicted.celltype.major_fine) %>%
  dplyr::summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

p_subtypes <- p_subtypes +
  ggrepel::geom_label_repel(
    data = centroids,
    aes(x = UMAP1, y = UMAP2, label = predicted.celltype.major_fine),
    inherit.aes = FALSE, size = 5.2, fontface = "bold",
    color = "black", fill = "white", label.size = 0.4,
    box.padding = 0.6, point.padding = 0.5,
    segment.color = "grey40", max.overlaps = Inf
  )

save_fig(p_subtypes, "TCD8_UMAP_4subtypes.png", 7.8, 6.8)

# =============================================================================
# SECTION 7: Monocle3 Trajectory (Core-Naive Root)
# =============================================================================
message("\n>>> Monocle3 trajectory analysis")

Idents(obj) <- CLUSTER_COL
cds <- build_cds_from_seurat(obj, umap_name)

res_root  <- select_root_cells(cds, root_label = "CD8 Naive-like", p_root = 0.95)
cds       <- res_root$cds
cds       <- order_cells(cds, root_cells = res_root$root_cells)

pt     <- pseudotime(cds)
part   <- partitions(cds)
main_p <- names(sort(table(part), decreasing = TRUE))[1]
message("[INFO] main partition: ", main_p)

pt_main <- pt
pt_main[part != main_p] <- NA
pt_main[is.infinite(pt_main)] <- NA

colData(cds)$pseudotime_main     <- pt_main
colData(cds)$pseudotime_main_log <- log1p(pt_main)

# Transfer pseudotime to Seurat
obj$monocle3_pseudotime_main     <- colData(cds)$pseudotime_main[colnames(obj)]
obj$monocle3_pseudotime_main_log <- colData(cds)$pseudotime_main_log[colnames(obj)]

# Trajectory plots
p_traj <- plot_cells(
  cds, reduction_method = "UMAP",
  color_cells_by = "pseudotime_main",
  label_cell_groups = FALSE, label_leaves = FALSE,
  label_branch_points = FALSE, label_roots = FALSE,
  graph_label_size = 0, show_trajectory_graph = TRUE
) +
  ggtitle("CD8T lineage trajectory") +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 13),
    legend.text     = element_text(size = 12),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title      = element_text(face = "bold", size = 14),
    axis.text       = element_text(face = "bold", size = 10.5),
    axis.line       = element_line(linewidth = 1.6)
  )

save_fig(p_traj, "TCD8_Monocle3_pseudotime.png", 5, 6)

p_traj_state <- plot_cells(
  cds, reduction_method = "UMAP",
  color_cells_by = "CD8_state4_conf",
  label_cell_groups = FALSE, label_leaves = FALSE,
  label_branch_points = FALSE, label_roots = FALSE,
  graph_label_size = 0, show_trajectory_graph = TRUE
) + ggtitle("Trajectory: CD8 state")

save_fig(p_traj_state, "TCD8_Monocle3_state4.png", 7.2, 6.2)

p_traj_clust <- plot_cells(
  cds, reduction_method = "UMAP",
  color_cells_by = "seurat_clusters",
  label_cell_groups = TRUE, label_leaves = FALSE,
  label_branch_points = FALSE, label_roots = FALSE,
  graph_label_size = 0, show_trajectory_graph = TRUE
) + ggtitle("Trajectory: Seurat clusters")

save_fig(p_traj_clust, "TCD8_Monocle3_clusters.png", 7.2, 6.2)

# =============================================================================
# SECTION 8: Cell Proportion – Clusters 1 & 3 (Paired Boxplot)
# =============================================================================
message("\n>>> Cluster 1 paired proportion boxplot")

meta_TCD8 <- obj@meta.data %>%
  as.data.frame() %>%
  dplyr::mutate(
    orig.ident      = as.character(orig.ident),
    seurat_clusters = as.character(unlist(seurat_clusters)),
    Group           = as.character(unlist(Group)),
    Prognosis       = as.character(unlist(Prognosis)),
    PatientID       = stringr::str_remove(orig.ident, "-d[17]$")
  )

sample_meta <- meta_TCD8 %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(Prognosis = Prognosis[1], Group = Group[1],
                   PatientID = PatientID[1], total_n = dplyr::n(),
                   .groups = "drop")

cluster_counts <- meta_TCD8 %>%
  dplyr::filter(seurat_clusters %in% c("1")) %>%
  dplyr::group_by(orig.ident, seurat_clusters) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

cluster_prop_c1 <- sample_meta %>%
  dplyr::select(orig.ident, Prognosis, Group, PatientID, total_n) %>%
  tidyr::crossing(seurat_clusters = "1") %>%
  dplyr::left_join(cluster_counts, by = c("orig.ident","seurat_clusters")) %>%
  dplyr::mutate(
    n_cells         = ifelse(is.na(n_cells), 0L, n_cells),
    percent         = (n_cells / total_n) * 100,
    Prognosis_Group = factor(paste0(Prognosis, "_", Group), levels = GROUP_LEVELS)
  )

paired_ids_c1 <- cluster_prop_c1 %>%
  dplyr::group_by(Prognosis, PatientID) %>%
  dplyr::summarise(has_d1 = "d1" %in% Group, has_d7 = "d7" %in% Group,
                   .groups = "drop") %>%
  dplyr::filter(has_d1 & has_d7) %>%
  dplyr::pull(PatientID)

paired_c1 <- cluster_prop_c1 %>% dplyr::filter(PatientID %in% paired_ids_c1)

run_paired_ttest <- function(df) {
  d1 <- df %>% dplyr::filter(Group == "d1") %>% dplyr::arrange(PatientID)
  d7 <- df %>% dplyr::filter(Group == "d7") %>% dplyr::arrange(PatientID)
  if (nrow(d1) > 1 && nrow(d1) == nrow(d7))
    t.test(d1$percent, d7$percent, paired = TRUE)$p.value
  else NA
}

ttest_c1 <- paired_c1 %>%
  dplyr::group_by(seurat_clusters, Prognosis) %>%
  dplyr::summarise(p_value = run_paired_ttest(dplyr::cur_data()), .groups = "drop") %>%
  dplyr::mutate(
    sig_label = dplyr::case_when(
      is.na(p_value)   ~ "ns", p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**", p_value <= 0.05  ~ "*", TRUE ~ "ns"
    )
  )

p_prop_c1 <- ggplot(paired_c1,
                    aes(Prognosis_Group, percent, fill = Prognosis_Group)) +
  geom_boxplot(width = 0.7, alpha = 0.85, outlier.shape = NA, color = "gray30") +
  geom_line(aes(group = interaction(PatientID, Prognosis)),
            color = "gray60", alpha = 0.5, linewidth = 0.5) +
  geom_point(aes(group = interaction(PatientID, Prognosis)),
             size = 1.5, alpha = 0.85, color = "black") +
  facet_wrap(~ seurat_clusters, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = my_pub_colors, drop = FALSE) +
  labs(title = "CD8 T: Cluster 1", x = NULL, y = "Cell proportion (%)") +
  theme_classic(base_size = 11) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "gray95", color = NA),
        strip.text       = element_text(size = 9),
        plot.title       = element_text(hjust = 0.5, size = 12, face = "bold"))

pvalue_labels_c1 <- ttest_c1 %>%
  dplyr::left_join(
    paired_c1 %>% dplyr::group_by(seurat_clusters) %>%
      dplyr::summarise(y_pos = max(percent, na.rm = TRUE) * 1.07, .groups = "drop"),
    by = "seurat_clusters"
  )

p_prop_c1 <- p_prop_c1 +
  geom_text(
    data = pvalue_labels_c1,
    aes(x = ifelse(Prognosis == "Poor", 1.5, 3.5),
        y = y_pos, label = paste(Prognosis, sig_label)),
    size = 3.6, color = "black", inherit.aes = FALSE
  )

png("TCD8_Clusters_1_Paired_BoxPlot.png", width = 1200, height = 900, res = 300)
print(p_prop_c1); dev.off()
message("[Saved] TCD8_Clusters_1_Paired_BoxPlot.png")

# =============================================================================
# SECTION 9: Cell Proportion – Poor vs Good (All Clusters, Unpaired)
# =============================================================================
message("\n>>> Unpaired cell proportion: Poor vs Good by cluster")

T_use <- T_CD8
T_use <- safe_join_layers(T_use)
meta  <- T_use@meta.data %>% as.data.frame() %>%
  dplyr::mutate(
    orig.ident      = as.character(orig.ident),
    seurat_clusters = as.character(unlist(seurat_clusters)),
    Prognosis_Group = if ("Prognosis_Group" %in% colnames(.)) {
      as.character(unlist(Prognosis_Group))
    } else {
      paste0(as.character(unlist(Prognosis)), "_", as.character(unlist(Group)))
    }
  )

sample_meta_all <- meta %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(Prognosis_Group = Prognosis_Group[1],
                   total_n = dplyr::n(), .groups = "drop")

all_cl <- as.character(sort(unique(suppressWarnings(as.integer(meta$seurat_clusters)))))
if (any(is.na(suppressWarnings(as.integer(meta$seurat_clusters)))))
  all_cl <- sort(unique(meta$seurat_clusters))

df_full_all <- sample_meta_all %>%
  tidyr::crossing(seurat_clusters = all_cl) %>%
  dplyr::left_join(
    meta %>% dplyr::group_by(orig.ident, seurat_clusters) %>%
      dplyr::summarise(n_cells = dplyr::n(), .groups = "drop"),
    by = c("orig.ident","seurat_clusters")
  ) %>%
  dplyr::mutate(n_cells = ifelse(is.na(n_cells), 0L, n_cells),
                percent = (n_cells / total_n) * 100,
                seurat_clusters = factor(seurat_clusters, levels = all_cl))

make_unpaired_cluster_plot <- function(df, g1, g2, title_txt, out_png,
                                       nrow_facet = 2,
                                       width = 14, height = 6) {
  df_sub <- df %>%
    dplyr::filter(Prognosis_Group %in% c(g1, g2)) %>%
    dplyr::mutate(Prognosis_Group = factor(Prognosis_Group, levels = c(g1, g2)))
  
  pvals <- df_sub %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(
      n_g1    = sum(Prognosis_Group == g1),
      n_g2    = sum(Prognosis_Group == g2),
      p.value = if (sum(Prognosis_Group == g1) >= 2 & sum(Prognosis_Group == g2) >= 2) {
        t.test(percent[Prognosis_Group == g1],
               percent[Prognosis_Group == g2], paired = FALSE)$p.value
      } else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      signif = dplyr::case_when(
        is.na(p.value)  ~ "NA",  p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",  p.value < 0.05  ~ "*", TRUE ~ "ns"
      )
    ) %>%
    dplyr::left_join(
      df_sub %>% dplyr::group_by(seurat_clusters) %>%
        dplyr::summarise(y.position = {
          mx <- suppressWarnings(max(percent, na.rm = TRUE))
          if (is.finite(mx)) mx * 1.15 + 0.2 else 0.5
        }, .groups = "drop"),
      by = "seurat_clusters"
    )
  
  p <- ggplot(df_sub, aes(Prognosis_Group, percent, fill = Prognosis_Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.7, color = "gray30") +
    geom_point(shape = 21, color = "black", size = 2.1, alpha = 0.85) +
    facet_wrap(~ seurat_clusters, scales = "free_y", nrow = nrow_facet) +
    scale_fill_manual(values = my_pub_colors, breaks = c(g1, g2), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    geom_text(
      data = pvals,
      aes(x = 1.5, y = y.position,
          label = ifelse(is.na(p.value), "n<2",
                         paste0("p=", signif(p.value, 3), " (", signif, ")"))),
      inherit.aes = FALSE, size = 5.2
    ) +
    labs(x = NULL, y = "Cell proportion (%)", title = title_txt) +
    theme_classic(base_size = 16) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 15,
                                          color = "black"),
          axis.text.y      = element_text(size = 12, color = "black"),
          strip.text       = element_text(face = "bold", size = 14),
          strip.background = element_rect(color = "black", fill = "grey90",
                                          linewidth = 0.8),
          legend.position  = "none",
          plot.title       = element_text(hjust = 0.5, face = "bold"))
  
  save_fig(p, out_png, width, height)
  invisible(p)
}

make_unpaired_cluster_plot(
  df_full_all, "Poor_d1","Good_d1",
  "CD8 T: Poor vs Good (Day 1)",
  "CD8T_Poor_vs_Good_cluster_proportion_pre.png",
  nrow_facet = 2, width = 14, height = 6
)

make_unpaired_cluster_plot(
  df_full_all, "Poor_d7","Good_d7",
  "CD8 T: Poor vs Good (Day 7)",
  "CD8T_Poor_vs_Good_cluster_proportion_post.png",
  nrow_facet = 2, width = 14, height = 6
)

# =============================================================================
# SECTION 10: CD8 TCM Proportion (Paired Boxplot, Poor & Good)
# =============================================================================
message("\n>>> CD8 TCM proportion paired boxplot")

obj_tcm <- T_CD8
obj_tcm <- safe_join_layers(obj_tcm)
obj_tcm$Prognosis_Group <- factor(as.character(obj_tcm$Prognosis_Group),
                                  levels = GROUP_LEVELS)

md_tcm <- obj_tcm@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::transmute(
    cell, orig.ident = as.character(orig.ident),
    Prognosis_Group  = as.character(Prognosis_Group),
    fineT = trimws(as.character(predicted.celltype.major_fine))
  ) %>%
  dplyr::filter(!is.na(orig.ident), Prognosis_Group %in% GROUP_LEVELS) %>%
  dplyr::mutate(
    Prognosis_Group = factor(Prognosis_Group, levels = GROUP_LEVELS),
    Patient_ID      = sub("-d[17]$", "", orig.ident)
  )

denom_tcm <- md_tcm %>%
  dplyr::group_by(orig.ident, Prognosis_Group, Patient_ID) %>%
  dplyr::summarise(total_cells = dplyr::n(), .groups = "drop")

num_tcm <- md_tcm %>%
  dplyr::filter(fineT == "CD8 TCM") %>%
  dplyr::group_by(orig.ident, Prognosis_Group, Patient_ID, fineT) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

prop_tcm <- denom_tcm %>%
  tidyr::crossing(fineT = "CD8 TCM") %>%
  dplyr::left_join(num_tcm, by = c("orig.ident","Prognosis_Group","Patient_ID","fineT")) %>%
  dplyr::mutate(n_cells = tidyr::replace_na(n_cells, 0L),
                percent = ifelse(total_cells > 0, n_cells / total_cells * 100, 0),
                fineT   = factor(fineT, levels = "CD8 TCM"))

make_paired_plot <- function(prop_df, prognosis, out_png) {
  groups_use <- if (prognosis == "Poor") c("Poor_d1","Poor_d7") else
    c("Good_d1","Good_d7")
  g1 <- groups_use[1]; g2 <- groups_use[2]
  
  df_sub <- prop_df %>%
    dplyr::filter(Prognosis_Group %in% groups_use) %>%
    dplyr::mutate(Prognosis_Group = factor(as.character(Prognosis_Group),
                                           levels = groups_use))
  
  wide <- df_sub %>%
    dplyr::select(fineT, Patient_ID, Prognosis_Group, percent) %>%
    tidyr::pivot_wider(names_from = Prognosis_Group, values_from = percent)
  
  pvals <- wide %>%
    dplyr::group_by(fineT) %>%
    dplyr::summarise(
      n_pairs = sum(!is.na(.data[[g1]]) & !is.na(.data[[g2]])),
      p.value = if (sum(!is.na(.data[[g1]]) & !is.na(.data[[g2]])) >= 2) {
        t.test(.data[[g2]], .data[[g1]], paired = TRUE)$p.value
      } else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      signif = dplyr::case_when(
        is.na(p.value)  ~ "NA",  p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",  p.value < 0.05  ~ "*", TRUE ~ "ns"
      )
    ) %>%
    dplyr::left_join(
      df_sub %>% dplyr::group_by(fineT) %>%
        dplyr::summarise(y.position = {
          mx <- suppressWarnings(max(percent, na.rm = TRUE))
          if (is.finite(mx)) mx * 1.10 else 1
        }, .groups = "drop"),
      by = "fineT"
    )
  
  p <- ggplot(df_sub,
              aes(Prognosis_Group, percent, fill = Prognosis_Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_line(aes(group = Patient_ID), color = "gray40",
              linewidth = 0.6, alpha = 0.7) +
    geom_point(shape = 21, color = "black", size = 2) +
    facet_wrap(~ fineT, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = my_pub_colors, breaks = groups_use) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    geom_text(
      data = pvals,
      aes(x = 1.5, y = y.position,
          label = ifelse(is.na(p.value), "n<2",
                         paste0("p=", signif(p.value, 3), " (", signif, ")"))),
      inherit.aes = FALSE, size = 5
    ) +
    labs(x = NULL, y = "Percentage of cells (%)", title = "") +
    theme_classic(base_size = 15) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 12,
                                          color = "black"),
          strip.text       = element_text(face = "bold", size = 14),
          strip.background = element_rect(color = "black", fill = "grey90",
                                          linewidth = 1),
          legend.position  = "none")
  
  save_fig(p, out_png, 3, 4)
  invisible(p)
}

make_paired_plot(prop_tcm, "Poor", "T_CD8_TCM_proportion_Poor.png")
make_paired_plot(prop_tcm, "Good", "T_CD8_TCM_proportion_Good.png")

# =============================================================================
# SECTION 11: Cluster 1 vs Cluster 3 DEG Quadrant Plot
# =============================================================================
message("\n>>> Cluster 1 vs Cluster 3 DEG quadrant plot")

Idents(obj) <- CLUSTER_COL

message("Calculating Cluster 1 DEGs...")
deg1 <- FindMarkers(obj, ident.1 = "1", only.pos = FALSE,
                    logfc.threshold = 0, min.pct = 0.1) %>%
  as.data.frame() %>% tibble::rownames_to_column("gene")

message("Calculating Cluster 3 DEGs...")
deg3 <- FindMarkers(obj, ident.1 = "3", only.pos = FALSE,
                    logfc.threshold = 0, min.pct = 0.1) %>%
  as.data.frame() %>% tibble::rownames_to_column("gene")

# Remove housekeeping genes
rm_genes <- function(df) {
  df %>% dplyr::filter(!grepl("^(MT-|RPS|RPL)", gene),
                       !gene %in% c("MALAT1","NEAT1","MTRNR2L12"))
}
deg1 <- rm_genes(deg1)
deg3 <- rm_genes(deg3)

m <- dplyr::full_join(
  deg1 %>% dplyr::select(gene, fc1 = avg_log2FC),
  deg3 %>% dplyr::select(gene, fc3 = avg_log2FC),
  by = "gene"
) %>%
  dplyr::mutate(fc1 = tidyr::replace_na(fc1, 0),
                fc3 = tidyr::replace_na(fc3, 0))

FC_LIMIT <- 0.5
m_outside <- m %>% dplyr::filter(!(abs(fc1) < FC_LIMIT & abs(fc3) < FC_LIMIT))

# Genes of interest
q1_genes <- c("CARD11","ABCA2","NBEAL2","MEF2D","RNF213","ARID1B")
q3_genes <- c("CCR7","IL7R","TCF7","LEF1","SELL","LTB","MAL","TRAC")

lab_q1 <- m_outside %>%
  dplyr::filter(gene %in% q1_genes, fc1 >= FC_LIMIT, fc3 >= FC_LIMIT) %>%
  dplyr::mutate(group = "Shared backbone (Q1)")

lab_q3 <- m_outside %>%
  dplyr::filter(gene %in% q3_genes, fc1 <= -FC_LIMIT, fc3 <= -FC_LIMIT) %>%
  dplyr::mutate(group = "Naive-like (Q3)")

lab_final <- dplyr::bind_rows(lab_q1, lab_q3)
cat("[INFO] Labeled genes:\n"); print(lab_final$gene)

DEG_COLORS <- c("Shared backbone (Q1)" = "darkred", "Naive-like (Q3)" = "darkblue")

p_quad <- ggplot(m_outside, aes(x = fc1, y = fc3)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  geom_vline(xintercept = c(-FC_LIMIT, FC_LIMIT),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-FC_LIMIT, FC_LIMIT),
             linetype = "dashed", color = "grey50") +
  geom_point(color = "grey75", size = 0.6, alpha = 0.6) +
  geom_point(data = lab_final, aes(color = group), size = 1) +
  ggrepel::geom_text_repel(data = lab_final,
                           aes(label = gene, color = group),
                           size = 2.7, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = DEG_COLORS) +
  labs(title = "CD8 T: Cluster 1 vs Cluster 3 DEG",
       x = "Log2 FC (Cluster 1 vs Others)",
       y = "Log2 FC (Cluster 3 vs Others)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank(), legend.position = "bottom")

save_fig(p_quad, "CD8T_cluster1_3_DEG_quadrant.png", 3.5, 5.5)

# =============================================================================
# SECTION 12: CD8 TEM Heatmap
# =============================================================================
message("\n>>> CD8 TEM gene expression heatmap")

obj_tem <- subset(T_CD8, subset = predicted.celltype.major_fine == "CD8 TEM")
obj_tem <- safe_join_layers(obj_tem)

obj_tem$FourGroup <- factor(
  paste0(obj_tem$Prognosis, "_", obj_tem$Group),
  levels = c("Poor_d1","Poor_d7","Good_d1","Good_d7")
)

# Ordered gene groups (narrative pattern)
gene_groups_tem <- list(
  Poor_Up_d7   = c("JUN","JUNB","IRF1","HSPA5","LDHA","GAPDH",
                   "EEF1A1","EEF1B2","TPT1","RBM3"),
  Poor_High_d1 = c("IFI44L","IFI27","ISG15","MX1","OAS1","OAS2",
                   "IRF7","ZBP1","XAF1","SLFN5"),
  Good_High_d1 = c("CD2","SLAMF6","PRF1","KLRD1"),
  Good_Up_d7   = c("GZMB","GZMA","GZMK","NKG7","CD3D","CD3E","CD27",
                   "LIME1","IL32","KLRC2","KIR3DL2","KIR2DL3","CX3CR1",
                   "TGFB1","B2M","HLA-A","HLA-C","IFNGR1")
)

heat_genes_tem <- intersect(unique(unlist(gene_groups_tem)), rownames(obj_tem))

avg_tem <- AggregateExpression(obj_tem, group.by = "FourGroup",
                               assays = "RNA", features = heat_genes_tem,
                               return.seurat = FALSE, slot = "data")$RNA

col_order <- c("Poor-d1","Poor-d7","Good-d1","Good-d7")
col_order <- col_order[col_order %in% colnames(avg_tem)]
avg_tem   <- avg_tem[, col_order, drop = FALSE]

mat_z   <- t(scale(t(as.matrix(avg_tem))))
mat_z[is.na(mat_z)] <- 0

gene_cnt <- sapply(gene_groups_tem, function(g) sum(g %in% rownames(mat_z)))
gap_rows <- cumsum(gene_cnt)[-length(gene_cnt)]

png("CD8_TEM_Heatmap.png", width = 1600, height = 2000, res = 220)
pheatmap::pheatmap(
  mat_z, cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(-2, 2, length.out = 101),
  border_color = NA, main = "CD8 TEM: Pattern-ordered expression",
  fontsize_row = 10, fontsize_col = 12, angle_col = 45,
  gaps_row = gap_rows, cellwidth = 60
)
dev.off()
message("[Saved] CD8_TEM_Heatmap.png")

# =============================================================================
# SECTION 13: CD8 TEM ROS Gene DEG (Poor d7 vs d1, Good d7 vs d1)
# =============================================================================
message("\n>>> CD8 TEM ROS gene DEG analysis")

obj_tem2 <- subset(T_CD8, subset = predicted.celltype.major_fine == "CD8 TEM")
obj_tem2 <- safe_join_layers(obj_tem2)
Idents(obj_tem2) <- "Prognosis_Group"

ros_production <- c("NOX1","NOX2","NOX4","CYBB","NCF1","NCF2","NCF4",
                    "NOS2","PTGS2")
ros_scavenging <- c("SOD1","SOD2","SOD3","CAT","GPX1","GPX4",
                    "PRDX1","PRDX2","PRDX3","TXN","TXNRD1",
                    "GSR","GCLC","GCLM","NFE2L2")
ros_response   <- c("HMOX1","HSPA1A","HSPA1B","DDIT4","CDKN1A",
                    "TP53","ATF3","ATF4")

deg_poor_tem <- FindMarkers(obj_tem2, ident.1 = "Poor_d7", ident.2 = "Poor_d1",
                            logfc.threshold = 0, min.pct = 0.05,
                            test.use = "wilcox") %>%
  tibble::rownames_to_column("gene")

deg_good_tem <- FindMarkers(obj_tem2, ident.1 = "Good_d7", ident.2 = "Good_d1",
                            logfc.threshold = 0, min.pct = 0.05,
                            test.use = "wilcox") %>%
  tibble::rownames_to_column("gene")

print_ros <- function(deg, label) {
  cat("\n========", label, "========\n")
  cat("Positive = higher in d7; Negative = higher in d1\n")
  for (cat_name in c("ROS production","ROS scavenging","ROS response")) {
    genes_cat <- switch(cat_name,
                        "ROS production" = ros_production,
                        "ROS scavenging" = ros_scavenging,
                        "ROS response"   = ros_response)
    cat("\n---", cat_name, "---\n")
    deg %>% dplyr::filter(gene %in% genes_cat) %>%
      dplyr::select(gene, avg_log2FC, p_val_adj) %>%
      dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
      as.data.frame() %>% print()
  }
}

print_ros(deg_poor_tem, "Poor: d7 vs d1")
print_ros(deg_good_tem, "Good: d7 vs d1")

# =============================================================================
# SECTION 14: Cytotoxicity & Exhaustion Module Scores
#             (CD8 TCM+TEM and CD4 CTL)
# =============================================================================
message("\n>>> Cytotoxicity and exhaustion module scores")

SUBTYPE_COL  <- "predicted.celltype.l2"
GROUP_COL    <- "Prognosis_Group"
KEEP_CD8     <- c("CD8 TCM","CD8 TEM")

# CD8 subset
T_CD8_sub <- subset(T_CD8,
                    subset = predicted.celltype.l2 %in% KEEP_CD8)
T_CD8_sub$Prognosis_Group <- factor(
  trimws(as.character(T_CD8_sub$Prognosis_Group)), levels = GROUP_LEVELS
)

res8     <- add_module_scores_safe(T_CD8_sub, assay = "RNA", prefix = "Tcell")
T_CD8_sub <- res8$obj

df8_cyto_dot <- make_dot_df(T_CD8_sub, GROUP_COL, SUBTYPE_COL, res8$cyto_col) %>%
  dplyr::mutate(group   = factor(group, levels = GROUP_LEVELS),
                subtype = factor(subtype, levels = KEEP_CD8))
df8_exh_dot  <- make_dot_df(T_CD8_sub, GROUP_COL, SUBTYPE_COL, res8$exh_col) %>%
  dplyr::mutate(group   = factor(group, levels = GROUP_LEVELS),
                subtype = factor(subtype, levels = KEEP_CD8))

df8_cyto_box <- T_CD8_sub@meta.data %>% as.data.frame() %>%
  dplyr::transmute(group   = factor(as.character(Prognosis_Group),
                                    levels = GROUP_LEVELS),
                   subtype = as.character(.data[[SUBTYPE_COL]]),
                   score   = as.numeric(.data[[res8$cyto_col]])) %>%
  dplyr::filter(!is.na(group), subtype %in% KEEP_CD8, is.finite(score)) %>%
  dplyr::mutate(subtype = factor(subtype, levels = KEEP_CD8))

df8_exh_box <- T_CD8_sub@meta.data %>% as.data.frame() %>%
  dplyr::transmute(group   = factor(as.character(Prognosis_Group),
                                    levels = GROUP_LEVELS),
                   subtype = as.character(.data[[SUBTYPE_COL]]),
                   score   = as.numeric(.data[[res8$exh_col]])) %>%
  dplyr::filter(!is.na(group), subtype %in% KEEP_CD8, is.finite(score)) %>%
  dplyr::mutate(subtype = factor(subtype, levels = KEEP_CD8))

p8_dot_cyto <- plot_dot_score(df8_cyto_dot, "CD8 (TCM+TEM): Cytotoxicity")
p8_dot_exh  <- plot_dot_score(df8_exh_dot,  "CD8 (TCM+TEM): Exhaustion")
p8_box_cyto <- plot_box_score(df8_cyto_box, "CD8 (TCM+TEM): Cytotoxicity", "subtype")
p8_box_exh  <- plot_box_score(df8_exh_box,  "CD8 (TCM+TEM): Exhaustion",   "subtype")

save_fig(p8_dot_cyto, "CD8_TCM_TEM_Cytotoxicity_dot.png", 8.0, 3.2)
save_fig(p8_dot_exh,  "CD8_TCM_TEM_Exhaustion_dot.png",   8.0, 3.2)
save_fig(p8_box_cyto, "CD8_TCM_TEM_Cytotoxicity_box.png", 8.0, 4.8)
save_fig(p8_box_exh,  "CD8_TCM_TEM_Exhaustion_box.png",   8.0, 4.8)

# CD4 CTL subset
T_CD4_CTL <- subset(T_CD4, subset = predicted.celltype.l2 == "CD4 CTL")
T_CD4_CTL$Prognosis_Group <- factor(
  trimws(as.character(T_CD4_CTL$Prognosis_Group)), levels = GROUP_LEVELS
)

res4      <- add_module_scores_safe(T_CD4_CTL, assay = "RNA", prefix = "Tcell")
T_CD4_CTL <- res4$obj

df4_cyto_dot <- make_dot_df(T_CD4_CTL, GROUP_COL, SUBTYPE_COL, res4$cyto_col) %>%
  dplyr::mutate(group = factor(group, levels = GROUP_LEVELS))

df4_cyto_box <- T_CD4_CTL@meta.data %>% as.data.frame() %>%
  dplyr::transmute(group   = factor(as.character(Prognosis_Group),
                                    levels = GROUP_LEVELS),
                   subtype = as.character(.data[[SUBTYPE_COL]]),
                   score   = as.numeric(.data[[res4$cyto_col]])) %>%
  dplyr::filter(!is.na(group), is.finite(score))

p4_dot_cyto <- plot_dot_score(df4_cyto_dot, "CD4 CTL: Cytotoxicity")
p4_box_cyto <- plot_box_score(df4_cyto_box, "CD4 CTL: Cytotoxicity")

save_fig(p4_dot_cyto, "CD4_CTL_Cytotoxicity_dot.png", 7.2, 3.0)
save_fig(p4_box_cyto, "CD4_CTL_Cytotoxicity_box.png", 6.4, 4.2)

message("\n>>> CD4/CD8 T cell analysis complete.")