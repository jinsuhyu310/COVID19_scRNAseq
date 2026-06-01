# =============================================================================
# 05_Monocyte_Integration_Analysis.R
# Monocyte Integration: SCT integration → IL6/TNF/Glycolysis/Hypoxia module
# scores → UMAP visualization → cluster characterization → violin plots →
# heatmap → volcano → cluster overlap
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
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(rstatix)
library(KEGGREST)
library(AnnotationDbi)
library(future)

future::plan("sequential")

# -----------------------------------------------------------------------------
# 0. Shared settings
# -----------------------------------------------------------------------------
GROUP_LEVELS <- c("Poor_d1", "Poor_d7", "Good_d1", "Good_d7")
PADJ_CUT     <- 0.05

my_pub_colors <- c(
  "Good_d1" = "#80B1D3", "Good_d7" = "#483D8B",
  "Poor_d1" = "#FDB462", "Poor_d7" = "#D95F02"
)

COLOR_LOW  <- "#08306B"
COLOR_MID  <- "#FFFFFF"
COLOR_HIGH <- "#99000D"

get_md_vec <- function(seu, col) {
  x <- seu@meta.data[[col]]
  if (is.data.frame(x)) x <- x[[1]]
  trimws(as.character(x))
}

# Auto-detect UMAP reduction
get_umap_name <- function(obj) {
  if ("ref.umap" %in% names(obj@reductions)) "ref.umap"
  else if ("umap" %in% names(obj@reductions)) "umap"
  else stop("No UMAP reduction found (ref.umap / umap).")
}

# Auto-detect patient ID column
get_patient_col <- function(obj) {
  candidates <- c("Patient_ID", "patient", "Patient", "orig.ident",
                  "SampleID", "donor", "Donor", "ID", "Sample")
  col <- candidates[candidates %in% colnames(obj@meta.data)][1]
  if (is.na(col)) stop("No patient ID column found in meta.data.")
  col
}

# Auto-detect seurat cluster column
get_cluster_col <- function(md) {
  if ("seurat_clusters" %in% colnames(md)) "seurat_clusters"
  else if ("seurat_cluster" %in% colnames(md)) "seurat_cluster"
  else stop("No seurat cluster column found.")
}

# =============================================================================
# SECTION 1: Helper Functions
# =============================================================================

# -----------------------------------------------------------------------------
# 1-A. IL6 counts extraction across all count layers (pre-integration)
# -----------------------------------------------------------------------------
get_il6_counts <- function(seu, gene = "IL6") {
  DefaultAssay(seu) <- "RNA"
  lyr           <- tryCatch(SeuratObject::Layers(seu[["RNA"]]),
                            error = function(e) character(0))
  counts_layers <- lyr[grepl("^counts", lyr)]
  
  out <- rep(0, ncol(seu)); names(out) <- colnames(seu)
  
  if (length(counts_layers) == 0) {
    mat0 <- GetAssayData(seu, assay = "RNA", slot = "counts")
    if (!gene %in% rownames(mat0)) stop("Gene ", gene, " not in RNA counts.")
    v <- as.numeric(mat0[gene, ]); names(v) <- colnames(seu)
    out[names(v)] <- v
    return(out)
  }
  for (lay in counts_layers) {
    mat <- SeuratObject::LayerData(seu, assay = "RNA", layer = lay)
    if (!gene %in% rownames(mat)) next
    v <- as.numeric(mat[gene, ]); names(v) <- colnames(mat)
    out[names(v)] <- pmax(out[names(v)], v)
  }
  out
}

# -----------------------------------------------------------------------------
# 1-B. Module score UMAP (raw gradient, 4-panel split)
# -----------------------------------------------------------------------------
build_umap_md <- function(obj) {
  umap_name <- get_umap_name(obj)
  emb <- Embeddings(obj, reduction = umap_name)
  md  <- as.data.frame(obj@meta.data)
  md$cell   <- rownames(md)
  md$UMAP_1 <- emb[, 1]
  md$UMAP_2 <- emb[, 2]
  md$Prognosis_Group <- factor(as.character(md$Prognosis_Group),
                               levels = GROUP_LEVELS)
  md
}

plot_umap_module_4groups <- function(md, feature, title_txt, out_png,
                                     ncol = 2, width = 7, height = 7) {
  xlim_use <- range(md$UMAP_1, na.rm = TRUE)
  ylim_use <- range(md$UMAP_2, na.rm = TRUE)
  
  v <- suppressWarnings(as.numeric(md[[feature]]))
  z <- as.numeric(scale(v))
  md$Value <- pmax(-1.5, pmin(1.5, z))
  
  p <- ggplot(md, aes(UMAP_1, UMAP_2, color = Value)) +
    geom_point(size = 0.35, alpha = 0.95) +
    facet_wrap(~ Prognosis_Group, ncol = ncol) +
    coord_fixed(xlim = xlim_use, ylim = ylim_use) +
    scale_color_viridis_c(option = "magma", direction = -1,
                          limits = c(-1.5, 1.5),
                          oob = scales::squish, name = "z-score") +
    theme_classic(base_size = 16) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.border = element_rect(color = "black", fill = NA,
                                      linewidth = 0.6),
          legend.title = element_text(face = "bold"),
          plot.title   = element_text(face = "bold", hjust = 0.5)) +
    labs(x = "UMAP 1", y = "UMAP 2", title = title_txt)
  
  ggsave(out_png, p, width = width, height = height, dpi = 300)
  message("Saved: ", out_png)
  invisible(p)
}

# -----------------------------------------------------------------------------
# 1-C. Delta-score UMAP (d7 - d1 per cluster, Poor vs Good side by side)
# -----------------------------------------------------------------------------
calc_delta_by_cluster <- function(md, prog_label,
                                  group_col   = "Prognosis_Group",
                                  cluster_col = "seurat_clusters",
                                  feature) {
  g1 <- paste0(prog_label, "_d1")
  g7 <- paste0(prog_label, "_d7")
  
  df <- md %>%
    dplyr::filter(.data[[group_col]] %in% c(g1, g7)) %>%
    dplyr::mutate(Time    = ifelse(.data[[group_col]] == g1, "d1", "d7"),
                  cluster = as.factor(.data[[cluster_col]]))
  
  delta_tbl <- df %>%
    dplyr::group_by(cluster, Time) %>%
    dplyr::summarise(mean_score = mean(.data[[feature]], na.rm = TRUE),
                     .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Time, values_from = mean_score) %>%
    dplyr::mutate(delta = d7 - d1)
  
  df2 <- df %>%
    dplyr::left_join(delta_tbl %>% dplyr::select(cluster, delta),
                     by = "cluster")
  
  centers <- df2 %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE),
                     UMAP_2 = median(UMAP_2, na.rm = TRUE),
                     .groups = "drop")
  
  list(md_prog = df2, delta_by_cluster = delta_tbl, cluster_centers = centers)
}

make_delta_umap_plot <- function(md_prog, cluster_centers,
                                 subtitle_txt = NULL,
                                 xlim_use, ylim_use,
                                 show_legend = TRUE) {
  ggplot(md_prog, aes(UMAP_1, UMAP_2, color = delta)) +
    geom_point(size = 0.5, stroke = 0) +
    coord_fixed(xlim = xlim_use, ylim = ylim_use) +
    scale_color_gradient2(low = COLOR_LOW, mid = COLOR_MID, high = COLOR_HIGH,
                          midpoint = 0, oob = squish,
                          name = expression(bold(Delta ~ "Score"))) +
    ggrepel::geom_label_repel(
      data = cluster_centers,
      aes(UMAP_1, UMAP_2, label = cluster),
      inherit.aes = FALSE, size = 4, fontface = "bold",
      label.size = 0.2, color = "black",
      fill = alpha("white", 0.85), box.padding = 0.5,
      max.overlaps = Inf, seed = 42
    ) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(), axis.title = element_blank(),
          plot.subtitle     = element_text(hjust = 0.5, size = 13,
                                           color = "grey30"),
          legend.position   = if (show_legend) "right" else "none",
          legend.title      = element_text(size = 13),
          legend.text       = element_text(size = 11),
          legend.key.height = unit(1.5, "cm"),
          panel.border = element_rect(color = "black", fill = NA,
                                      linewidth = 1)) +
    labs(subtitle = subtitle_txt)
}

save_delta_pair_plot <- function(res_poor, res_good,
                                 xlim_use, ylim_use,
                                 main_title, out_png,
                                 width = 10, height = 5) {
  p_poor <- make_delta_umap_plot(res_poor$md_prog, res_poor$cluster_centers,
                                 "Poor Prognosis (Day 7 vs. Day 1)",
                                 xlim_use, ylim_use, show_legend = TRUE)
  p_good <- make_delta_umap_plot(res_good$md_prog, res_good$cluster_centers,
                                 "Good Prognosis (Day 7 vs. Day 1)",
                                 xlim_use, ylim_use, show_legend = TRUE)
  
  p_combined <- (p_poor | (p_good + guides(color = "none"))) +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "right")
  
  p_combined <- p_combined +
    plot_annotation(
      title = main_title,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 18))
    )
  
  ggsave(out_png, p_combined, width = width, height = height,
         dpi = 300, bg = "white")
  message("Saved: ", out_png)
  invisible(p_combined)
}

# -----------------------------------------------------------------------------
# 1-D. Cluster violin plot (z-score, paired patient lines, Wilcoxon stats)
# -----------------------------------------------------------------------------
make_cluster_violin <- function(obj, cluster_id, score_col, module_label,
                                comparisons, out_png,
                                title_txt = NULL, width = 6.2, height = 5.2) {
  
  patient_col <- get_patient_col(obj)
  fill_map    <- my_pub_colors[GROUP_LEVELS]
  
  obj_sub <- subset(obj, subset = seurat_clusters == cluster_id)
  
  df_long <- obj_sub@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::transmute(
      cell_id,
      group       = as.character(Prognosis_Group),
      patient_raw = as.character(.data[[patient_col]]),
      score       = suppressWarnings(as.numeric(.data[[score_col]]))
    ) %>%
    dplyr::filter(!is.na(group), !is.na(patient_raw), is.finite(score)) %>%
    dplyr::mutate(
      group        = factor(group, levels = GROUP_LEVELS),
      Prognosis    = sub("_(d1|d7)$", "", group),
      Timepoint    = sub("^.*_(d1|d7)$", "\\1", group),
      Prognosis    = factor(Prognosis, levels = c("Poor","Good")),
      Timepoint    = factor(Timepoint, levels = c("d1","d7")),
      patient_base = sub("-(d1|d7)$", "", patient_raw),
      mu   = mean(score, na.rm = TRUE),
      sig  = sd(score,  na.rm = TRUE),
      z    = ifelse(is.na(sig) | sig == 0, NA_real_,
                    (score - mu) / sig)
    ) %>%
    dplyr::filter(is.finite(z)) %>%
    dplyr::mutate(z_cap  = pmax(-1, pmin(1, z)),
                  module = module_label)
  
  df_pat <- df_long %>%
    dplyr::group_by(Prognosis, Timepoint, patient_base) %>%
    dplyr::summarise(pat_mean = mean(z_cap, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(is.finite(pat_mean)) %>%
    dplyr::mutate(
      group  = factor(paste0(Prognosis, "_", Timepoint), levels = GROUP_LEVELS),
      module = module_label
    )
  
  # Wilcoxon (unpaired) for each comparison
  get_p <- function(g1, g2) {
    x <- df_long$z_cap[df_long$group == g1]; x <- x[is.finite(x)]
    y <- df_long$z_cap[df_long$group == g2]; y <- y[is.finite(y)]
    if (length(x) < 2 || length(y) < 2) return(NA_real_)
    suppressWarnings(wilcox.test(x, y, paired = FALSE)$p.value)
  }
  
  stat_tbl <- do.call(rbind, lapply(comparisons, function(cc) {
    data.frame(group1 = cc[1], group2 = cc[2],
               p = get_p(cc[1], cc[2]), stringsAsFactors = FALSE)
  })) %>%
    dplyr::mutate(
      p.adj = p.adjust(p, method = "BH"),
      p.signif = dplyr::case_when(
        is.na(p.adj)    ~ "NA",
        p.adj <= 1e-4   ~ "****",
        p.adj <= 1e-3   ~ "***",
        p.adj <= 1e-2   ~ "**",
        p.adj <= 5e-2   ~ "*",
        TRUE            ~ "ns"
      ),
      module = module_label
    )
  
  ymax <- max(df_long$z_cap, na.rm = TRUE)
  stat_tbl <- stat_tbl %>%
    dplyr::mutate(y.position = ymax + 1e-6 + (seq_len(dplyr::n()) * 0.13))
  ylim_top <- max(stat_tbl$y.position) + 0.01
  
  p <- ggplot(df_long, aes(x = group, y = z_cap, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.9) +
    geom_boxplot(width = 0.15, outlier.size = 0.1, alpha = 0.9,
                 fill = "white") +
    geom_line(data = df_pat,
              aes(x = group, y = pat_mean,
                  group = interaction(Prognosis, patient_base)),
              inherit.aes = FALSE, color = "grey40",
              linewidth = 0.5, alpha = 0.7) +
    geom_point(data = df_pat, aes(x = group, y = pat_mean),
               inherit.aes = FALSE, color = "grey20", size = 1.2, alpha = 0.8) +
    ggpubr::stat_pvalue_manual(stat_tbl, label = "p.signif",
                               xmin = "group1", xmax = "group2",
                               y.position = "y.position",
                               tip.length = 0.005, bracket.size = 0.3,
                               size = 4.5, vjust = 0.4) +
    scale_fill_manual(values = fill_map) +
    facet_wrap(~ module, ncol = 1, scales = "fixed") +
    coord_cartesian(ylim = c(-1.0, ylim_top)) +
    labs(x = NULL, y = "Module score (z-score)",
         title = title_txt %||% paste0("Monocytes (C", cluster_id, ")")) +
    theme_bw(base_size = 14) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1,
                                          color = "black"),
          axis.text.y      = element_text(color = "black"),
          strip.background = element_rect(fill = "grey95"),
          strip.text       = element_text(face = "bold"),
          plot.title       = element_text(face = "bold", hjust = 0.5),
          panel.grid       = element_blank(),
          legend.position  = "none")
  
  ggsave(out_png, p, width = width, height = height, dpi = 300, bg = "white")
  message("Saved: ", out_png)
  invisible(p)
}

# `%||%` operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------------------------------------------------------
# 1-E. Ensure module score column exists
# -----------------------------------------------------------------------------
ensure_module_score <- function(obj, col_name, gene_list, prefix) {
  if (col_name %in% colnames(obj@meta.data)) {
    message("'", col_name, "' already exists. Skipping.")
    return(obj)
  }
  genes_use <- intersect(gene_list, rownames(obj))
  cat(sprintf("[%s] Genes in object: %d / %d\n",
              col_name, length(genes_use), length(gene_list)))
  if (length(genes_use) < 5) stop("Too few genes (<5) for ", col_name)
  tmp_name <- paste0(prefix, "_tmp_")
  obj <- AddModuleScore(obj, features = list(genes_use),
                        name = tmp_name, assay = "RNA", search = FALSE)
  colnames(obj@meta.data)[
    colnames(obj@meta.data) == paste0(tmp_name, "1")
  ] <- col_name
  message("Module score '", col_name, "' added.")
  obj
}

# Helper to get KEGG pathway gene symbols
get_kegg_symbols <- function(kegg_id) {
  kk <- tryCatch(KEGGREST::keggLink("hsa", paste0("path:", kegg_id)),
                 error = function(e) NULL)
  if (is.null(kk)) return(character(0))
  entrez <- unique(sub("^hsa:", "", c(names(kk), unname(kk))))
  entrez <- entrez[grepl("^\\d+$", entrez)]
  syms <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = entrez,
                                column = "SYMBOL", keytype = "ENTREZID",
                                multiVals = "first")
  unique(na.omit(unname(syms)))
}

# =============================================================================
# SECTION 2: Monocyte Subset & SCT Integration
# =============================================================================
message("\n>>> Subsetting monocytes and running SCT integration")

obj <- merged_int
DefaultAssay(obj) <- "RNA"
suppressWarnings({ obj <- JoinLayers(obj, assay = "RNA") })
obj@meta.data[["Prognosis_Group"]] <- factor(
  get_md_vec(obj, "Prognosis_Group"), levels = GROUP_LEVELS
)

mono <- subset(obj,
               subset = predicted.celltype.major %in% c("CD14 Mono", "CD16 Mono"))
DefaultAssay(mono) <- "RNA"
suppressWarnings({ mono <- JoinLayers(mono, assay = "RNA") })

message("Monocyte cell counts by group:")
print(table(mono@meta.data[["Prognosis_Group"]], useNA = "ifany"))

# IL6 counts before integration (for later mapping)
il6_counts_pre <- get_il6_counts(mono, "IL6")
message("Total IL6+ (counts>0) cells (pre-integration): ",
        sum(il6_counts_pre > 0, na.rm = TRUE))

# SCT integration per sample
vars_regress <- intersect(c("percent.mt", "nCount_RNA"), colnames(mono@meta.data))
mono_list <- SplitObject(mono, split.by = "orig.ident")
mono_list <- mono_list[sapply(mono_list, ncol) > 0]
stopifnot(length(mono_list) >= 2)

message("Split sizes by orig.ident:")
print(sapply(mono_list, ncol))

mono_list <- lapply(mono_list, function(x) {
  DefaultAssay(x) <- "RNA"
  SCTransform(x, assay = "RNA", new.assay.name = "SCT",
              vars.to.regress = if (length(vars_regress) > 0) vars_regress else NULL,
              verbose = FALSE)
})

features  <- SelectIntegrationFeatures(object.list = mono_list, nfeatures = 2000)
mono_list <- PrepSCTIntegration(object.list = mono_list,
                                anchor.features = features, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = mono_list,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  dims = 1:20, reduction = "rpca",
                                  k.anchor = 20)

n_anchor     <- tryCatch(nrow(anchors@anchors), error = function(e) NA_integer_)
k_weight_use <- max(20, min(100, floor(n_anchor / 2)))
message("Using k.weight = ", k_weight_use)

mono_int <- IntegrateData(anchorset = anchors,
                          normalization.method = "SCT",
                          dims = 1:20, k.weight = k_weight_use)

# Dimensionality reduction & clustering
DefaultAssay(mono_int) <- "integrated"
set.seed(1234)
mono_int <- ScaleData(mono_int, verbose = FALSE)
mono_int <- RunPCA(mono_int, npcs = 50, verbose = FALSE)
mono_int <- FindNeighbors(mono_int, dims = 1:20)
mono_int <- FindClusters(mono_int, resolution = 0.6)
mono_int <- RunUMAP(mono_int, dims = 1:20, reduction = "pca", verbose = FALSE)
DefaultAssay(mono_int) <- "RNA"
suppressWarnings({ mono_int <- JoinLayers(mono_int, assay = "RNA") })

message("Integration complete. Cluster distribution:")
print(table(mono_int$seurat_clusters))

# =============================================================================
# SECTION 3: IL6+ Cell UMAP
# =============================================================================
message("\n>>> IL6+ cell UMAP")

md <- build_umap_md(mono_int)
md$IL6_counts <- il6_counts_pre[md$cell]
md$IL6_pos    <- is.finite(md$IL6_counts) & (md$IL6_counts > 0)

message("IL6+ monocyte cells by group:")
print(table(md$group[md$IL6_pos], useNA = "ifany"))

xlim_use <- range(md$UMAP_1, na.rm = TRUE)
ylim_use <- range(md$UMAP_2, na.rm = TRUE)

p_il6_umap <- ggplot() +
  geom_point(data = md, aes(UMAP_1, UMAP_2),
             color = "grey70", alpha = 0.55, size = 0.55, shape = 16) +
  geom_point(data = md %>% dplyr::filter(IL6_pos),
             aes(UMAP_1, UMAP_2),
             color = "#DC0000", alpha = 0.95, size = 0.95, shape = 16) +
  facet_wrap(~ Prognosis_Group, ncol = 2, drop = FALSE) +
  labs(title = expression(bold("Monocytes: ") * italic("IL6")^"+" *
                            bold(" cells (counts>0)")),
       x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 16) +
  theme(plot.title       = element_text(hjust = 0.5, face = "bold",
                                        size = 14, margin = margin(b = 10)),
        strip.text       = element_text(face = "bold", size = 11),
        strip.background = element_blank(),
        axis.title       = element_text(face = "bold", size = 11),
        axis.text        = element_text(color = "black", size = 10),
        panel.grid       = element_blank(),
        legend.position  = "none")

ggsave("Monocytes_IL6pos_UMAP.png", p_il6_umap,
       width = 5, height = 6.5, dpi = 600)
message("Saved: Monocytes_IL6pos_UMAP.png")

# =============================================================================
# SECTION 4: Module Score Calculation
# =============================================================================
message("\n>>> Calculating module scores")

msig_h  <- msigdbr(species = "Homo sapiens", category = "H")
msig_go <- msigdbr(species = "Homo sapiens", category = "C5",
                   subcategory = "GO:BP")

get_hallmark <- function(gs) {
  msig_h %>% dplyr::filter(gs_name == gs) %>%
    dplyr::pull(gene_symbol) %>% unique()
}
get_gobp <- function(gs) {
  msig_go %>% dplyr::filter(gs_name == gs) %>%
    dplyr::pull(gene_symbol) %>% unique()
}

# IL6-JAK-STAT3
mono_int <- ensure_module_score(
  mono_int, "IL6_JAK_STAT1",
  get_hallmark("HALLMARK_IL6_JAK_STAT3_SIGNALING"),
  "IL6_JAK_STAT"
)

# TNFa-NFkB (AddModuleScore appends "1" to name)
if (!"TNFA_NFKB1" %in% colnames(mono_int@meta.data)) {
  mono_int <- ensure_module_score(
    mono_int, "TNFA_NFKB1",
    get_hallmark("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
    "TNFA_NFKB"
  )
}

# Glycolysis (KEGG hsa00010)
mono_int <- ensure_module_score(
  mono_int, "Glycolysis_KEGG_hsa00010",
  get_kegg_symbols("hsa00010"),
  "Glycolysis_KEGG_hsa00010"
)

# Hypoxia
mono_int <- ensure_module_score(
  mono_int, "Hypoxia_HALLMARK",
  get_hallmark("HALLMARK_HYPOXIA"),
  "Hypoxia_HALLMARK"
)

cat("Module scores added to meta.data.\n")

# =============================================================================
# SECTION 5: Module Score UMAPs (z-score, 4-panel)
# =============================================================================
message("\n>>> Module score UMAPs")

md <- build_umap_md(mono_int)

plot_umap_module_4groups(md, "IL6_JAK_STAT1",
                         "Monocytes: IL6-JAK-STAT3 signaling",
                         "mono_int_UMAP_IL6_JAK_STAT3.png")

plot_umap_module_4groups(md, "TNFA_NFKB1",
                         "Monocytes: TNF\u03b1-NF-\u03baB signaling",
                         "mono_int_UMAP_TNFA_NFKB1.png")

plot_umap_module_4groups(md, "Glycolysis_KEGG_hsa00010",
                         "Monocytes: Glycolysis",
                         "mono_int_UMAP_Glycolysis.png")

# =============================================================================
# SECTION 6: Delta Score UMAPs (d7 - d1 per cluster)
# =============================================================================
message("\n>>> Delta score UMAPs")

md <- build_umap_md(mono_int)
cluster_col <- get_cluster_col(md)
xlim_use    <- range(md$UMAP_1, na.rm = TRUE)
ylim_use    <- range(md$UMAP_2, na.rm = TRUE)

for (feat_info in list(
  list(feature = "IL6_JAK_STAT1",          title = "Changes in IL6-JAK-STAT3 Signaling",
       out = "Monocyte_UMAP_IL6_Delta_Score_Poor_vs_Good.png"),
  list(feature = "TNFA_NFKB1",             title = "Changes in TNF\u03b1-NF-\u03baB Signaling",
       out = "Monocyte_UMAP_TNF_Delta_Score_Poor_vs_Good.png"),
  list(feature = "Glycolysis_KEGG_hsa00010", title = "Temporal Dynamics: Glycolysis (KEGG hsa00010)",
       out = "Monocyte_UMAP_Glycolysis_Delta_Score_Poor_vs_Good.png"),
  list(feature = "Hypoxia_HALLMARK",        title = "Temporal Dynamics: Hypoxia (HALLMARK)",
       out = "Monocyte_UMAP_Hypoxia_Delta_Score_Poor_vs_Good.png")
)) {
  res_poor <- calc_delta_by_cluster(md, "Poor", "Prognosis_Group",
                                    cluster_col, feat_info$feature)
  res_good <- calc_delta_by_cluster(md, "Good", "Prognosis_Group",
                                    cluster_col, feat_info$feature)
  
  save_delta_pair_plot(res_poor, res_good, xlim_use, ylim_use,
                       feat_info$title, feat_info$out)
  
  cat("\n=== Poor delta by cluster (", feat_info$feature, ") ===\n", sep = "")
  print(res_poor$delta_by_cluster %>% dplyr::arrange(dplyr::desc(delta)) %>%
          dplyr::select(cluster, d1, d7, delta), n = 30)
  cat("\n=== Good delta by cluster ===\n")
  print(res_good$delta_by_cluster %>% dplyr::arrange(dplyr::desc(delta)) %>%
          dplyr::select(cluster, d1, d7, delta), n = 30)
}

# Top cluster showing greatest IL6 signaling decrease in Good prognosis
cl_delta_good <- calc_delta_by_cluster(md, "Good", "Prognosis_Group",
                                       cluster_col, "IL6_JAK_STAT1")
top1 <- cl_delta_good$delta_by_cluster %>%
  dplyr::filter(pmin(d1, d7, na.rm = TRUE) > 0) %>%
  dplyr::arrange(delta) %>% head(1) %>% dplyr::pull(cluster)
message("Top IL6-decreasing cluster (Good): ", top1)

mono_int$top_cluster_highlight <- factor(
  ifelse(as.character(mono_int$seurat_clusters) == as.character(top1),
         "Top1_cluster", "Other"),
  levels = c("Other", "Top1_cluster")
)

p_top1 <- DimPlot(mono_int, reduction = get_umap_name(mono_int),
                  group.by = "top_cluster_highlight",
                  cols = c(Other = "grey90", Top1_cluster = "#4A148C"),
                  pt.size = 0.12) +
  ggtitle(paste0("Top IL6 signaling-decreasing cluster\n(Good prognosis): C", top1)) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none")

ggsave("mono_int_GoodTop1DeltaCluster_UMAP.png", p_top1,
       width = 5, height = 6, dpi = 300, bg = "white")

# =============================================================================
# SECTION 7: Cluster 7 GO BP Dotplot
# =============================================================================
message("\n>>> GO BP dotplot: Cluster 7 (Poor_d7)")

mono_int$cluster_group <- paste0(mono_int$seurat_clusters, "_",
                                 mono_int$Prognosis_Group)
Idents(mono_int) <- "cluster_group"

target_c7     <- "7_Poor_d7"
other_c7      <- setdiff(levels(Idents(mono_int)), target_c7)

deg_c7 <- FindMarkers(mono_int, ident.1 = target_c7, ident.2 = other_c7,
                      min.pct = 0.1, logfc.threshold = 0.25,
                      test.use = "wilcox") %>%
  tibble::rownames_to_column("gene")

universe_c7  <- deg_c7$gene
genes_up_c7  <- deg_c7 %>% dplyr::filter(p_val_adj < PADJ_CUT, avg_log2FC > 1) %>%
  dplyr::pull(gene)
genes_dn_c7  <- deg_c7 %>% dplyr::filter(p_val_adj < PADJ_CUT, avg_log2FC < -1) %>%
  dplyr::pull(gene)

run_go_bp_with_simplify <- function(genes, universe, simplify_cut = 0.7) {
  if (length(genes) == 0) return(NULL)
  ego <- tryCatch(
    enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
             keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
             pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
    error = function(e) NULL
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  tryCatch(clusterProfiler::simplify(ego, cutoff = simplify_cut,
                                     by = "p.adjust", select_fun = min),
           error = function(e) ego)
}

ego_up_c7 <- run_go_bp_with_simplify(genes_up_c7, universe_c7)
ego_dn_c7 <- run_go_bp_with_simplify(genes_dn_c7, universe_c7)

get_top_go_terms <- function(ego, n, reg) {
  if (is.null(ego)) return(NULL)
  as.data.frame(ego) %>% dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust) %>% head(n) %>%
    dplyr::mutate(regulation = reg, Description = as.character(Description))
}

go_c7 <- dplyr::bind_rows(
  get_top_go_terms(ego_up_c7, 15, "up"),
  get_top_go_terms(ego_dn_c7,  5, "down")
)

if (!is.null(go_c7) && nrow(go_c7) > 0) {
  term_order <- go_c7 %>%
    dplyr::group_by(Description, regulation) %>%
    dplyr::summarise(min_p = min(p.adjust), .groups = "drop") %>%
    dplyr::arrange(regulation, min_p) %>% dplyr::pull(Description)
  
  go_c7 <- go_c7 %>%
    dplyr::mutate(
      Term_wrapped = factor(stringr::str_wrap(Description, 40),
                            levels = rev(stringr::str_wrap(term_order, 40))),
      regulation   = factor(regulation, levels = c("down","up"))
    )
  
  p_go_c7 <- ggplot(go_c7, aes(x = regulation, y = Term_wrapped)) +
    geom_point(aes(size = Count, colour = regulation)) +
    scale_colour_manual(values = c(down = "red", up = "blue"), name = "") +
    scale_size_continuous(range = c(2.5, 8), name = "Count") +
    labs(x = NULL, y = NULL,
         title = "GO BP: Monocyte Cluster 7 in Poor_d7 (p<0.05)") +
    theme_bw(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 17, face = "bold"),
          axis.text.x = element_text(size = 17, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank())
  
  out_go_c7 <- "MONO_C7_Poor_d7_vs_all_GO_BP_dotplot.png"
  png(out_go_c7, width = 3500, height = 3300, res = 300)
  print(p_go_c7); dev.off()
  message("Saved: ", out_go_c7)
}

# =============================================================================
# SECTION 8: Curated Inflammation Panel Dotplot (all clusters)
# =============================================================================
message("\n>>> Curated inflammation panel dotplot")

panel <- tibble::tribble(
  ~set_id,                                        ~set_type, ~label,
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",             "H", "IL6-JAK-STAT3",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",             "H", "TNFa-NFkB",
  "GOBP_INFLAMMATORY_RESPONSE",                   "G", "Inflammatory Response",
  "GOBP_INNATE_IMMUNE_RESPONSE",                  "G", "Innate Immune Response",
  "GOBP_CYTOKINE_PRODUCTION",                     "G", "Cytokine Production",
  "GOBP_CHEMOKINE_MEDIATED_SIGNALING_PATHWAY",    "G", "Chemokine Signaling",
  "GOBP_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",    "G", "TLR Signaling",
  "GOBP_NEUTROPHIL_ACTIVATION",                   "G", "Neutrophil Activation",
  "GOBP_COMPLEMENT_ACTIVATION",                   "G", "Complement Activation",
  "GOBP_RESPONSE_TO_TYPE_I_INTERFERON",           "G", "Type I IFN Response"
)

get_geneset <- function(sid, stype) {
  if (stype == "H") msig_h %>% dplyr::filter(gs_name == sid) %>%
    dplyr::pull(gene_symbol) %>% unique()
  else msig_go %>% dplyr::filter(gs_name == sid) %>%
    dplyr::pull(gene_symbol) %>% unique()
}

features_list <- setNames(
  lapply(seq_len(nrow(panel)),
         function(i) get_geneset(panel$set_id[i], panel$set_type[i])),
  panel$label
)

# Filter gene sets with enough genes
keep <- sapply(features_list, length) >= 10
features_list <- features_list[keep]
panel         <- panel[keep, ]

score_prefix <- "CurPanel_"
Idents(mono_int) <- "seurat_clusters"
mono_int <- AddModuleScore(mono_int, features = features_list,
                           name = score_prefix)

score_cols <- grep(paste0("^", score_prefix, "\\d+$"),
                   colnames(mono_int@meta.data), value = TRUE)
score_map  <- tibble::tibble(score_col = score_cols,
                             label = names(features_list))

md_panel <- mono_int@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::select(cell, seurat_clusters, dplyr::all_of(score_cols)) %>%
  dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(score_cols),
                      names_to = "score_col", values_to = "score") %>%
  dplyr::left_join(score_map, by = "score_col") %>%
  dplyr::group_by(seurat_clusters, label) %>%
  dplyr::summarise(mean_score = mean(score, na.rm = TRUE),
                   pct_pos    = mean(score > 0, na.rm = TRUE),
                   .groups = "drop")

label_levels <- score_map$label
md_panel$label <- factor(md_panel$label, levels = rev(label_levels))

cl_levels <- sort(unique(md_panel$seurat_clusters))
suppressWarnings({
  cl_num <- as.numeric(cl_levels)
  if (all(!is.na(cl_num))) cl_levels <- as.character(sort(cl_num))
})
md_panel$seurat_clusters <- factor(md_panel$seurat_clusters, levels = cl_levels)

p_panel <- ggplot(md_panel, aes(x = seurat_clusters, y = label)) +
  geom_point(aes(size = pct_pos, colour = mean_score), alpha = 0.95) +
  scale_size_continuous(name = "% Cells (Score > 0)",
                        labels = scales::percent_format(accuracy = 1),
                        range = c(1.5, 9), limits = c(0, 1)) +
  scale_colour_viridis_c(option = "magma", direction = -1,
                         name = "Mean Module Score") +
  labs(x = "Seurat Cluster", y = NULL,
       title = "Curated Inflammation Panel: Monocyte Clusters") +
  theme_classic(base_size = 16) +
  theme(plot.title  = element_text(face = "bold", hjust = 0.5, size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        legend.title = element_text(face = "bold"))

plot_h <- max(2400, 180 * nlevels(md_panel$label) + 1200)
png("MONO_curatedInflammationPanel_by_cluster_dotplot.png",
    width = 4200, height = plot_h, res = 300)
print(p_panel); dev.off()
message("Saved: MONO_curatedInflammationPanel_by_cluster_dotplot.png")

# =============================================================================
# SECTION 9: IL6 Expression Violin (IL6+ monocytes only)
# =============================================================================
message("\n>>> IL6 expression violin (IL6+ cells)")

# Rebuild il6_counts from mono_int layers
mono_int@meta.data[["Prognosis_Group"]] <- factor(
  trimws(as.character(mono_int@meta.data[["Prognosis_Group"]])),
  levels = GROUP_LEVELS
)

il6_counts_int <- get_il6_counts(mono_int, "IL6")

dfA <- FetchData(mono_int, vars = c("IL6", "Prognosis_Group",
                                    "predicted.celltype.major")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    group    = factor(trimws(Prognosis_Group), levels = GROUP_LEVELS),
    celltype = trimws(predicted.celltype.major),
    IL6      = suppressWarnings(as.numeric(IL6)),
    IL6_counts = il6_counts_int[cell],
    IL6_pos    = is.finite(IL6_counts) & IL6_counts > 0
  ) %>%
  dplyr::filter(!is.na(group),
                celltype %in% c("CD14 Mono","CD16 Mono"),
                IL6_pos)

cat("\nMonocyte IL6+ cell counts by group:\n")
print(table(dfA$group, useNA = "ifany"))

comparisons_il6 <- list(c("Poor_d1","Poor_d7"),
                        c("Poor_d1","Good_d1"),
                        c("Poor_d7","Good_d1"))

stat_il6 <- dfA %>%
  rstatix::pairwise_wilcox_test(IL6 ~ group, p.adjust.method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::inner_join(
    tibble::tibble(group1 = sapply(comparisons_il6, `[`, 1),
                   group2 = sapply(comparisons_il6, `[`, 2),
                   order  = seq_along(comparisons_il6)),
    by = c("group1","group2")
  ) %>%
  dplyr::arrange(order) %>%
  dplyr::mutate(y.position = max(dfA$IL6, na.rm = TRUE) +
                  (0.08 + (order - 1) * 0.10) *
                  diff(range(dfA$IL6, na.rm = TRUE)))

pA <- ggplot(dfA, aes(x = group, y = IL6, fill = group)) +
  geom_violin(trim = TRUE, scale = "width", linewidth = 0.25) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.35,
               alpha = 1, fill = "white", color = "black") +
  scale_fill_manual(values = my_pub_colors, drop = FALSE) +
  ggpubr::stat_pvalue_manual(stat_il6, label = "p.adj.signif",
                             xmin = "group1", xmax = "group2",
                             y.position = "y.position",
                             tip.length = 0.01, bracket.size = 0.6,
                             size = 5.5, hide.ns = FALSE) +
  labs(title = "Monocytes: IL6 expression (IL6+ cells only)",
       x = NULL, y = "Expression (normalized)") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y  = element_text(size = 14),
        plot.title   = element_text(face = "bold", hjust = 0.5, size = 16))

ggsave("Monocytes_IL6pos_IL6expr_4groups_violin.png", pA,
       width = 7.5, height = 6.2, dpi = 300)
message("Saved: Monocytes_IL6pos_IL6expr_4groups_violin.png")

# IL6+ fraction by patient (paired, Poor and Good separately)
patient_col <- get_patient_col(mono_int)
df_mono <- FetchData(mono_int, vars = c("Prognosis_Group",
                                        "predicted.celltype.major",
                                        patient_col)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    group       = factor(trimws(.data[["Prognosis_Group"]]), levels = GROUP_LEVELS),
    celltype    = trimws(.data[["predicted.celltype.major"]]),
    patient_raw = trimws(.data[[patient_col]]),
    IL6_counts  = il6_counts_int[cell],
    IL6_pos     = is.finite(IL6_counts) & IL6_counts > 0,
    Prognosis   = ifelse(grepl("^Poor", as.character(group)), "Poor", "Good"),
    Timepoint   = factor(ifelse(grepl("_d7$", as.character(group)),
                                "d7", "d1"), levels = c("d1","d7")),
    patient_base = gsub("(-d1|-d7|_d1|_d7|d1$|d7$)", "", patient_raw)
  ) %>%
  dplyr::filter(!is.na(group), !is.na(patient_raw),
                celltype %in% c("CD14 Mono","CD16 Mono"))

df_pat_il6 <- df_mono %>%
  dplyr::group_by(Prognosis, patient_base, Timepoint) %>%
  dplyr::summarise(n_mono = dplyr::n(), n_IL6pos = sum(IL6_pos, na.rm = TRUE),
                   pct_IL6pos = 100 * n_IL6pos / n_mono, .groups = "drop") %>%
  dplyr::filter(is.finite(pct_IL6pos), n_mono > 0)

plot_il6_fraction <- function(prog, out_png) {
  df_sub <- df_pat_il6 %>% dplyr::filter(Prognosis == prog)
  keep_pat <- df_sub %>%
    dplyr::count(patient_base, Timepoint) %>%
    tidyr::pivot_wider(names_from = Timepoint, values_from = n, values_fill = 0) %>%
    dplyr::filter(d1 > 0 & d7 > 0) %>% dplyr::pull(patient_base)
  df_sub2 <- df_sub %>% dplyr::filter(patient_base %in% keep_pat)
  
  p <- ggplot(df_sub2, aes(x = Timepoint, y = pct_IL6pos)) +
    geom_line(aes(group = patient_base), linewidth = 0.4, alpha = 0.35) +
    geom_boxplot(width = 0.45, fill = "white", color = "black",
                 outlier.shape = NA, linewidth = 0.35) +
    geom_point(size = 2.2, alpha = 0.9) +
    labs(title = paste0(prog, " Prognosis: IL6+ Monocyte Fraction"),
         x = NULL, y = "IL6+ Monocytes (%)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(face = "bold", size = 13),
          axis.text.y = element_text(face = "bold", size = 10),
          plot.title  = element_text(hjust = 0.5, face = "bold", size = 16),
          panel.grid  = element_blank())
  
  if (length(keep_pat) >= 1 && nrow(df_sub2) >= 2) {
    stat_frac <- df_sub2 %>%
      rstatix::wilcox_test(pct_IL6pos ~ Timepoint, paired = TRUE) %>%
      rstatix::add_significance("p") %>%
      dplyr::mutate(group1 = "d1", group2 = "d7",
                    y.position = max(df_sub2$pct_IL6pos, na.rm = TRUE) +
                      max(0.10 * diff(range(df_sub2$pct_IL6pos, na.rm = TRUE)), 1))
    p <- p + ggpubr::stat_pvalue_manual(stat_frac, label = "p.signif",
                                        xmin = "group1", xmax = "group2",
                                        y.position = "y.position",
                                        tip.length = 0.01, bracket.size = 0.35,
                                        size = 5, hide.ns = FALSE)
  }
  ggsave(out_png, p, width = 5, height = 8, dpi = 300)
  message("Saved: ", out_png)
}

plot_il6_fraction("Poor", "Monocytes_IL6pos_fraction_Poor_d1_vs_d7.png")
plot_il6_fraction("Good", "Monocytes_IL6pos_fraction_Good_d1_vs_d7.png")

# =============================================================================
# SECTION 10: Cluster 7 Violin Plots (IL6, TNF, Glycolysis)
# =============================================================================
message("\n>>> Cluster 7 violin plots (IL6 / TNF / Glycolysis)")

COMP_C7 <- list(c("Poor_d1","Poor_d7"),
                c("Good_d1","Good_d7"),
                c("Poor_d7","Good_d7"))

make_cluster_violin(mono_int, cluster_id = 7,
                    score_col = "IL6_JAK_STAT1",
                    module_label = "IL6-JAK-STAT3 signaling",
                    comparisons = COMP_C7,
                    out_png = "Cluster7_IL6module_zscore_violin.png")

make_cluster_violin(mono_int, cluster_id = 7,
                    score_col = "TNFA_NFKB1",
                    module_label = "TNF\u03b1-NF-\u03baB signaling",
                    comparisons = COMP_C7,
                    out_png = "Cluster7_TNFmodule_zscore_violin.png")

make_cluster_violin(mono_int, cluster_id = 7,
                    score_col = "Glycolysis_KEGG_hsa00010",
                    module_label = "Glycolysis activity",
                    comparisons = COMP_C7,
                    out_png = "Cluster7_Glycolysis_zscore_violin.png")

# Hypoxia – 3 pathways combined
COMP_HYP <- list(c("Poor_d1","Poor_d7"),
                 c("Good_d1","Good_d7"),
                 c("Poor_d1","Good_d1"))

hypoxia_pathways <- list(
  HIF1A_Signaling  = "GOBP_HYPOXIA_INDUCIBLE_FACTOR_1ALPHA_SIGNALING_PATHWAY",
  Cellular_Hypoxia = "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
  Myeloid_Hypoxia  = "GSE22282_HYPOXIA_VS_NORMOXIA_MYELOID_DC_UP"
)

all_msig <- msigdbr(species = "Homo sapiens")
obj_c7   <- subset(mono_int, subset = seurat_clusters == "7")
DefaultAssay(obj_c7) <- "RNA"
suppressWarnings({ obj_c7 <- JoinLayers(obj_c7, assay = "RNA") })

p_hyp_list <- lapply(names(hypoxia_pathways), function(nm) {
  gs <- all_msig %>% dplyr::filter(gs_name == hypoxia_pathways[[nm]]) %>%
    dplyr::pull(gene_symbol) %>% unique()
  obj_c7 <- ensure_module_score(obj_c7, nm, gs, paste0(nm, "_tmp"))
  make_cluster_violin(obj_c7, cluster_id = 7, score_col = nm,
                      module_label = nm, comparisons = COMP_HYP,
                      out_png = paste0("C7_", nm, "_violin.png"))
})

p_hyp_combined <- (p_hyp_list[[1]] | p_hyp_list[[2]] | p_hyp_list[[3]]) +
  plot_annotation(title = "Monocyte Cluster 7: Hypoxia Pathway Activity",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          face = "bold",
                                                          size = 16)))
ggsave("C7_Hypoxia_3pathways_combined.png", p_hyp_combined,
       width = 16, height = 5, dpi = 300, bg = "white")
message("Saved: C7_Hypoxia_3pathways_combined.png")

# =============================================================================
# SECTION 11: Inflammatory Gene Heatmap Across Clusters
# =============================================================================
message("\n>>> Inflammatory gene heatmap across monocyte clusters")

Idents(mono_int) <- "seurat_clusters"

genes_heatmap <- unique(c(
  # Chemokines
  "CCL2","CCL3","CCL4","CCL7","CXCL1","CXCL2","CXCL8",
  # Cytokines & inflammation
  "TNF","IL6","IL1B","PTGS2","HMOX1","CYBB","IRF1","THBD"
))
genes_heatmap <- intersect(genes_heatmap, rownames(mono_int))

avg_exp <- AverageExpression(mono_int, features = genes_heatmap,
                             group.by = "seurat_clusters",
                             assays = "RNA", layer = "data")$RNA
colnames(avg_exp) <- sub("^g", "", colnames(avg_exp))
avg_z   <- t(scale(t(avg_exp)))

df_heat <- avg_z %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "cluster", values_to = "zscore") %>%
  dplyr::mutate(
    cluster = factor(as.integer(cluster),
                     levels = as.character(sort(unique(as.integer(cluster))))),
    gene    = factor(gene, levels = rev(genes_heatmap))
  )

z_lim <- min(max(abs(df_heat$zscore), na.rm = TRUE), 2.5)

p_heat <- ggplot(df_heat, aes(x = cluster, y = gene, fill = zscore)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#08306B", mid = "white", high = "#99000D",
                       midpoint = 0, limits = c(-z_lim, z_lim),
                       oob = squish, name = "z-score") +
  labs(title = "Monocytes: Inflammatory Gene Expression",
       x = "Cluster", y = NULL) +
  theme_bw(base_size = 14) +
  theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 14.5),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold.italic", size = 11),
        panel.grid  = element_blank(),
        legend.position = "right")

ggsave("Mono_cluster_inflammatory_heatmap.png", p_heat,
       width = 10, height = 6, dpi = 300, bg = "white")
message("Saved: Mono_cluster_inflammatory_heatmap.png")

# =============================================================================
# SECTION 12: Cluster 7 Volcano – Poor_d7 vs All Other Groups
# =============================================================================
message("\n>>> Cluster 7 volcano: Poor_d7 vs all")

obj_c7_deg <- subset(mono_int, subset = seurat_clusters == "7")
DefaultAssay(obj_c7_deg) <- "RNA"
suppressWarnings({ obj_c7_deg <- JoinLayers(obj_c7_deg, assay = "RNA") })
Idents(obj_c7_deg) <- "Prognosis_Group"

deg_c7_volcano <- FindMarkers(
  obj_c7_deg, ident.1 = "Poor_d7",
  ident.2 = c("Poor_d1","Good_d1","Good_d7"),
  test.use = "wilcox", logfc.threshold = 0.1,
  min.pct = 0.1, assay = "RNA"
) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(
    log10p = -log10(p_val_adj + 1e-300),
    sig    = dplyr::case_when(
      p_val_adj < 0.05 & avg_log2FC >  0.25 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down",
      TRUE ~ "NS"
    ),
    sig = factor(sig, levels = c("Up","Down","NS"))
  )

cat("[C7 volcano] Up:", sum(deg_c7_volcano$sig == "Up"),
    "| Down:", sum(deg_c7_volcano$sig == "Down"), "\n")

# Pathway gene categories for labeling
cat_genes <- list(
  "Glycolysis"    = c("PKM","ENO1","LDHA","PGK1","PGD","ALDOA","GAPDH",
                      "TPI1","HK1","HK2","GPI","PFKL","PFKP","PFKFB3","SLC2A1"),
  "TNF Signaling" = c("TNF","TNFRSF1A","TNFRSF1B","TRAF2","RIPK1","IKBKB",
                      "NFKB1","RELA","NFKBIA","TNFAIP3","TNFAIP6","BIRC3",
                      "ICAM1","CXCL2","CCL2","CCL4"),
  "Inflammation"  = c("IL1B","IL1A","IL6","IL18","CXCL8","CXCL10","PTGS2",
                      "S100A8","S100A9","MIF","NLRP3","PYCARD","CASP1","GSDMD"),
  "IL6/STAT3"     = c("IL6","IL6R","IL6ST","JAK1","JAK2","STAT3",
                      "SOCS1","SOCS3","PIM1","MYC","MCL1")
)
cat_colors <- c("Glycolysis" = "#8B0000", "TNF Signaling" = "#E64B35",
                "Inflammation" = "#FF8C00", "IL6/STAT3" = "darkgreen")

up_genes_c7 <- deg_c7_volcano %>% dplyr::filter(sig == "Up") %>%
  dplyr::pull(gene)

label_cat <- do.call(rbind, lapply(names(cat_genes), function(cat) {
  hits <- intersect(cat_genes[[cat]], up_genes_c7)
  if (length(hits) == 0) return(NULL)
  data.frame(gene = hits, category = cat, stringsAsFactors = FALSE)
})) %>% dplyr::distinct(gene, .keep_all = TRUE)

deg_plot_c7 <- deg_c7_volcano %>%
  dplyr::left_join(label_cat, by = "gene")

down_label_c7 <- deg_c7_volcano %>%
  dplyr::filter(sig == "Down") %>%
  dplyr::arrange(p_val_adj, avg_log2FC) %>% head(7)

x_max_c7 <- max(abs(deg_c7_volcano$avg_log2FC), na.rm = TRUE) * 1.05

p_c7_volcano <- ggplot() +
  geom_point(data = deg_plot_c7 %>% dplyr::filter(sig == "NS"),
             aes(avg_log2FC, log10p), color = "grey92", size = 1.0, alpha = 0.5) +
  geom_point(data = deg_plot_c7 %>% dplyr::filter(sig == "Down"),
             aes(avg_log2FC, log10p), color = "grey60", size = 1.2, alpha = 0.7) +
  geom_point(data = deg_plot_c7 %>%
               dplyr::filter(sig == "Up", is.na(category)),
             aes(avg_log2FC, log10p), color = "grey60", size = 1.2, alpha = 0.7) +
  geom_point(data = deg_plot_c7 %>% dplyr::filter(!is.na(category)),
             aes(avg_log2FC, log10p, color = category), size = 1.8, alpha = 0.95) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed",
             color = "grey80", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey80", linewidth = 0.4) +
  ggrepel::geom_text_repel(data = down_label_c7,
                           aes(avg_log2FC, log10p, label = gene),
                           color = "grey50", size = 2.8, fontface = "bold",
                           box.padding = 0.3, max.overlaps = 20,
                           segment.color = "grey60", segment.size = 0.3,
                           show.legend = FALSE) +
  ggrepel::geom_text_repel(data = deg_plot_c7 %>% dplyr::filter(!is.na(category)),
                           aes(avg_log2FC, log10p, label = gene, color = category),
                           size = 3.5, fontface = "bold", box.padding = 0.4,
                           max.overlaps = 60, segment.size = 0.3,
                           show.legend = FALSE) +
  scale_color_manual(values = cat_colors, name = "Pathway") +
  scale_x_continuous(limits = c(-x_max_c7, x_max_c7)) +
  annotate("text", x =  x_max_c7 * 0.95,
           y = max(deg_c7_volcano$log10p) * 0.97,
           label = "Up in Poor_d7", color = "#E64B35",
           size = 4, fontface = "bold", hjust = 1) +
  annotate("text", x = -x_max_c7 * 0.95,
           y = max(deg_c7_volcano$log10p) * 0.97,
           label = "Down in Poor_d7", color = "#3182BD",
           size = 4, fontface = "bold", hjust = 0) +
  labs(title = "Monocytes (C7): Poor_d7 vs All Other Groups",
       x = expression(Log[2] ~ FC), y = expression(-log[10] ~ P)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 11))

ggsave("C7_volcano_PoorD7_vs_others.png", p_c7_volcano,
       width = 8, height = 6.5, dpi = 300, bg = "white")
message("Saved: C7_volcano_PoorD7_vs_others.png")

# =============================================================================
# SECTION 13: CD14 C4 vs Mono C7 Cluster Overlap UMAP
# =============================================================================
message("\n>>> CD14 C4 vs Mono C7 cluster overlap UMAP")

cd14_c4_cells <- rownames(cd14_int@meta.data)[
  as.character(cd14_int@meta.data$seurat_clusters) == "4"
]
mono_c7_cells <- rownames(mono_int@meta.data)[
  as.character(mono_int@meta.data$seurat_clusters) == "7"
]
overlap_cells <- intersect(cd14_c4_cells, mono_c7_cells)

cat("CD14 C4 cells:", length(cd14_c4_cells),
    "| Mono C7 cells:", length(mono_c7_cells),
    "| Overlap:", length(overlap_cells), "\n")

md_ov <- build_umap_md(mono_int) %>%
  dplyr::mutate(
    cell_group = dplyr::case_when(
      cell %in% overlap_cells                             ~ "Both (CD14-C4 & Mono-C7)",
      cell %in% cd14_c4_cells & !cell %in% mono_c7_cells ~ "CD14 Cluster 4 only",
      cell %in% mono_c7_cells & !cell %in% cd14_c4_cells ~ "Mono Cluster 7 only",
      TRUE                                                ~ "Other"
    ),
    cell_group = factor(cell_group,
                        levels = c("Other","CD14 Cluster 4 only",
                                   "Mono Cluster 7 only",
                                   "Both (CD14-C4 & Mono-C7)"))
  )

color_ov <- c("Other" = "grey85", "CD14 Cluster 4 only" = "#E64B35",
              "Mono Cluster 7 only" = "#4DBBD5",
              "Both (CD14-C4 & Mono-C7)" = "#7B2D8B")
size_ov  <- c("Other" = 0.3, "CD14 Cluster 4 only" = 1.5,
              "Mono Cluster 7 only" = 1.5, "Both (CD14-C4 & Mono-C7)" = 2.5)
alpha_ov <- c("Other" = 0.3, "CD14 Cluster 4 only" = 0.9,
              "Mono Cluster 7 only" = 0.9, "Both (CD14-C4 & Mono-C7)" = 1.0)

md_ov_ordered <- dplyr::bind_rows(
  md_ov %>% dplyr::filter(cell_group == "Other"),
  md_ov %>% dplyr::filter(cell_group == "CD14 Cluster 4 only"),
  md_ov %>% dplyr::filter(cell_group == "Mono Cluster 7 only"),
  md_ov %>% dplyr::filter(cell_group == "Both (CD14-C4 & Mono-C7)")
)

p_overlap <- ggplot(md_ov_ordered,
                    aes(UMAP_1, UMAP_2, color = cell_group,
                        size = cell_group, alpha = cell_group)) +
  geom_point(stroke = 0) +
  scale_color_manual(values = color_ov, name = NULL) +
  scale_size_manual(values  = size_ov,  name = NULL) +
  scale_alpha_manual(values = alpha_ov, name = NULL) +
  coord_fixed() +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = paste0("CD14 C4 only: ",
                          sum(md_ov$cell_group == "CD14 Cluster 4 only"), "\n",
                          "Mono C7 only: ",
                          sum(md_ov$cell_group == "Mono Cluster 7 only"), "\n",
                          "Both:         ",
                          sum(md_ov$cell_group == "Both (CD14-C4 & Mono-C7)")),
           size = 3.5, color = "grey30") +
  labs(title = "CD14 Cluster 4 vs Monocyte Cluster 7") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text  = element_blank(), axis.ticks = element_blank(),
        axis.line  = element_blank(), axis.title = element_blank(),
        legend.position = "right", legend.text = element_text(size = 11),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)),
         size = "none", alpha = "none")

ggsave("CD14_C4_vs_Mono_C7_overlap_UMAP.png", p_overlap,
       width = 8, height = 7, dpi = 300, bg = "white")
message("Saved: CD14_C4_vs_Mono_C7_overlap_UMAP.png")

message("\n>>> All monocyte integration analyses complete.")