# =============================================================================
# 04_CD14_Mono_Analysis.R
# CD14 Monocyte Subcluster Analysis:
#   Cell proportion → DEG volcano → GO BP dotplot → GSEA KEGG → S100 violin
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(rstatix)
library(ggpubr)

# -----------------------------------------------------------------------------
# 0. Shared settings
# -----------------------------------------------------------------------------
obj <- cd14_int
DefaultAssay(obj) <- "RNA"
suppressWarnings({ obj <- JoinLayers(obj, assay = "RNA") })

GROUP_LEVELS <- c("Poor_d1", "Poor_d7", "Good_d1", "Good_d7")
PADJ_CUT     <- 0.05

my_pub_colors <- c(
  "Good_d1" = "#80B1D3",
  "Good_d7" = "#483D8B",
  "Poor_d1" = "#FDB462",
  "Poor_d7" = "#D95F02"
)

comp_list <- list(
  c("Poor_d1", "Poor_d7"), c("Poor_d1", "Good_d1"),
  c("Poor_d1", "Good_d7"), c("Poor_d7", "Good_d1"),
  c("Poor_d7", "Good_d7"), c("Good_d1", "Good_d7")
)

# =============================================================================
# SECTION 1: Helper Functions
# =============================================================================

# Cluster color palette (auto-scaled to number of clusters)
make_cluster_colors <- function(cluster_levels) {
  n <- length(cluster_levels)
  cols <- if (n <= 12) {
    RColorBrewer::brewer.pal(n, "Set3")
  } else {
    grDevices::colorRampPalette(
      c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
        "#FF7F00", "#A65628", "#F781BF", "#999999",
        "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
    )(n)
  }
  setNames(cols, cluster_levels)
}

# Comparison key for pairwise ordering
make_key <- function(a, b, lvls) {
  ia <- match(a, lvls); ib <- match(b, lvls)
  paste0(pmin(ia, ib), "_", pmax(ia, ib))
}

# Paired t-test on d1 vs d7 within a group
run_paired_ttest <- function(df) {
  d1 <- df %>% filter(Group == "d1") %>% arrange(PatientID)
  d7 <- df %>% filter(Group == "d7") %>% arrange(PatientID)
  if (nrow(d1) > 1 && nrow(d1) == nrow(d7))
    t.test(d1$percent, d7$percent, paired = TRUE)$p.value
  else NA
}

# Pairwise Wilcoxon + significance bracket data for a single gene
make_stat_df <- function(df_long, gene_name) {
  comp_keys <- vapply(comp_list,
                      function(x) make_key(x[1], x[2], GROUP_LEVELS),
                      character(1))
  
  stat_df <- df_long %>%
    dplyr::group_by(gene) %>%
    rstatix::pairwise_wilcox_test(expr ~ group, p.adjust.method = "BH") %>%
    dplyr::ungroup() %>%
    rstatix::add_significance("p.adj")
  
  stat_df$key        <- mapply(make_key, stat_df$group1, stat_df$group2,
                               MoreArgs = list(lvls = GROUP_LEVELS))
  stat_df$comp_order <- match(stat_df$key, comp_keys)
  stat_df <- stat_df %>%
    dplyr::filter(!is.na(comp_order)) %>%
    dplyr::arrange(gene, comp_order)
  
  range_df <- df_long %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(ymax = max(expr, na.rm = TRUE),
                     rng  = pmax(max(expr, na.rm = TRUE) - min(expr, na.rm = TRUE), 0.1),
                     .groups = "drop") %>%
    dplyr::mutate(step = pmax(0.15 * rng, 0.15),
                  base = ymax + step)
  
  stat_df <- dplyr::left_join(stat_df, range_df, by = "gene") %>%
    dplyr::mutate(y.position = base + (comp_order - 1) * step)
  
  stat_df
}

# Violin + boxplot with pairwise Wilcoxon brackets
plot_gene_violin <- function(obj, gene, group_col = "Prognosis_Group",
                             title_txt, out_png) {
  df <- Seurat::FetchData(obj, vars = c(group_col, gene)) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")
  
  df$group <- factor(trimws(as.character(df[[group_col]])),
                     levels = GROUP_LEVELS)
  
  df_long <- tidyr::pivot_longer(df, cols = dplyr::all_of(gene),
                                 names_to = "gene", values_to = "expr") %>%
    dplyr::mutate(gene = factor(gene, levels = gene))
  
  stat_df <- make_stat_df(df_long, gene)
  
  p <- ggplot2::ggplot(df_long,
                       ggplot2::aes(x = group, y = expr, fill = group)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.3) +
    ggplot2::geom_boxplot(width = 0.18, outlier.size = 0.15, alpha = 1,
                          linewidth = 0.25, fill = "white", color = "black") +
    ggplot2::scale_fill_manual(values = my_pub_colors, drop = FALSE) +
    ggpubr::stat_pvalue_manual(
      stat_df,
      label        = "p.adj.signif",
      xmin         = "group1",
      xmax         = "group2",
      y.position   = "y.position",
      tip.length   = 0.01,
      bracket.size = 0.35,
      size         = 4,
      hide.ns      = FALSE
    ) +
    ggplot2::labs(title = title_txt, x = NULL,
                  y = "Expression (normalized)") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1,
                                 face = "bold", size = 13),
      axis.text.y = element_text(face = "bold", size = 10),
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 17),
      panel.grid  = element_blank()
    )
  
  ggplot2::ggsave(out_png, p, width = 6, height = 7, dpi = 300)
  message("Saved: ", out_png)
}

# GSEA KEGG dotplot
run_gsea_kegg_dotplot <- function(deg, prefix,
                                  title_txt,
                                  n_up = 10, n_dn = 5,
                                  xlim_range = c(-3.5, 3.5),
                                  exclude_terms = c("Parkinson disease",
                                                    "Huntington disease",
                                                    "Prion disease",
                                                    "Alzheimer disease",
                                                    "Amyotrophic lateral sclerosis")) {
  out_file <- paste0(prefix, "_GSEA_KEGG.png")
  
  ranked <- deg %>%
    dplyr::filter(!is.na(avg_log2FC)) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
  
  gene_list <- setNames(ranked$avg_log2FC, ranked$gene)
  
  gene_df <- bitr(names(gene_list), fromType = "SYMBOL",
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  gene_list_entrez <- gene_list[names(gene_list) %in% gene_df$SYMBOL]
  names(gene_list_entrez) <- gene_df$ENTREZID[
    match(names(gene_list_entrez), gene_df$SYMBOL)]
  gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
  
  gsea_res <- gseKEGG(geneList      = gene_list_entrez,
                      organism      = "hsa",
                      minGSSize     = 10,
                      maxGSSize     = 500,
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      verbose       = FALSE)
  gsea_res <- setReadable(gsea_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  cat("KEGG significant pathways:", nrow(as.data.frame(gsea_res)), "\n")
  
  dir_labels <- c(up = paste0("UP in ",   prefix),
                  dn = paste0("DOWN in ", prefix))
  
  df_kegg <- as.data.frame(gsea_res) %>%
    dplyr::filter(!Description %in% exclude_terms)
  
  df_up <- df_kegg %>% dplyr::filter(NES > 0, p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust) %>% head(n_up) %>%
    dplyr::mutate(direction = dir_labels["up"])
  
  df_dn <- df_kegg %>% dplyr::filter(NES < 0, p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust) %>% head(n_dn) %>%
    dplyr::mutate(direction = dir_labels["dn"])
  
  df_plot <- dplyr::bind_rows(df_up, df_dn) %>%
    dplyr::mutate(
      Description      = stringr::str_wrap(Description, width = 40),
      neg_log10_padj   = -log10(p.adjust),
      direction        = factor(direction, levels = rev(unname(dir_labels)))
    )
  
  term_order <- df_plot %>% dplyr::arrange(NES) %>% dplyr::pull(Description)
  df_plot$Description <- factor(df_plot$Description, levels = unique(term_order))
  
  p <- ggplot(df_plot, aes(x = NES, y = Description)) +
    geom_point(aes(size = setSize, colour = neg_log10_padj)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    scale_colour_gradient(
      low    = "blue", high = "red",
      name   = expression(-log[10]~"(p.adjust)"),
      breaks = scales::pretty_breaks(n = 4),
      labels = scales::number_format(accuracy = 0.1)
    ) +
    scale_size_continuous(range  = c(3, 9), name = "Gene set size",
                          breaks = c(50, 100, 200, 300, 400)) +
    scale_x_continuous(limits = xlim_range,
                       breaks = seq(xlim_range[1], xlim_range[2], 1),
                       expand = expansion(mult = c(0.05, 0.05))) +
    coord_cartesian(xlim = xlim_range) +
    facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL,
         title = title_txt) +
    theme_bw(base_size = 18) +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 17),
      axis.text.y        = element_text(size = 13, face = "bold"),
      axis.text.x        = element_text(size = 13),
      strip.text         = element_text(size = 13, face = "bold"),
      strip.background   = element_rect(fill = "grey92", color = "grey70"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title       = element_text(size = 12),
      legend.text        = element_text(size = 11),
      legend.key.height  = unit(0.8, "cm")
    )
  
  png(out_file, width = 3000, height = 2300, res = 300)
  print(p); dev.off()
  message("Saved: ", out_file)
  invisible(gsea_res)
}

# =============================================================================
# SECTION 2: Cell Proportion – Stacked Bar (Prognosis_Group x Cluster)
# =============================================================================
message("\n>>> Cell proportion stacked bar plot")

cluster_levels <- levels(cd14_int$seurat_clusters)
cluster_colors <- make_cluster_colors(cluster_levels)

cluster_prop <- cd14_int@meta.data %>%
  as.data.frame() %>%
  dplyr::group_by(Prognosis_Group, seurat_clusters) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Prognosis_Group) %>%
  dplyr::mutate(percent = n_cells / sum(n_cells) * 100) %>%
  dplyr::ungroup()

p_bar <- ggplot(cluster_prop,
                aes(x = Prognosis_Group, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = cluster_colors) +
  labs(x = "Prognosis Group", y = "Percentage (%)", fill = "Seurat Cluster") +
  theme_classic(base_size = 14) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.text     = element_text(size = 10),
        legend.title    = element_text(size = 11, face = "bold"))

png("CD14Mono_StackBar_Cluster.png", width = 4000, height = 3000, res = 400)
print(p_bar); dev.off()
message("Saved: CD14Mono_StackBar_Cluster.png")

# =============================================================================
# SECTION 3: Paired Boxplot – All Clusters (d1 vs d7, by Prognosis)
# =============================================================================
message("\n>>> Paired boxplot across all clusters")

meta <- cd14_int@meta.data %>%
  as.data.frame() %>%
  dplyr::mutate(PatientID = stringr::str_remove(orig.ident, "-d[17]$"))

cluster_prop_sample <- meta %>%
  dplyr::group_by(Prognosis, Group, PatientID, seurat_clusters) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Prognosis, Group, PatientID) %>%
  dplyr::mutate(percent = n_cells / sum(n_cells) * 100) %>%
  dplyr::ungroup()

# Keep only patients with both d1 and d7
paired_ids <- cluster_prop_sample %>%
  dplyr::group_by(Prognosis, PatientID) %>%
  dplyr::summarise(has_d1 = any(Group == "d1"),
                   has_d7 = any(Group == "d7"),
                   .groups = "drop") %>%
  dplyr::filter(has_d1 & has_d7) %>%
  dplyr::pull(PatientID)

cat("Paired patient IDs:\n"); print(paired_ids)

# Add "C" prefix to cluster labels
add_c_prefix <- function(df, col = "seurat_clusters", ref_df = cluster_prop_sample) {
  df %>% dplyr::mutate(
    !!col := factor(paste0("C", .data[[col]]),
                    levels = paste0("C", levels(ref_df[[col]])))
  )
}

paired_df <- cluster_prop_sample %>%
  dplyr::filter(PatientID %in% paired_ids) %>%
  dplyr::mutate(Prognosis_Group = factor(
    paste0(Prognosis, "_", Group), levels = GROUP_LEVELS
  )) %>%
  add_c_prefix()

cluster_prop_sample_c <- cluster_prop_sample %>% add_c_prefix()

# Paired t-test results
ttest_results <- paired_df %>%
  dplyr::group_by(seurat_clusters, Prognosis) %>%
  dplyr::do(tibble::tibble(p_value = run_paired_ttest(.))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    sig_label = dplyr::case_when(
      is.na(p_value)    ~ "ns",
      p_value <= 0.001  ~ "***",
      p_value <= 0.01   ~ "**",
      p_value <= 0.05   ~ "*",
      TRUE              ~ "ns"
    )
  )

ypos_tbl <- cluster_prop_sample_c %>%
  dplyr::group_by(seurat_clusters, Prognosis) %>%
  dplyr::summarise(y_max = max(percent, na.rm = TRUE), .groups = "drop")

sig_labels <- ttest_results %>%
  dplyr::filter(sig_label != "ns") %>%
  dplyr::left_join(ypos_tbl, by = c("seurat_clusters", "Prognosis")) %>%
  dplyr::mutate(
    x_left  = ifelse(Prognosis == "Poor", 1, 3),
    x_right = ifelse(Prognosis == "Poor", 2, 4),
    x_mid   = ifelse(Prognosis == "Poor", 1.5, 3.5),
    y_bar   = y_max * 1.08,
    y_label = y_max * 1.14
  )

p_paired <- ggplot(paired_df,
                   aes(x = Prognosis_Group, y = percent, fill = Prognosis_Group)) +
  geom_boxplot(width = 0.7, alpha = 0.85,
               outlier.shape = NA, color = "gray30") +
  geom_point(aes(group = PatientID), size = 1.5, alpha = 0.6,
             position = position_jitter(width = 0.1)) +
  geom_line(aes(group = interaction(PatientID, Prognosis)),
            color = "gray60", alpha = 0.4) +
  # Significance bracket (horizontal bar)
  geom_segment(data = sig_labels,
               aes(x = x_left, xend = x_right, y = y_bar, yend = y_bar),
               color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  # Bracket ticks
  geom_segment(data = sig_labels,
               aes(x = x_left,  xend = x_left,
                   y = y_bar * 0.995, yend = y_bar),
               color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_segment(data = sig_labels,
               aes(x = x_right, xend = x_right,
                   y = y_bar * 0.995, yend = y_bar),
               color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = sig_labels,
            aes(x = x_mid, y = y_label, label = sig_label),
            size = 4.5, color = "black", inherit.aes = FALSE) +
  facet_wrap(~ seurat_clusters, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = my_pub_colors, drop = FALSE) +
  labs(title = "", x = NULL, y = "Cell proportion (%)") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x          = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.background     = element_rect(fill = "gray95", color = NA),
    strip.text           = element_text(size = 10, face = "bold"),
    legend.position      = c(0.98, 0.005),
    legend.justification = c(1, 0),
    legend.background    = element_rect(fill = "white", color = "gray80",
                                        linewidth = 0.3),
    legend.margin        = margin(4, 6, 4, 6),
    legend.title         = element_text(size = 10, face = "bold"),
    legend.text          = element_text(size = 9)
  )

png("CD14Mono_AllClusters_Paired_BoxPlot.png",
    width = 4000, height = 2800, res = 400)
print(p_paired); dev.off()
message("Saved: CD14Mono_AllClusters_Paired_BoxPlot.png")

# =============================================================================
# SECTION 4: Cluster 4 DEG – Poor_d7 vs All Other Cluster-Groups
# =============================================================================
message("\n>>> DEG: Cluster 4 (Poor_d7) vs all other cluster-groups")

cd14_int$cluster_group <- paste0(cd14_int$seurat_clusters, "_",
                                 cd14_int$Prognosis_Group)
DefaultAssay(cd14_int) <- "RNA"
suppressWarnings({ cd14_int <- JoinLayers(cd14_int, assay = "RNA") })
Idents(cd14_int) <- "cluster_group"

target_group <- "4_Poor_d7"
other_groups <- setdiff(levels(Idents(cd14_int)), target_group)

deg_c4_pg <- FindMarkers(
  cd14_int,
  ident.1         = target_group,
  ident.2         = other_groups,
  min.pct         = 0.1,
  logfc.threshold = 0,
  test.use        = "wilcox"
) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(dplyr::desc(avg_log2FC))

write.csv(deg_c4_pg, "Cluster4_Poor_d7_vs_all_DEG.csv", row.names = FALSE)
message("Saved: Cluster4_Poor_d7_vs_all_DEG.csv")

# =============================================================================
# SECTION 5: Volcano Plots – Cluster 4 (two gene-set themes)
# =============================================================================

build_volcano <- function(deg_df,
                          goi_up, goi_down,
                          out_png,
                          title_txt,
                          padj_cut = PADJ_CUT,
                          fc_cut   = 1.0,
                          n_up     = 20,
                          n_dn     = 25) {
  
  deg_plot <- deg_df %>%
    dplyr::mutate(p_val_adj_plot = ifelse(p_val_adj == 0, 1e-300, p_val_adj))
  rownames(deg_plot) <- deg_plot$gene
  
  exclude_pat <- "^AC|^MT|^G0"
  
  label_up <- deg_plot %>%
    dplyr::filter(p_val_adj < padj_cut, avg_log2FC >= fc_cut,
                  !grepl(exclude_pat, gene)) %>%
    dplyr::arrange(p_val_adj) %>% head(n_up) %>% dplyr::pull(gene)
  
  label_dn <- deg_plot %>%
    dplyr::filter(p_val_adj < padj_cut, avg_log2FC <= -fc_cut,
                  !grepl(exclude_pat, gene)) %>%
    dplyr::arrange(p_val_adj) %>% head(n_dn) %>% dplyr::pull(gene)
  
  label_genes <- unique(c(
    intersect(c(goi_up, goi_down), deg_plot$gene),
    label_up, label_dn
  ))
  
  col_vec <- rep("grey70", nrow(deg_plot))
  names(col_vec) <- rownames(deg_plot)
  col_vec[deg_plot$avg_log2FC >=  fc_cut & deg_plot$p_val_adj < padj_cut] <- "red"
  col_vec[deg_plot$avg_log2FC <= -fc_cut & deg_plot$p_val_adj < padj_cut] <- "blue"
  
  x_max     <- max(abs(deg_plot$avg_log2FC), na.rm = TRUE)
  ymax      <- max(-log10(deg_plot$p_val_adj_plot), na.rm = TRUE)
  y_lim_use <- c(0, ymax + 0.5)
  
  png(out_png, width = 2200, height = 2500, res = 300)
  
  p <- EnhancedVolcano(
    deg_plot,
    lab              = deg_plot$gene,
    selectLab        = label_genes,
    x                = "avg_log2FC",
    y                = "p_val_adj_plot",
    xlab             = bquote(~Log[2]~ "FC"),
    ylab             = expression(-log[10]~"P"["adj"]),
    title            = title_txt,
    subtitle         = "",
    pCutoff          = padj_cut,
    FCcutoff         = fc_cut,
    labSize          = 3.5,
    pointSize        = 1.5,
    colAlpha         = 0.8,
    drawConnectors   = TRUE,
    widthConnectors  = 0.5,
    max.overlaps     = 30,
    legendPosition   = "right",
    legendLabSize    = 6,
    legendIconSize   = 3.0,
    boxedLabels      = TRUE,
    caption          = NULL,
    ylim             = y_lim_use,
    xlim             = c(-x_max, x_max),
    col              = rep("grey70", 4),
    colCustom        = col_vec
  ) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    theme(
      legend.position  = "none",
      plot.title       = element_text(size = 17, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Legend annotation
  y_legend <- y_lim_use[2] * 0.98
  p <- p +
    annotate("point", x = -x_max * 0.90, y = y_legend,
             colour = "blue", size = 2.5) +
    annotate("text",  x = -x_max * 0.84, y = y_legend,
             label = paste0("Down in ", target_group),
             hjust = 0, size = 4.5) +
    annotate("point", x =  x_max * 0.40, y = y_legend,
             colour = "red",  size = 2.5) +
    annotate("text",  x =  x_max * 0.46, y = y_legend,
             label = paste0("Up in ",   target_group),
             hjust = 0, size = 4.5)
  
  print(p); dev.off()
  message("Saved: ", out_png)
  invisible(p)
}

# --- Theme 1: Immune / Inflammatory genes ---
goi_up_immune   <- c("ATP5F1A", "IDH2", "SHMT2", "PHB", "PRDX3",
                     "MKI67", "BIRC5", "PCLAF", "CDC123", "KIF22")
goi_down_immune <- c("IFI27", "IFITM1", "IFI6", "IFI44L", "ISG15",
                     "MX1", "MX2", "IFI44", "XAF1",
                     "CXCL8", "CXCL16", "CCL5", "C5AR1", "DUSP2",
                     "GPR183", "SIGLEC1", "CLEC4E", "ITGAX", "VSIG4",
                     "LGALS3BP", "LGALS2", "VEGFA")

build_volcano(
  deg_df   = deg_c4_pg,
  goi_up   = goi_up_immune,
  goi_down = goi_down_immune,
  out_png  = "Cluster4_Poor_d7_vs_all_Volcano_Immune.png",
  title_txt = "CD14 Monocytes (C4): Poor_d7 vs all cluster groups",
  fc_cut   = 1.0, n_up = 20, n_dn = 25
)

# --- Theme 2: Biosynthetic / Metabolic genes ---
goi_up_meta <- c("FBL", "EBNA1BP2", "RUVBL1", "RUVBL2", "SNRNP25",
                 "EIF3L", "EIF2D",
                 "PSMD14", "PRDX4", "CACYBP", "USP5",
                 "GPI", "PGK1", "PGAM1", "LDHA", "LDHB")

build_volcano(
  deg_df    = deg_c4_pg,
  goi_up    = goi_up_meta,
  goi_down  = character(0),
  out_png   = "Cluster4_Poor_d7_vs_all_Volcano_Metabolic.png",
  title_txt = "CD14 Monocytes (C4): Poor_d7 vs all cluster groups",
  fc_cut    = 0.75, n_up = 20, n_dn = 25
)

# =============================================================================
# SECTION 6: GO BP Dotplot – Cluster 4 (Poor_d7 vs all cluster-groups)
# =============================================================================
message("\n>>> GO BP dotplot: Cluster 4 (Poor_d7) vs all")

# Re-run DEG with logfc.threshold = 0.25 for GO input
deg_c4_go <- FindMarkers(
  cd14_int,
  ident.1         = target_group,
  ident.2         = other_groups,
  min.pct         = 0.1,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
) %>%
  tibble::rownames_to_column("gene")

universe_genes <- deg_c4_go$gene

run_go_bp <- function(genes, universe, simplify_cut = 0.7) {
  if (length(genes) == 0) return(NULL)
  ego <- tryCatch(
    enrichGO(gene = genes, universe = universe,
             OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
             ont = "BP", pAdjustMethod = "BH",
             pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
    error = function(e) NULL
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  tryCatch(
    clusterProfiler::simplify(ego, cutoff = simplify_cut,
                              by = "p.adjust", select_fun = min),
    error = function(e) ego
  )
}

get_top_go <- function(ego, n = 7, reg) {
  if (is.null(ego)) return(NULL)
  as.data.frame(ego) %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust) %>%
    head(n) %>%
    dplyr::mutate(regulation = reg)
}

genes_up_go   <- deg_c4_go %>%
  dplyr::filter(p_val_adj < PADJ_CUT, avg_log2FC > 0.75) %>%
  dplyr::pull(gene)
genes_down_go <- deg_c4_go %>%
  dplyr::filter(p_val_adj < PADJ_CUT, avg_log2FC < -0.75) %>%
  dplyr::pull(gene)

ego_up   <- run_go_bp(genes_up_go,   universe_genes)
ego_down <- run_go_bp(genes_down_go, universe_genes)

go_comb <- dplyr::bind_rows(
  get_top_go(ego_up,   n = 7, reg = "up"),
  get_top_go(ego_down, n = 7, reg = "down")
)

if (!is.null(go_comb) && nrow(go_comb) > 0) {
  term_order <- go_comb %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(min_p = min(p.adjust), .groups = "drop") %>%
    dplyr::arrange(min_p) %>%
    dplyr::pull(Description)
  
  go_comb <- go_comb %>%
    dplyr::mutate(
      Term_wrapped = factor(stringr::str_wrap(Description, 40),
                            levels = rev(stringr::str_wrap(term_order, 40))),
      regulation   = factor(regulation, levels = c("down", "up"))
    )
  
  p_go <- ggplot(go_comb, aes(x = regulation, y = Term_wrapped)) +
    geom_point(aes(size = Count, colour = regulation)) +
    scale_colour_manual(values = c(down = "blue", up = "red"), name = "") +
    scale_size_continuous(range = c(2.5, 8), name = "Count") +
    labs(x = NULL, y = NULL,
         title = "GO BP: CD14 Monocytes (C4) in Poor_d7") +
    theme_bw(base_size = 20) +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold"),
      axis.text.y        = element_text(size = 17, face = "bold"),
      axis.text.x        = element_text(size = 17, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
  out_go <- "CD14_C4_Poor_d7_vs_all_GO_BP_dotplot.png"
  png(out_go, width = 3500, height = 2500, res = 300)
  print(p_go); dev.off()
  message("Saved: ", out_go)
}

# =============================================================================
# SECTION 7: GSEA KEGG Dotplots
# =============================================================================
message("\n>>> GSEA KEGG: Cluster 4 (Poor_d7) vs all cluster-groups")
run_gsea_kegg_dotplot(
  deg       = deg_c4_pg,
  prefix    = "CD14_C4_Poor_d7_vs_all",
  title_txt = "GSEA KEGG: CD14 Monocytes (C4) Poor_d7 vs all"
)

message("\n>>> GSEA KEGG: Cluster 4 vs all other clusters (seurat_clusters)")
Idents(cd14_int) <- "seurat_clusters"
deg_c4_clust <- FindMarkers(
  cd14_int,
  ident.1         = "4",
  min.pct         = 0.1,
  logfc.threshold = 0,
  test.use        = "wilcox"
) %>%
  tibble::rownames_to_column("gene")

run_gsea_kegg_dotplot(
  deg       = deg_c4_clust,
  prefix    = "Cluster4_vs_All",
  title_txt = "GSEA KEGG: Cluster 4 vs All Other Clusters"
)

# Restore cluster_group ident for downstream use
Idents(cd14_int) <- "cluster_group"

# =============================================================================
# SECTION 8: S100A8 / S100A9 / S100A12 Violin Plots in CD14 Monocytes
# =============================================================================
message("\n>>> S100 violin plots")

for (gene in c("S100A8", "S100A9", "S100A12")) {
  plot_gene_violin(
    obj       = cd14_int,
    gene      = gene,
    title_txt = paste0("CD14 Monocytes: ", gene),
    out_png   = paste0("CD14_", gene, ".png")
  )
}

message("\n>>> CD14 monocyte analysis complete.")