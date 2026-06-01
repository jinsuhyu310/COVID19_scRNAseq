# =============================================================================
# 07_B_Cell_Analysis.R
# B Cell Analysis:
#   Cell proportion (whole PBMC & within B) → UMAP → DEG Volcano →
#   GO enrichment → TNF/IL6 module score violin → IL6 DotPlot
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(patchwork)
library(ggrepel)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ggpubr)

# -----------------------------------------------------------------------------
# 0. Shared settings
# -----------------------------------------------------------------------------
GROUP_LEVELS  <- c("Poor_d1", "Poor_d7", "Good_d1", "Good_d7")
B_SUB_LEVELS  <- c("naive B", "memory B", "Plasmablast")
PADJ_CUT      <- 0.05
FC_CUT        <- 1.0

my_pub_colors <- c(
  "Good_d1" = "#80B1D3", "Good_d7" = "#483D8B",
  "Poor_d1" = "#FDB462", "Poor_d7" = "#D95F02"
)

my_b_colors <- c(
  "naive B"     = "#B2DF8A",
  "memory B"    = "#33A02C",
  "Plasmablast" = "#E7298A"
)

# Genes of interest for volcano plots
GENES_UP_POOR_D7 <- c(
  "NR4A1","NR4A2","NR4A3",
  "BCL3","NFKBID","MAP3K8","RIPK2","TRAF4","IRF1","IL10RA",
  "GADD45B","DDIT4","PMAIP1","AREG",
  "SLC2A3","SLC7A5","SREBF2","CYP51A1",
  "CDKN1A","DBF4","KDM6B",
  "IL6","IL6R","STAT3","TNF","NFKB1","RELA"
)

GENES_DOWN_POOR_D7 <- c(
  "IFI44L","IFI27","MX1","MX2","IFITM1","IFITM3","ISG15",
  "OAS1","OAS2","IRF7","STAT1",
  "TCL1A","CCR7","SELL","CD38","CD180","BTLA","FCMR"
)

# =============================================================================
# SECTION 1: Helper Functions
# =============================================================================

# Safe layer join for Seurat v5
safe_join_layers <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  suppressWarnings(tryCatch(JoinLayers(obj, assay = assay),
                            error = function(e) obj))
}

# Auto-detect patient ID column
get_patient_col <- function(obj) {
  cands <- c("Patient_ID","patient","Patient","orig.ident","SampleID","donor","Donor")
  col   <- cands[cands %in% colnames(obj@meta.data)][1]
  if (is.na(col)) stop("No patient ID column found in meta.data.")
  col
}

# Paired t-test p-value (d1 vs d7)
paired_ttest_pval <- function(df, g1, g2, val_col = "percent") {
  x <- df[[val_col]][df$Prognosis_Group == g1]
  y <- df[[val_col]][df$Prognosis_Group == g2]
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) >= 2) t.test(y[ok], x[ok], paired = TRUE)$p.value else NA_real_
}

signif_label <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "NA", p < 0.001 ~ "***",
    p < 0.01  ~ "**", p < 0.05  ~ "*", TRUE ~ "ns"
  )
}

# Generic significance bracket (horizontal bar + ticks + label)
add_sig_bracket <- function(p, pvals_df,
                            x_col = "x_mid", y_bar_col = "y_bar",
                            y_label_col = "y_label", label_col = "signif") {
  p +
    geom_segment(data = pvals_df,
                 aes(x = x_left, xend = x_right,
                     y = .data[[y_bar_col]], yend = .data[[y_bar_col]]),
                 inherit.aes = FALSE, color = "black", linewidth = 0.5) +
    geom_segment(data = pvals_df,
                 aes(x = x_left, xend = x_left,
                     y = .data[[y_bar_col]] * 0.995,
                     yend = .data[[y_bar_col]]),
                 inherit.aes = FALSE, color = "black", linewidth = 0.5) +
    geom_segment(data = pvals_df,
                 aes(x = x_right, xend = x_right,
                     y = .data[[y_bar_col]] * 0.995,
                     yend = .data[[y_bar_col]]),
                 inherit.aes = FALSE, color = "black", linewidth = 0.5) +
    geom_text(data = pvals_df,
              aes(x = .data[[x_col]], y = .data[[y_label_col]],
                  label = .data[[label_col]]),
              inherit.aes = FALSE, size = 4.5, color = "black")
}

# Paired boxplot + line for one prognosis group
make_paired_prop_plot <- function(df_sub, groups_use,
                                  fill_col, group_col = "Prognosis_Group",
                                  facet_col, y_lab, pvals_df = NULL) {
  p <- ggplot(df_sub, aes(x = .data[[group_col]], y = percent,
                          fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_line(aes(group = Patient_ID), color = "gray40",
              linewidth = 0.6, alpha = 0.7) +
    geom_point(shape = 21, color = "black", size = 2) +
    facet_wrap(as.formula(paste("~", facet_col)),
               scales = "free_y", nrow = 1) +
    scale_fill_manual(values = fill_col, breaks = groups_use) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(x = NULL, y = y_lab, title = "") +
    theme_classic(base_size = 12) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 10,
                                          color = "black"),
          axis.text.y      = element_text(size = 10, color = "black"),
          axis.title.y     = element_text(size = 11),
          axis.line        = element_line(colour = "black", linewidth = 0.5),
          strip.text       = element_text(face = "bold", size = 12,
                                          color = "black"),
          strip.background = element_rect(color = "black", fill = "grey90",
                                          linewidth = 1),
          legend.position  = "none")
  if (!is.null(pvals_df)) {
    p <- p +
      geom_text(data = pvals_df,
                aes(x = 1.5, y = y.position,
                    label = ifelse(is.na(p.value), "n<2",
                                   paste0("p=", signif(p.value, 3),
                                          " (", signif_label, ")"))),
                inherit.aes = FALSE, size = 4)
  }
  p
}

# GO enrichment with simplify
run_go_enrichment <- function(genes, comp_name, ont = "BP",
                              p_cutoff = 0.05, q_cutoff = 0.2,
                              simplify_cutoff = 0.5,
                              max_count = 150) {
  if (length(genes) < 3) {
    warning("[", comp_name, "] Too few genes (n=", length(genes), "). Skipping.")
    return(NULL)
  }
  ego <- tryCatch(
    enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
             ont = ont, pAdjustMethod = "BH",
             pvalueCutoff = p_cutoff, qvalueCutoff = q_cutoff,
             minGSSize = 5, maxGSSize = 500, readable = TRUE),
    error = function(e) NULL
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  
  ego <- tryCatch(
    clusterProfiler::simplify(ego, cutoff = simplify_cutoff,
                              by = "p.adjust", select_fun = min),
    error = function(e) ego
  )
  df <- as.data.frame(ego) %>%
    dplyr::filter(Count <= max_count) %>%
    dplyr::arrange(p.adjust, Count)
  if (nrow(df) == 0) return(NULL)
  ego@result <- df
  ego
}

# Bidirectional GO barplot (up = right, down = left)
save_go_bidirectional <- function(ego_up, ego_down, comp_name, out_png,
                                  top_n_up = 9, top_n_down = 9,
                                  width = 17, height = 11) {
  if (is.null(ego_up) && is.null(ego_down)) {
    warning("[", comp_name, "] Both ego objects are NULL. Skipping.")
    return(invisible(NULL))
  }
  get_df <- function(ego, n, direction) {
    if (is.null(ego)) return(NULL)
    df <- as.data.frame(ego) %>% dplyr::arrange(p.adjust) %>% head(n)
    if (nrow(df) == 0) return(NULL)
    df$direction <- direction; df
  }
  df_up   <- get_df(ego_up,   top_n_up,   "Up-regulated in Poor_d7")
  df_down <- get_df(ego_down, top_n_down, "Down-regulated in Poor_d7")
  
  plot_df <- dplyr::bind_rows(df_up, df_down) %>%
    dplyr::mutate(
      log10_p = -log10(p.adjust),
      value   = ifelse(direction == "Up-regulated in Poor_d7",
                       log10_p, -log10_p),
      Description_wrapped = stringr::str_wrap(Description, 60)
    )
  
  up_terms   <- plot_df %>% dplyr::filter(direction == "Up-regulated in Poor_d7") %>%
    dplyr::arrange(p.adjust) %>% dplyr::pull(Description_wrapped)
  down_terms <- plot_df %>% dplyr::filter(direction == "Down-regulated in Poor_d7") %>%
    dplyr::arrange(p.adjust) %>% dplyr::pull(Description_wrapped)
  
  plot_df$Description_wrapped <- factor(plot_df$Description_wrapped,
                                        levels = c(up_terms, rev(down_terms)))
  
  x_lim <- max(abs(plot_df$value), na.rm = TRUE) * 1.1
  brks  <- pretty(c(0, x_lim)); brks <- brks[brks > 0]
  
  p <- ggplot(plot_df, aes(x = value, y = Description_wrapped,
                           fill = direction)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = 0, colour = "black") +
    scale_fill_manual(
      values = c("Down-regulated in Poor_d7" = "#4575b4",
                 "Up-regulated in Poor_d7"   = "#b2182b"),
      name = NULL
    ) +
    scale_x_continuous(limits = c(-x_lim, x_lim),
                       breaks = c(-rev(brks), brks),
                       labels = function(x) abs(x),
                       name   = "-log10(adj. p-value)") +
    labs(y = "", title = comp_name) +
    theme_bw() +
    theme(plot.title           = element_text(hjust = 0.5, face = "bold",
                                              size = 23),
          axis.text.x          = element_text(size = 17),
          axis.text.y          = element_text(size = 20, face = "bold"),
          axis.title.x         = element_text(size = 18, face = "bold"),
          legend.position      = c(0.55, 0.9),
          legend.justification = c("left","top"),
          legend.text          = element_text(size = 18, face = "bold"),
          panel.grid           = element_blank())
  
  ggsave(out_png, p, width = width, height = height, dpi = 300)
  message("Saved: ", out_png)
  invisible(p)
}

# Volcano plot (EnhancedVolcano wrapper)
build_volcano <- function(deg_df, title_txt, out_png,
                          genes_up, genes_down,
                          fc_cut = 0.5, padj_cut = 0.05,
                          width = 1800, height = 2200, res = 300) {
  deg_plot <- deg_df %>%
    dplyr::filter(!grepl("^MT|^AC|^XIST", gene)) %>%
    dplyr::mutate(
      p_val_adj_plot = ifelse(p_val_adj == 0, 1e-300, p_val_adj),
      sig_fc         = abs(avg_log2FC) >= fc_cut & p_val_adj < padj_cut,
      color_custom   = dplyr::case_when(
        sig_fc & avg_log2FC > 0 ~ "red",
        sig_fc & avg_log2FC < 0 ~ "royalblue3",
        TRUE                    ~ "grey70"
      )
    )
  
  col_vec        <- deg_plot$color_custom
  names(col_vec) <- deg_plot$gene
  
  sig_genes   <- deg_plot %>% dplyr::filter(sig_fc) %>% dplyr::pull(gene)
  label_genes <- intersect(c(genes_up, genes_down), sig_genes)
  
  x_max     <- max(abs(deg_plot$avg_log2FC), na.rm = TRUE)
  ymax      <- max(-log10(deg_plot$p_val_adj_plot), na.rm = TRUE)
  y_lim_use <- c(0, ymax + 1)
  
  png(out_png, width = width, height = height, res = res)
  p <- EnhancedVolcano(
    deg_plot,
    lab             = deg_plot$gene,
    selectLab       = label_genes,
    x               = "avg_log2FC",
    y               = "p_val_adj_plot",
    xlab            = bquote(~log[2]~"FC"),
    ylab            = bquote(-log[10]~italic(P)),
    title           = title_txt, subtitle = "",
    pCutoff         = padj_cut, FCcutoff = fc_cut,
    labSize         = 3.5, pointSize = 1.5, colAlpha = 0.8,
    drawConnectors  = TRUE, widthConnectors = 0.5,
    max.overlaps    = 20, legendPosition = "none",
    boxedLabels     = TRUE, caption = NULL,
    ylim = y_lim_use, xlim = c(-x_max, x_max),
    col      = rep("grey70", 4),
    colCustom = col_vec
  ) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    theme(plot.title       = element_text(size = 15, face = "bold",
                                          hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  y_legend <- y_lim_use[2] * 0.95
  p <- p +
    annotate("point", x = -x_max * 0.92, y = y_legend,
             colour = "royalblue3", size = 2.5) +
    annotate("text",  x = -x_max * 0.87, y = y_legend,
             label = "Down in Poor_d7", hjust = 0, size = 4) +
    annotate("point", x =  x_max * 0.55, y = y_legend,
             colour = "red", size = 2.5) +
    annotate("text",  x =  x_max * 0.60, y = y_legend,
             label = "Up in Poor_d7", hjust = 0, size = 4)
  
  print(p); dev.off()
  message("Saved: ", out_png)
  invisible(p)
}

# Violin + boxplot with Wilcoxon stats (z-score, paired patient lines)
make_violin_module <- function(obj, score_col, plot_title, out_png,
                               comparisons = list(
                                 c("Poor_d1","Poor_d7"),
                                 c("Good_d1","Good_d7"),
                                 c("Poor_d7","Good_d7")
                               ),
                               y_step = 0.13, width = 6.2, height = 5.2) {
  patient_col  <- get_patient_col(obj)
  fill_map     <- my_pub_colors[GROUP_LEVELS]
  pair_levels  <- sapply(comparisons, function(x) paste(x, collapse = "__"))
  
  df_long <- obj@meta.data %>%
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
      mu           = mean(score, na.rm = TRUE),
      sig          = sd(score, na.rm = TRUE),
      z            = ifelse(is.na(sig) | sig == 0, NA_real_,
                            (score - mu) / sig)
    ) %>%
    dplyr::filter(is.finite(z)) %>%
    dplyr::mutate(z_cap = pmax(-1, pmin(1, z)))
  
  df_pat <- df_long %>%
    dplyr::group_by(Prognosis, Timepoint, patient_base) %>%
    dplyr::summarise(pat_mean = mean(z_cap, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(is.finite(pat_mean)) %>%
    dplyr::mutate(group = factor(paste0(Prognosis, "_", Timepoint),
                                 levels = GROUP_LEVELS))
  
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
      p.adj    = p.adjust(p, method = "BH"),
      p.signif = dplyr::case_when(
        is.na(p.adj)  ~ "NA", p.adj <= 1e-4 ~ "****",
        p.adj <= 1e-3 ~ "***", p.adj <= 1e-2 ~ "**",
        p.adj <= 5e-2 ~ "*",  TRUE           ~ "ns"
      ),
      pair_key = factor(paste(group1, group2, sep = "__"),
                        levels = pair_levels)
    ) %>%
    dplyr::arrange(pair_key) %>%
    dplyr::mutate(y.position = max(df_long$z_cap, na.rm = TRUE) +
                    1e-6 + seq_len(dplyr::n()) * y_step)
  
  ylim_top <- max(stat_tbl$y.position) + 0.05
  
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
    coord_cartesian(ylim = c(-1.0, ylim_top)) +
    labs(x = NULL, y = "Module Score (z-score)", title = plot_title) +
    theme_bw(base_size = 14) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y  = element_text(color = "black"),
          plot.title   = element_text(face = "bold", hjust = 0.5),
          panel.grid   = element_blank(),
          legend.position = "none")
  
  ggsave(out_png, p, width = width, height = height, dpi = 300, bg = "white")
  message("Saved: ", out_png)
  invisible(p)
}

# =============================================================================
# SECTION 2: B Cell Proportion within Whole PBMC
# =============================================================================
message("\n>>> B cell subset proportion (whole PBMC)")

DefaultAssay(merged_int) <- "RNA"
merged_int <- safe_join_layers(merged_int)
merged_int$Prognosis_Group <- factor(merged_int$Prognosis_Group,
                                     levels = GROUP_LEVELS)

md_pbmc <- merged_int@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::transmute(
    cell,
    orig.ident      = as.character(orig.ident),
    Prognosis_Group = as.character(Prognosis_Group),
    fineB = trimws(as.character(predicted.celltype.major_fineB))
  ) %>%
  dplyr::filter(!is.na(orig.ident), Prognosis_Group %in% GROUP_LEVELS) %>%
  dplyr::mutate(
    Prognosis_Group = factor(Prognosis_Group, levels = GROUP_LEVELS),
    Patient_ID      = sub("-d[17]$", "", orig.ident)
  )

denom_pbmc <- md_pbmc %>%
  dplyr::group_by(orig.ident, Prognosis_Group, Patient_ID) %>%
  dplyr::summarise(total_cells = dplyr::n(), .groups = "drop")

num_pbmc <- md_pbmc %>%
  dplyr::filter(fineB %in% B_SUB_LEVELS) %>%
  dplyr::group_by(orig.ident, Prognosis_Group, Patient_ID, fineB) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

prop_pbmc <- dplyr::left_join(num_pbmc, denom_pbmc,
                              by = c("orig.ident","Prognosis_Group","Patient_ID")) %>%
  dplyr::mutate(percent = n_cells / total_cells * 100,
                fineB   = factor(fineB, levels = B_SUB_LEVELS))

make_paired_whole_plot <- function(prop_df, prognosis, out_png) {
  groups_use <- if (prognosis == "Poor") c("Poor_d1","Poor_d7") else
    c("Good_d1","Good_d7")
  g1 <- groups_use[1]; g2 <- groups_use[2]
  
  df_sub <- prop_df %>%
    dplyr::filter(Prognosis_Group %in% groups_use) %>%
    dplyr::mutate(Prognosis_Group = factor(as.character(Prognosis_Group),
                                           levels = groups_use))
  
  wide <- df_sub %>%
    dplyr::select(fineB, Patient_ID, Prognosis_Group, percent) %>%
    tidyr::pivot_wider(names_from = Prognosis_Group, values_from = percent)
  
  pvals <- wide %>%
    dplyr::group_by(fineB) %>%
    dplyr::summarise(
      n_pairs = sum(!is.na(.data[[g1]]) & !is.na(.data[[g2]])),
      p.value = if (sum(!is.na(.data[[g1]]) & !is.na(.data[[g2]])) >= 2)
        t.test(.data[[g2]], .data[[g1]], paired = TRUE)$p.value else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      signif_label = signif_label(p.value),
      y.position   = NA_real_
    )
  
  ypos <- df_sub %>%
    dplyr::group_by(fineB) %>%
    dplyr::summarise(y.position = {
      mx <- suppressWarnings(max(percent, na.rm = TRUE))
      if (is.finite(mx)) mx * 1.10 else 1
    }, .groups = "drop")
  pvals <- dplyr::left_join(pvals, ypos, by = "fineB") %>%
    dplyr::mutate(y.position = y.position.y) %>%
    dplyr::select(-y.position.x, -y.position.y)
  
  p <- ggplot(df_sub, aes(x = Prognosis_Group, y = percent,
                          fill = Prognosis_Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_line(aes(group = Patient_ID), color = "gray40",
              linewidth = 0.6, alpha = 0.7) +
    geom_point(shape = 21, color = "black", size = 2) +
    facet_wrap(~ fineB, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = my_pub_colors, breaks = groups_use) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    geom_text(data = pvals,
              aes(x = 1.5, y = y.position,
                  label = ifelse(is.na(p.value), "n<2",
                                 paste0("p=", signif(p.value, 3),
                                        " (", signif_label, ")"))),
              inherit.aes = FALSE, size = 4) +
    labs(x = NULL, y = "Percentage of cells (%)", title = "") +
    theme_classic(base_size = 12) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y      = element_text(size = 10),
          axis.title.y     = element_text(size = 11),
          strip.text       = element_text(face = "bold", size = 12),
          strip.background = element_rect(color = "black", fill = "grey90",
                                          linewidth = 1),
          legend.position  = "none")
  
  ggsave(out_png, p, width = 7, height = 3.5, dpi = 300)
  message("Saved: ", out_png)
  invisible(p)
}

make_paired_whole_plot(prop_pbmc, "Poor", "B_subset_whole_proportion_Poor.png")
make_paired_whole_plot(prop_pbmc, "Good", "B_subset_whole_proportion_Good.png")

# Memory B proportion: Poor vs Good (stacked, 2 panels)
df_memB <- prop_pbmc %>%
  dplyr::filter(fineB == "memory B") %>%
  dplyr::mutate(Prognosis_Group = factor(as.character(Prognosis_Group),
                                         levels = GROUP_LEVELS))

make_memB_panel <- function(df, groups_use) {
  g1 <- groups_use[1]; g2 <- groups_use[2]
  df_sub <- df %>%
    dplyr::filter(Prognosis_Group %in% groups_use) %>%
    dplyr::mutate(Prognosis_Group = factor(as.character(Prognosis_Group),
                                           levels = groups_use))
  
  wide <- df_sub %>%
    dplyr::select(Patient_ID, Prognosis_Group, percent) %>%
    tidyr::pivot_wider(names_from = Prognosis_Group, values_from = percent)
  n_pairs <- sum(!is.na(wide[[g1]]) & !is.na(wide[[g2]]))
  pval    <- if (n_pairs >= 2) t.test(wide[[g2]], wide[[g1]], paired = TRUE)$p.value else NA_real_
  ymax    <- max(df_sub$percent, na.rm = TRUE)
  
  ggplot(df_sub, aes(x = Prognosis_Group, y = percent, fill = Prognosis_Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_line(aes(group = Patient_ID), color = "gray40",
              linewidth = 0.6, alpha = 0.7) +
    geom_point(shape = 21, color = "black", size = 2) +
    annotate("segment", x = 1, xend = 2,
             y = ymax * 1.08, yend = ymax * 1.08,
             color = "black", linewidth = 0.5) +
    annotate("segment", x = 1, xend = 1,
             y = ymax * 1.08 * 0.995, yend = ymax * 1.08,
             color = "black", linewidth = 0.5) +
    annotate("segment", x = 2, xend = 2,
             y = ymax * 1.08 * 0.995, yend = ymax * 1.08,
             color = "black", linewidth = 0.5) +
    annotate("text", x = 1.5, y = ymax * 1.16,
             label = signif_label(pval), size = 4.5, color = "black") +
    scale_fill_manual(values = my_pub_colors, breaks = GROUP_LEVELS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(x = NULL, y = "% of memory B\namong PBMCs", title = "") +
    facet_wrap(~ "memory B") +
    theme_classic(base_size = 12) +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y     = element_text(size = 10),
          strip.text      = element_text(face = "bold", size = 12),
          strip.background = element_rect(fill = "grey90", color = "black",
                                          linewidth = 1),
          panel.border    = element_rect(color = "black", fill = NA,
                                         linewidth = 0.8),
          legend.position = "none")
}

p_memB <- make_memB_panel(df_memB, c("Poor_d1","Poor_d7")) /
  make_memB_panel(df_memB, c("Good_d1","Good_d7")) +
  plot_layout(ncol = 1, nrow = 2)

ggsave("MemoryB_proportion_stacked.png", p_memB,
       width = 2.8, height = 6, dpi = 300, bg = "white")
message("Saved: MemoryB_proportion_stacked.png")

# =============================================================================
# SECTION 3: B Cell Subset Proportion within B Cells
# =============================================================================
message("\n>>> B cell subset proportion (within B cells)")

B$Prognosis_Group <- factor(as.character(B$Prognosis_Group),
                            levels = GROUP_LEVELS)

B_prop <- B@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    fine            = trimws(as.character(predicted.celltype.major_fineB)),
    orig.ident      = as.character(orig.ident),
    Prognosis_Group = as.character(Prognosis_Group)
  ) %>%
  dplyr::filter(!is.na(fine), fine %in% B_SUB_LEVELS) %>%
  dplyr::group_by(orig.ident, Prognosis_Group, fine) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(orig.ident, Prognosis_Group) %>%
  dplyr::mutate(percent = n_cells / sum(n_cells) * 100) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Prognosis_Group = factor(Prognosis_Group, levels = GROUP_LEVELS),
    fine            = factor(fine, levels = B_SUB_LEVELS),
    Patient_ID      = sub("-d[17]$", "", orig.ident)
  )

make_paired_B_plot <- function(B_prop, prognosis, out_png) {
  groups_use <- if (prognosis == "Poor") c("Poor_d1","Poor_d7") else
    c("Good_d1","Good_d7")
  g1 <- groups_use[1]; g2 <- groups_use[2]
  
  df_sub <- B_prop %>%
    dplyr::filter(Prognosis_Group %in% groups_use) %>%
    dplyr::mutate(Prognosis_Group = factor(as.character(Prognosis_Group),
                                           levels = groups_use))
  
  wide <- df_sub %>%
    dplyr::select(fine, Patient_ID, Prognosis_Group, percent) %>%
    tidyr::pivot_wider(names_from = Prognosis_Group, values_from = percent)
  
  pvals <- wide %>%
    dplyr::group_by(fine) %>%
    dplyr::summarise(
      n_pairs = sum(!is.na(.data[[g1]]) & !is.na(.data[[g2]])),
      p.value = if (sum(!is.na(.data[[g1]]) & !is.na(.data[[g2]])) >= 2)
        t.test(.data[[g2]], .data[[g1]], paired = TRUE)$p.value else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(signif_label = signif_label(p.value))
  
  ypos <- df_sub %>%
    dplyr::group_by(fine) %>%
    dplyr::summarise(
      y_bar   = {mx <- suppressWarnings(max(percent, na.rm=TRUE));
      if (is.finite(mx)) mx*1.08 else 1},
      y_label = {mx <- suppressWarnings(max(percent, na.rm=TRUE));
      if (is.finite(mx)) mx*1.15 else 1},
      .groups = "drop"
    )
  pvals <- dplyr::left_join(pvals, ypos, by = "fine") %>%
    dplyr::mutate(x_left = 1, x_right = 2, x_mid = 1.5)
  
  p <- ggplot(df_sub, aes(x = Prognosis_Group, y = percent,
                          fill = Prognosis_Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_line(aes(group = Patient_ID), color = "gray40",
              linewidth = 0.6, alpha = 0.7) +
    geom_point(shape = 21, color = "black", size = 2) +
    facet_wrap(~ fine, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = my_pub_colors, breaks = groups_use) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20)))
  
  p <- add_sig_bracket(p, pvals %>% dplyr::rename(signif = signif_label),
                       y_bar_col = "y_bar", y_label_col = "y_label",
                       label_col = "signif")
  
  p <- p +
    labs(x = NULL, y = "Percentage within B cells (%)", title = "") +
    theme_classic(base_size = 12) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y      = element_text(size = 10),
          strip.text       = element_text(face = "bold", size = 12),
          strip.background = element_rect(color = "black", fill = "grey90",
                                          linewidth = 1),
          legend.position  = "none")
  
  ggsave(out_png, p, width = 7, height = 3.5, dpi = 300)
  message("Saved: ", out_png)
  invisible(p)
}

make_paired_B_plot(B_prop, "Poor", "B_subset_proportion_Poor.png")
make_paired_B_plot(B_prop, "Good", "B_subset_proportion_Good.png")

# =============================================================================
# SECTION 4: B Cell UMAP
# =============================================================================
message("\n>>> B cell UMAP")

B$predicted.celltype.major_fineB <- factor(
  B$predicted.celltype.major_fineB, levels = B_SUB_LEVELS
)
B$Prognosis_Group <- factor(as.character(B$Prognosis_Group),
                            levels = GROUP_LEVELS)

# Overall UMAP
p_b_umap <- DimPlot(B, reduction = "umap",
                    group.by = "predicted.celltype.major_fineB",
                    cols = my_b_colors, label = FALSE,
                    pt.size = 0.6, raster = FALSE) +
  ggtitle("B") +
  theme_classic(base_size = 14) +
  theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
        panel.border    = element_blank(),
        axis.line       = element_line(color = "black", linewidth = 0.8),
        axis.text       = element_text(color = "black", size = 12),
        legend.title    = element_blank())

ggsave("B_UMAP.png", p_b_umap, width = 6, height = 5, dpi = 300)
message("Saved: B_UMAP.png")

# Split UMAP (2x2)
p_b_split <- DimPlot(B, reduction = "umap",
                     group.by = "predicted.celltype.major_fineB",
                     split.by = "Prognosis_Group", ncol = 2,
                     cols = my_b_colors, label = FALSE,
                     pt.size = 0.8, raster = FALSE) +
  ggtitle("") +
  theme_classic(base_size = 14) +
  theme(panel.border    = element_blank(),
        axis.line       = element_line(color = "black", linewidth = 0.8),
        axis.text       = element_text(color = "black", size = 12),
        strip.text      = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        legend.title    = element_blank(),
        legend.text     = element_text(size = 12))

ggsave("B_cell_UMAP_Split_2x2.png", p_b_split,
       width = 8, height = 7, dpi = 600)
message("Saved: B_cell_UMAP_Split_2x2.png")

# =============================================================================
# SECTION 5: B Cell DEG
# =============================================================================
message("\n>>> B cell DEG analysis")

OBJ_B <- B
DefaultAssay(OBJ_B) <- "RNA"
OBJ_B <- safe_join_layers(OBJ_B)
OBJ_B$Prognosis_Group <- factor(OBJ_B$Prognosis_Group, levels = GROUP_LEVELS)
Idents(OBJ_B) <- "Prognosis_Group"
print(table(Idents(OBJ_B)))

do_deg_pair <- function(obj, ident1, ident2) {
  FindMarkers(obj, ident.1 = ident1, ident.2 = ident2,
              logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox") %>%
    tibble::rownames_to_column("gene")
}

get_up_genes <- function(deg_df, lfc_cut = FC_CUT, padj_cut = PADJ_CUT) {
  deg_df <- if (!"gene" %in% colnames(deg_df))
    tibble::rownames_to_column(deg_df, "gene") else deg_df
  list(
    up_ident1 = deg_df %>%
      dplyr::filter(p_val_adj < padj_cut, avg_log2FC >=  lfc_cut) %>%
      dplyr::pull(gene),
    up_ident2 = deg_df %>%
      dplyr::filter(p_val_adj < padj_cut, avg_log2FC <= -lfc_cut) %>%
      dplyr::pull(gene)
  )
}

print_deg_summary <- function(title, id1, id2, genes_up1, genes_up2) {
  cat("\n==", title, "==\n")
  cat("  Up in", id1, "(n =", length(genes_up1), ")\n")
  if (length(genes_up1) > 0)
    cat("   ", paste(genes_up1, collapse = ", "), "\n")
  cat("  Up in", id2, "(n =", length(genes_up2), ")\n")
  if (length(genes_up2) > 0)
    cat("   ", paste(genes_up2, collapse = ", "), "\n\n")
}

# Run DEG for all comparisons
message("Running DEG comparisons...")
deg_poor_d7_d1   <- do_deg_pair(OBJ_B, "Poor_d7",  "Poor_d1")
deg_good_d7_d1   <- do_deg_pair(OBJ_B, "Good_d7",  "Good_d1")
deg_poor_d1_gd1  <- do_deg_pair(OBJ_B, "Poor_d1",  "Good_d1")
deg_poor_d7_gd7  <- do_deg_pair(OBJ_B, "Poor_d7",  "Good_d7")

sig_poor   <- get_up_genes(deg_poor_d7_d1)
sig_good   <- get_up_genes(deg_good_d7_d1)
sig_pd1_gd1 <- get_up_genes(deg_poor_d1_gd1)
sig_pd7_gd7 <- get_up_genes(deg_poor_d7_gd7)

cat("\n# [B cells] DEG summary (|log2FC|>=1, padj<0.05)\n")
print_deg_summary("Poor_d7 vs Poor_d1", "Poor_d7","Poor_d1",
                  sig_poor$up_ident1, sig_poor$up_ident2)
print_deg_summary("Good_d7 vs Good_d1", "Good_d7","Good_d1",
                  sig_good$up_ident1, sig_good$up_ident2)
print_deg_summary("Poor_d1 vs Good_d1", "Poor_d1","Good_d1",
                  sig_pd1_gd1$up_ident1, sig_pd1_gd1$up_ident2)
print_deg_summary("Poor_d7 vs Good_d7", "Poor_d7","Good_d7",
                  sig_pd7_gd7$up_ident1, sig_pd7_gd7$up_ident2)

# =============================================================================
# SECTION 6: B Cell Volcano Plots
# =============================================================================
message("\n>>> B cell volcano plots")

# 6-A: All B — Poor_d7 vs Poor_d1
build_volcano(deg_poor_d7_d1,
              "B cell: Poor_d7 vs Poor_d1",
              "Volcano_B_Poor_d7_vs_Poor_d1.png",
              GENES_UP_POOR_D7, GENES_DOWN_POOR_D7)

# 6-B: All B — Poor_d7 vs all other groups
deg_poor_d7_others <- FindMarkers(
  OBJ_B,
  ident.1 = "Poor_d7", ident.2 = c("Poor_d1","Good_d1","Good_d7"),
  logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox"
) %>% tibble::rownames_to_column("gene")

build_volcano(deg_poor_d7_others,
              "B cell: Poor_d7 vs All Other Groups",
              "Volcano_B_PoorD7_vs_Others.png",
              GENES_UP_POOR_D7, GENES_DOWN_POOR_D7)

# 6-C: Memory B — Poor_d7 vs all other memory B
OBJ_memB <- subset(OBJ_B, subset = predicted.celltype.major_fineB == "memory B")
cat("Memory B cell count:", ncol(OBJ_memB), "\n")
print(table(OBJ_memB$Prognosis_Group))

OBJ_memB$DEG_group <- ifelse(OBJ_memB$Prognosis_Group == "Poor_d7",
                             "Poor_d7_memB", "other_memB")
Idents(OBJ_memB) <- "DEG_group"

deg_memB <- FindMarkers(OBJ_memB,
                        ident.1 = "Poor_d7_memB", ident.2 = "other_memB",
                        logfc.threshold = 0, min.pct = 0.1,
                        test.use = "wilcox") %>%
  tibble::rownames_to_column("gene")

GENES_UP_MEMB <- c("NR4A1","NR4A2","NR4A3","STAT3","DDIT4","GADD45B",
                   "SLC2A3","AREG","JUND","TNFSF9","HES4","GPR183",
                   "CXCR5","ZFP36")
GENES_DOWN_MEMB <- c("IFI44L","IFI27","MX1","MX2","ISG15","OAS1","OAS2",
                     "IRF7","STAT1","XAF1","IFI6","SAMD9L",
                     "CCR7","SELL","FCMR")

build_volcano(deg_memB,
              "Memory B: Poor_d7 vs All Other Groups",
              "Volcano_MemB_PoorD7_vs_Others.png",
              GENES_UP_MEMB, GENES_DOWN_MEMB)

# =============================================================================
# SECTION 7: GO Enrichment – B cells (Poor_d7 vs Poor_d1)
# =============================================================================
message("\n>>> GO enrichment: B cells (Poor_d7 vs Poor_d1)")

ego_b_up   <- run_go_enrichment(sig_poor$up_ident1,
                                "B_DEG_Poor_d7_up")
ego_b_down <- run_go_enrichment(sig_poor$up_ident2,
                                "B_DEG_Poor_d1_up")

cat("GO terms Poor_d7 up:", if (is.null(ego_b_up))   0 else
  nrow(as.data.frame(ego_b_up)), "\n")
cat("GO terms Poor_d1 up:", if (is.null(ego_b_down))  0 else
  nrow(as.data.frame(ego_b_down)), "\n")

save_go_bidirectional(
  ego_b_up, ego_b_down,
  comp_name = "B: GO BP enrichment (Poor_d7 vs Poor_d1)",
  out_png   = "B_GO_BP_Poor_d7_vs_Poor_d1.png",
  top_n_up = 9, top_n_down = 9, width = 17, height = 11
)

# GO enrichment – Memory B (Poor_d7 vs others)
sig_memB <- get_up_genes(deg_memB)
ego_memB_up   <- run_go_enrichment(sig_memB$up_ident1,
                                   "MemoryB_PoorD7_up")
ego_memB_down <- run_go_enrichment(sig_memB$up_ident2,
                                   "MemoryB_others_up")

save_go_bidirectional(
  ego_memB_up, ego_memB_down,
  comp_name = "Memory B: GO BP enrichment (Poor_d7 vs All Other Groups)",
  out_png   = "MemoryB_GO_BP_PoorD7_vs_Others.png",
  top_n_up = 7, top_n_down = 7, width = 18, height = 10
)

# =============================================================================
# SECTION 8: TNF / IL6 Module Score Violin Plots
# =============================================================================
message("\n>>> TNF/IL6 module scores (all B & memory B)")

DefaultAssay(merged_int) <- "RNA"
merged_int <- safe_join_layers(merged_int)
obj_b <- subset(merged_int, subset = predicted.celltype.major == "B")
obj_b$Prognosis_Group <- factor(as.character(obj_b$Prognosis_Group),
                                levels = GROUP_LEVELS)

msig_h    <- msigdbr(species = "Homo sapiens", category = "H")
genes_tnf <- msig_h %>%
  dplyr::filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>%
  dplyr::pull(gene_symbol) %>% unique()
genes_il6 <- msig_h %>%
  dplyr::filter(gs_name == "HALLMARK_IL6_JAK_STAT3_SIGNALING") %>%
  dplyr::pull(gene_symbol) %>% unique()

obj_b <- AddModuleScore(obj_b, list(intersect(genes_tnf, rownames(obj_b))),
                        name = "TNFA_Score_")
obj_b <- AddModuleScore(obj_b, list(intersect(genes_il6, rownames(obj_b))),
                        name = "IL6_Score_")
obj_b$TNFA_Score <- obj_b$TNFA_Score_1
obj_b$IL6_Score  <- obj_b$IL6_Score_1

COMP_B <- list(c("Poor_d1","Poor_d7"),
               c("Good_d1","Good_d7"),
               c("Poor_d7","Good_d7"))

make_violin_module(obj_b, "TNFA_Score",
                   "B: TNF\u03b1-NF-\u03baB signaling",
                   "B_Cells_TNFA_Hallmark_violin.png",
                   comparisons = COMP_B)

make_violin_module(obj_b, "IL6_Score",
                   "B: IL6-JAK-STAT3 signaling",
                   "B_Cells_IL6_JAKSTAT_Hallmark_violin.png",
                   comparisons = COMP_B)

# Memory B subset
obj_memb <- subset(obj_b, subset = predicted.celltype.major_fineB == "memory B")
cat("Memory B count:", ncol(obj_memb), "\n")
print(table(obj_memb$Prognosis_Group))

make_violin_module(obj_memb, "TNFA_Score",
                   "Memory B: TNF\u03b1-NF-\u03baB signaling",
                   "MemoryB_TNFA_Hallmark_violin.png",
                   comparisons = COMP_B, width = 5.5, height = 4.8)

make_violin_module(obj_memb, "IL6_Score",
                   "Memory B: IL6-JAK-STAT3 signaling",
                   "MemoryB_IL6_JAKSTAT_Hallmark_violin.png",
                   comparisons = COMP_B, width = 5.5, height = 4.8)

# =============================================================================
# SECTION 9: IL6 Expression DotPlot (Cell type × Prognosis)
# =============================================================================
message("\n>>> IL6 expression DotPlot across cell types")

merged_int$Celltype_Prog <- paste(merged_int$predicted.celltype.major,
                                  merged_int$Prognosis, sep = "__")

p_dot_raw <- DotPlot(merged_int, features = "IL6",
                     group.by = "Celltype_Prog")

df_dot <- p_dot_raw$data %>%
  tidyr::separate(id, into = c("Celltype","Prognosis"), sep = "__") %>%
  dplyr::mutate(Prognosis = factor(Prognosis, levels = c("Poor","Good")))

p_il6_dot <- ggplot(df_dot, aes(x = Prognosis, y = Celltype,
                                size = pct.exp, color = avg.exp.scaled)) +
  geom_point(alpha = 0.9) +
  scale_color_viridis_c(option = "inferno", begin = 0.2, end = 0.95,
                        name = "Avg. expression\n(scaled)") +
  scale_size(range = c(1.2, 6), name = "% expressed") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        plot.title  = element_text(face = "bold", hjust = 0.5)) +
  labs(x = "Prognosis", y = NULL, title = "IL6 expression across cell types")

ggsave("IL6_DotPlot_Prognosis_x_Celltype.png", p_il6_dot,
       width = 5.5, height = 5, dpi = 300)
message("Saved: IL6_DotPlot_Prognosis_x_Celltype.png")

message("\n>>> B cell analysis complete.")