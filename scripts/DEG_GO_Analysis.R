# =============================================================================
# 03_DEG_GO_Analysis.R
# COVID-19 scRNA-seq: Whole-PBMC DEG → Volcano → GO Enrichment →
#                     Quadrant Plot → DEG Heatmap
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggrepel)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pheatmap)
library(msigdbr)
library(cowplot)
library(stringr)
library(SeuratExtend)

setwd("~/covid_scRNA/covid_R/NEW/Figure")
dir.create(".", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 0. Load data & shared settings
# -----------------------------------------------------------------------------
obj <- merged_int
DefaultAssay(obj) <- "RNA"
suppressWarnings({ obj <- JoinLayers(obj, assay = "RNA") })

group_levels <- c("Poor_d1", "Poor_d7", "Good_d1", "Good_d7")
obj$Prognosis_Group <- factor(obj$Prognosis_Group, levels = group_levels)
Idents(obj) <- "Prognosis_Group"

PADJ_CUT <- 0.05
FC_THR   <- 0.5

# =============================================================================
# SECTION 1: Helper Functions
# =============================================================================

# -----------------------------------------------------------------------------
# 1-A. DEG utility functions
# -----------------------------------------------------------------------------

standardize_deg_cols <- function(deg_df) {
  if (!"gene" %in% colnames(deg_df))
    deg_df <- deg_df %>% tibble::rownames_to_column("gene")
  if (!"avg_log2FC" %in% colnames(deg_df) && "log2FC" %in% colnames(deg_df))
    deg_df <- dplyr::rename(deg_df, avg_log2FC = log2FC)
  if (!"p_val_adj" %in% colnames(deg_df) && "p_val" %in% colnames(deg_df))
    deg_df$p_val_adj <- deg_df$p_val
  deg_df
}

prep_deg_for_volcano <- function(deg_df) {
  deg_df <- standardize_deg_cols(deg_df)
  deg_plot <- deg_df %>%
    dplyr::filter(!is.na(p_val_adj)) %>%
    dplyr::mutate(
      p_val_adj_plot = ifelse(p_val_adj == 0, 1e-300, p_val_adj),
      sig_fc         = (abs(avg_log2FC) >= FC_THR & p_val_adj < PADJ_CUT)
    )
  list(deg_all = deg_df, deg_plot = deg_plot)
}

run_deg_pair_safe <- function(seu,
                              ident1,
                              ident2,
                              group_col = "Prognosis_Group",
                              min.pct   = 0.1,
                              logfc.thr = 0,
                              test.use  = "wilcox") {
  ident1 <- trimws(as.character(ident1))
  ident2 <- trimws(as.character(ident2))
  
  grp     <- trimws(as.character(seu[[group_col]][, 1]))
  tmp_col <- "TMP_pair"
  
  seu[[tmp_col]] <- NA_character_
  seu[[tmp_col]][, 1] <- ifelse(grp == ident1, ident1,
                                ifelse(grp == ident2, ident2, NA))
  
  keep    <- !is.na(seu[[tmp_col]][, 1])
  seu_sub <- subset(seu, cells = rownames(seu@meta.data)[keep])
  
  tab <- table(seu_sub[[tmp_col]][, 1], useNA = "ifany")
  message("[DEG] ", ident1, " vs ", ident2, " => ",
          paste(names(tab), tab, sep = ":", collapse = "  "))
  
  Idents(seu_sub) <- tmp_col
  deg <- tryCatch({
    Seurat::FindMarkers(seu_sub,
                        ident.1         = ident1,
                        ident.2         = ident2,
                        min.pct         = min.pct,
                        logfc.threshold = logfc.thr,
                        test.use        = test.use) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene")
  }, error = function(e) {
    message("FindMarkers ERROR: ", conditionMessage(e)); NULL
  })
  
  if (!is.null(deg) && nrow(deg) > 0) {
    if (!"avg_log2FC" %in% colnames(deg) && "log2FC" %in% colnames(deg))
      deg <- dplyr::rename(deg, avg_log2FC = log2FC)
    if (!"p_val_adj" %in% colnames(deg) && "p_val" %in% colnames(deg))
      deg$p_val_adj <- deg$p_val
  }
  deg
}

print_deg_summary <- function(deg_df, comp_name) {
  deg_df <- standardize_deg_cols(deg_df)
  
  cat("\n==============================\n")
  cat("DEG Summary: ", comp_name, "\n", sep = "")
  cat("==============================\n")
  
  n_sig <- sum(!is.na(deg_df$p_val_adj) &
                 deg_df$p_val_adj < PADJ_CUT &
                 abs(deg_df$avg_log2FC) >= FC_THR)
  cat("Sig (padj<0.05 & |log2FC|>=0.5): ", n_sig, "\n", sep = "")
  
  up_tbl <- deg_df %>%
    dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.01, avg_log2FC >= 1) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
    dplyr::slice_head(n = 100)
  
  down_tbl <- deg_df %>%
    dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.01, avg_log2FC <= -1) %>%
    dplyr::arrange(avg_log2FC) %>%
    dplyr::slice_head(n = 100)
  
  cat("\n--- Top UP genes (padj<0.01 & log2FC>=1) ---\n")
  if (nrow(up_tbl) == 0) cat("(no genes)\n") else print(head(up_tbl, 10))
  
  cat("\n--- Top DOWN genes (padj<0.01 & log2FC<=-1) ---\n")
  if (nrow(down_tbl) == 0) cat("(no genes)\n") else print(head(down_tbl, 10))
  
  cat("\n[LIST] Top UP Genes:\n")
  if (nrow(up_tbl) > 0)
    cat(paste(paste0('"', up_tbl$gene, '"'), collapse = ", "), "\n")
  else cat("None\n")
  
  cat("\n[LIST] Top DOWN Genes:\n")
  if (nrow(down_tbl) > 0)
    cat(paste(paste0('"', down_tbl$gene, '"'), collapse = ", "), "\n")
  else cat("None\n")
  
  invisible(list(up_tbl = up_tbl, down_tbl = down_tbl))
}

# -----------------------------------------------------------------------------
# 1-B. Volcano plot function
# -----------------------------------------------------------------------------

plot_volcano <- function(deg_plot,
                         title_text,
                         out_png,
                         genes_interest_up,
                         genes_interest_down,
                         png_w, png_h,
                         dpi             = 300,
                         labSize         = 5.0,
                         pointSize       = 2.0,
                         legend_down_text,
                         legend_up_text,
                         n_top_label     = 20,
                         legend_pos_x    = 0.95,
                         legend_pos_y    = 0.95,
                         symmetric_x     = TRUE) {
  
  stopifnot(nrow(deg_plot) > 0)
  
  # Color assignment
  deg_plot$color_custom <- "grey80"
  deg_plot$color_custom[abs(deg_plot$avg_log2FC) > FC_THR &
                          deg_plot$p_val_adj < PADJ_CUT &
                          deg_plot$avg_log2FC > 0] <- "#B31B1B"
  deg_plot$color_custom[abs(deg_plot$avg_log2FC) > FC_THR &
                          deg_plot$p_val_adj < PADJ_CUT &
                          deg_plot$avg_log2FC < 0] <- "#0047AB"
  
  col_custom_vec <- setNames(deg_plot$color_custom, deg_plot$gene)
  
  # Label genes: manual + top auto
  manual_genes <- intersect(c(genes_interest_up, genes_interest_down),
                            deg_plot$gene)
  
  exclude_pat <- "^MT-|^MT|^AC|^XIST"
  top_up <- deg_plot %>%
    dplyr::filter(avg_log2FC > FC_THR, p_val_adj < PADJ_CUT,
                  !gene %in% manual_genes, !grepl(exclude_pat, gene)) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC))) %>%
    head(round(n_top_label / 2)) %>% dplyr::pull(gene)
  
  top_down <- deg_plot %>%
    dplyr::filter(avg_log2FC < -FC_THR, p_val_adj < PADJ_CUT,
                  !gene %in% manual_genes, !grepl(exclude_pat, gene)) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC))) %>%
    head(round(n_top_label / 2)) %>% dplyr::pull(gene)
  
  final_label_genes <- unique(c(manual_genes, top_up, top_down))
  
  # Axis limits
  if (symmetric_x) {
    limit <- max(abs(deg_plot$avg_log2FC), na.rm = TRUE)
    if (!is.finite(limit) || limit == 0) limit <- 1
    xlim_use <- c(-limit * 1.1, limit * 1.1)
  } else {
    min_fc <- min(deg_plot$avg_log2FC, na.rm = TRUE)
    max_fc <- max(deg_plot$avg_log2FC, na.rm = TRUE)
    if (max_fc < 1)  max_fc <- 1
    if (min_fc > -1) min_fc <- -1
    xlim_use <- c(min_fc * 1.1, max_fc * 1.1)
  }
  
  ymax <- max(-log10(deg_plot$p_val_adj_plot), na.rm = TRUE)
  if (!is.finite(ymax) || ymax == 0) ymax <- 1
  y_lim_use <- c(0, ymax * 1.05)
  
  png(out_png, width = png_w, height = png_h, res = dpi)
  
  p <- EnhancedVolcano(
    deg_plot,
    lab              = deg_plot$gene,
    selectLab        = final_label_genes,
    x                = "avg_log2FC",
    y                = "p_val_adj_plot",
    xlab             = bquote(~Log[2] ~ "Fold Change (Cutoff = " * .(FC_THR) * ")"),
    ylab             = bquote(~-Log[10] ~ italic(P)[adj]),
    title            = title_text,
    subtitle         = NULL,
    pCutoff          = PADJ_CUT,
    FCcutoff         = FC_THR,
    labSize          = labSize,
    pointSize        = pointSize,
    colAlpha         = 0.6,
    boxedLabels      = TRUE,
    drawConnectors   = TRUE,
    widthConnectors  = 0.6,
    typeConnectors   = "open",
    arrowheads       = FALSE,
    max.overlaps     = 150,
    legendPosition   = "none",
    caption          = NULL,
    ylim             = y_lim_use,
    xlim             = xlim_use,
    col              = rep("grey80", 4),
    colCustom        = col_custom_vec,
    gridlines.major  = FALSE,
    gridlines.minor  = FALSE,
    border           = "full",
    borderWidth      = 1.0,
    borderColour     = "black"
  ) +
    geom_vline(xintercept = c(-FC_THR, FC_THR),
               linetype = "longdash", color = "black", linewidth = 0.6) +
    theme(
      plot.title   = element_text(size = 24, face = "bold", hjust = 0.5,
                                  margin = margin(b = 15)),
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
      axis.text.x  = element_text(size = 18, color = "black"),
      axis.text.y  = element_text(size = 18, color = "black"),
      axis.ticks   = element_line(color = "black", linewidth = 0.8),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)
    )
  
  y_coord       <- y_lim_use[2] * legend_pos_y
  x_coord_left  <- xlim_use[1]  * legend_pos_x
  x_coord_right <- xlim_use[2]  * legend_pos_x
  
  p <- p +
    annotate("text", x = x_coord_left,  y = y_coord,
             label = legend_down_text,
             hjust = 0, size = 6, fontface = "bold", color = "#0047AB") +
    annotate("text", x = x_coord_right, y = y_coord,
             label = legend_up_text,
             hjust = 1, size = 6, fontface = "bold", color = "#B31B1B")
  
  print(p)
  dev.off()
  message("  -> Saved: ", out_png)
  invisible(p)
}

# -----------------------------------------------------------------------------
# 1-C. GO enrichment functions
# -----------------------------------------------------------------------------

run_go_enrich_bp <- function(genes,
                             ont             = "BP",
                             p_cutoff        = 0.05,
                             q_cutoff        = 0.2,
                             simplify_cutoff = 0.7) {
  genes <- unique(genes[!is.na(genes) & genes != ""])
  if (length(genes) < 3) return(NULL)
  
  ego <- tryCatch({
    clusterProfiler::enrichGO(
      gene          = genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff,
      readable      = TRUE
    )
  }, error = function(e) NULL)
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  
  ego <- tryCatch(enrichplot::pairwise_termsim(ego), error = function(e) ego)
  ego <- tryCatch(
    clusterProfiler::simplify(ego, cutoff = simplify_cutoff,
                              by = "p.adjust", select_fun = min),
    error = function(e) ego
  )
  if (nrow(as.data.frame(ego)) == 0) return(NULL)
  ego
}

save_go_bidirectional_barplot <- function(ego_up,
                                          ego_down,
                                          title_txt,
                                          out_png,
                                          top_n_up    = 10,
                                          top_n_down  = 10,
                                          wrap_width  = 45,
                                          sym_limit   = NULL,
                                          width       = 28,
                                          height      = 18) {
  if (is.null(ego_up) && is.null(ego_down)) {
    warning("[", title_txt, "] Both ego_up and ego_down are NULL. Skipping.")
    return(invisible(NULL))
  }
  
  df_up <- df_down <- NULL
  
  if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
    df_up <- as.data.frame(ego_up) %>%
      dplyr::arrange(p.adjust) %>% head(top_n_up) %>%
      dplyr::mutate(direction = "Up in D7",
                    score     = -log10(p.adjust))
  }
  
  if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
    df_down <- as.data.frame(ego_down) %>%
      dplyr::arrange(p.adjust) %>% head(top_n_down) %>%
      dplyr::mutate(direction = "Down in D7",
                    score     = log10(p.adjust))   # negative
  }
  
  df_all <- dplyr::bind_rows(df_down, df_up)
  if (is.null(df_all) || nrow(df_all) == 0) {
    warning("[", title_txt, "] No GO terms to plot. Skipping.")
    return(invisible(NULL))
  }
  
  df_all <- df_all %>%
    dplyr::mutate(Description_wrapped =
                    stringr::str_wrap(Description, width = wrap_width))
  
  up_terms   <- df_all %>% dplyr::filter(direction == "Up in D7") %>%
    dplyr::arrange(score) %>% dplyr::pull(Description_wrapped) %>% unique()
  down_terms <- df_all %>% dplyr::filter(direction == "Down in D7") %>%
    dplyr::arrange(score) %>% dplyr::pull(Description_wrapped) %>% unique()
  
  df_all$Description_wrapped <- factor(df_all$Description_wrapped,
                                       levels = c(up_terms, rev(down_terms)))
  
  max_abs <- max(abs(df_all$score), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  lim <- if (is.null(sym_limit)) ceiling(max_abs) else sym_limit
  
  p <- ggplot(df_all, aes(x = score, y = Description_wrapped, fill = direction)) +
    geom_col(width = 0.85) +
    scale_x_continuous(limits = c(-lim, lim),
                       name   = "-log10(adjusted P value)") +
    scale_fill_manual(
      values = c("Up in D7" = "#B2182B", "Down in D7" = "#2166AC"),
      name   = NULL
    ) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    labs(title = title_txt, y = NULL) +
    theme_bw(base_size = 14) +
    theme(
      text            = element_text(face = "bold"),
      panel.grid      = element_blank(),
      legend.position = c(0.25, 0.15)
    )
  
  png(out_png, width = width * 100, height = height * 100, res = 300)
  print(p)
  dev.off()
  message("Saved GO barplot: ", out_png, " (+-", lim, ")")
  invisible(df_all)
}

# Combined DEG + GO runner for a single group pair
run_go_for_pair <- function(ident_d7,
                            ident_d1,
                            out_png,
                            title_txt,
                            top_n_up   = 8,
                            top_n_down = 8,
                            padj_cut   = PADJ_CUT,
                            fc_thr     = 0.25) {
  deg <- run_deg_pair_safe(seu    = obj,
                           ident1 = ident_d7,
                           ident2 = ident_d1)
  
  if (is.null(deg) || nrow(deg) == 0) {
    warning("[", title_txt, "] No DEG results.")
    return(invisible(NULL))
  }
  
  genes_up   <- deg %>% dplyr::filter(p_val_adj < padj_cut,
                                      avg_log2FC > fc_thr) %>%
    dplyr::pull(gene)
  genes_down <- deg %>% dplyr::filter(p_val_adj < padj_cut,
                                      avg_log2FC < -fc_thr) %>%
    dplyr::pull(gene)
  
  cat(sprintf("[%s] Up: %d genes, Down: %d genes\n",
              title_txt, length(genes_up), length(genes_down)))
  
  ego_up   <- run_go_enrich_bp(genes_up)
  ego_down <- run_go_enrich_bp(genes_down)
  
  save_go_bidirectional_barplot(
    ego_up    = ego_up,
    ego_down  = ego_down,
    title_txt = title_txt,
    out_png   = out_png,
    top_n_up  = top_n_up,
    top_n_down = top_n_down
  )
  
  invisible(list(deg = deg, ego_up = ego_up, ego_down = ego_down,
                 genes_up = genes_up, genes_down = genes_down))
}

# =============================================================================
# SECTION 2: Run DEG (Poor and Good groups)
# =============================================================================
message("\n>>> Running DEG: Poor_d7 vs Poor_d1")
deg_poor <- FindMarkers(obj,
                        ident.1         = "Poor_d7",
                        ident.2         = "Poor_d1",
                        min.pct         = 0.1,
                        logfc.threshold = 0,
                        test.use        = "wilcox") %>%
  tibble::rownames_to_column("gene")

message(">>> Running DEG: Good_d7 vs Good_d1")
deg_good <- FindMarkers(obj,
                        ident.1         = "Good_d7",
                        ident.2         = "Good_d1",
                        min.pct         = 0.1,
                        logfc.threshold = 0,
                        test.use        = "wilcox") %>%
  tibble::rownames_to_column("gene")

print_deg_summary(deg_poor, "Poor_d7 vs Poor_d1")
print_deg_summary(deg_good, "Good_d7 vs Good_d1")

# =============================================================================
# SECTION 3: Volcano Plots
# =============================================================================

# --- Gene lists of interest ---
# Poor group: persistent inflammation (up) / steroid response failure (down)
genes_up_poor <- c("JUN", "FOSB", "EGR1", "IER3", "IER2", "ICAM1")
genes_down_poor <- c("FKBP5", "ZBTB16", "VSIG4", "IL1R2", "ADAMTS2")

# Good group: metabolic reprogramming & recovery (up) / inflammation resolution (down)
genes_up_good   <- c("PDK4", "SESN3", "SLC4A7", "CD3D", "CD27", "LTB", "IL32")
genes_down_good <- c("CXCL8", "S100A8", "IFI27", "ISG15", "OAS1",
                     "SIGLEC1", "FCGR1A", "RETN", "IL18")

# --- Poor volcano ---
pp <- prep_deg_for_volcano(deg_poor)
plot_volcano(
  deg_plot         = pp$deg_plot,
  title_text       = "Poor: Pre vs Post-treatment",
  out_png          = "PBMC_DEG_Poor.png",
  genes_interest_up   = genes_up_poor,
  genes_interest_down = genes_down_poor,
  png_w            = 2800, png_h = 3200, dpi = 300,
  labSize          = 5.0,  pointSize = 2.5,
  legend_down_text = paste0("Down in d7\n(n=",
                            sum(pp$deg_plot$avg_log2FC < -FC_THR &
                                  pp$deg_plot$p_val_adj < PADJ_CUT), ")"),
  legend_up_text   = paste0("Up in d7\n(n=",
                            sum(pp$deg_plot$avg_log2FC > FC_THR &
                                  pp$deg_plot$p_val_adj < PADJ_CUT), ")"),
  n_top_label      = 10,
  symmetric_x      = FALSE,
  legend_pos_x     = 0.9,
  legend_pos_y     = 0.45
)

# --- Good volcano ---
gg <- prep_deg_for_volcano(deg_good)
plot_volcano(
  deg_plot         = gg$deg_plot,
  title_text       = "Good: Pre vs Post-treatment",
  out_png          = "PBMC_DEG_Good.png",
  genes_interest_up   = genes_up_good,
  genes_interest_down = genes_down_good,
  png_w            = 2800, png_h = 3200, dpi = 300,
  labSize          = 5.0,  pointSize = 2.5,
  legend_down_text = paste0("Down in d7\n(n=",
                            sum(gg$deg_plot$avg_log2FC < -FC_THR &
                                  gg$deg_plot$p_val_adj < PADJ_CUT), ")"),
  legend_up_text   = paste0("Up in d7\n(n=",
                            sum(gg$deg_plot$avg_log2FC > FC_THR &
                                  gg$deg_plot$p_val_adj < PADJ_CUT), ")"),
  n_top_label      = 10,
  symmetric_x      = FALSE,
  legend_pos_x     = 0.9,
  legend_pos_y     = 0.45
)

# =============================================================================
# SECTION 4: GO Enrichment (whole PBMC, bidirectional barplot)
# =============================================================================
message("\n>>> GO Enrichment: Poor group")
genes_poor_up   <- deg_poor %>%
  dplyr::filter(p_val_adj < PADJ_CUT, avg_log2FC >= FC_THR) %>%
  dplyr::pull(gene) %>% unique()
genes_poor_down <- deg_poor %>%
  dplyr::filter(p_val_adj < PADJ_CUT, avg_log2FC <= -FC_THR) %>%
  dplyr::pull(gene) %>% unique()

ego_poor_up   <- run_go_enrich_bp(genes_poor_up)
ego_poor_down <- run_go_enrich_bp(genes_poor_down)

cat("\n--- Top 10 UP Pathways (Poor) ---\n")
if (!is.null(ego_poor_up))
  print(head(as.data.frame(ego_poor_up)[, c("ID", "Description", "p.adjust")], 10))
cat("\n--- Top 10 DOWN Pathways (Poor) ---\n")
if (!is.null(ego_poor_down))
  print(head(as.data.frame(ego_poor_down)[, c("ID", "Description", "p.adjust")], 10))

save_go_bidirectional_barplot(
  ego_poor_up, ego_poor_down,
  title_txt = "Poor: Pre vs Post-treatment (GO BP)",
  out_png   = "GO_BP_Poor_d7_vs_d1.png"
)

message("\n>>> GO Enrichment: Good group")
genes_good_up   <- deg_good %>%
  dplyr::filter(p_val_adj < PADJ_CUT, avg_log2FC >= FC_THR) %>%
  dplyr::pull(gene) %>% unique()
genes_good_down <- deg_good %>%
  dplyr::filter(p_val_adj < PADJ_CUT, avg_log2FC <= -FC_THR) %>%
  dplyr::pull(gene) %>% unique()

ego_good_up   <- run_go_enrich_bp(genes_good_up)
ego_good_down <- run_go_enrich_bp(genes_good_down)

cat("\n--- Top 10 UP Pathways (Good) ---\n")
if (!is.null(ego_good_up))
  print(head(as.data.frame(ego_good_up)[, c("ID", "Description", "p.adjust")], 10))
cat("\n--- Top 10 DOWN Pathways (Good) ---\n")
if (!is.null(ego_good_down))
  print(head(as.data.frame(ego_good_down)[, c("ID", "Description", "p.adjust")], 10))

save_go_bidirectional_barplot(
  ego_good_up, ego_good_down,
  title_txt = "Good: Pre vs Post-treatment (GO BP)",
  out_png   = "GO_BP_Good_d7_vs_d1.png"
)

# =============================================================================
# SECTION 5: Quadrant DEG Plot
# =============================================================================
message("\n>>> Building quadrant DEG data frame")

df_all <- full_join(
  deg_good %>% dplyr::rename(log2FC_good = avg_log2FC, padj_good = p_val_adj),
  deg_poor %>% dplyr::rename(log2FC_poor = avg_log2FC, padj_poor = p_val_adj),
  by = "gene"
) %>%
  dplyr::filter(!grepl("^AC|^C1orf|^XIST", gene)) %>%
  dplyr::mutate(
    padj_good   = ifelse(is.na(padj_good),   1, padj_good),
    padj_poor   = ifelse(is.na(padj_poor),   1, padj_poor),
    log2FC_good = ifelse(is.na(log2FC_good), 0, log2FC_good),
    log2FC_poor = ifelse(is.na(log2FC_poor), 0, log2FC_poor)
  )

# Keep only genes significant in both groups
df_pass <- df_all %>%
  dplyr::filter(padj_good < PADJ_CUT, padj_poor < PADJ_CUT) %>%
  dplyr::mutate(
    in_center_box  = (abs(log2FC_good) < FC_THR & abs(log2FC_poor) < FC_THR),
    outside_center = !in_center_box,
    quad = dplyr::case_when(
      log2FC_good >= 0 & log2FC_poor >= 0 ~ "Q1",
      log2FC_good <  0 & log2FC_poor >= 0 ~ "Q2",
      log2FC_good <  0 & log2FC_poor <  0 ~ "Q3",
      log2FC_good >= 0 & log2FC_poor <  0 ~ "Q4"
    ),
    quad      = factor(quad, levels = c("Q1","Q2","Q3","Q4")),
    quad_plot = ifelse(outside_center, as.character(quad), "NS"),
    quad_plot = factor(quad_plot, levels = c("Q1","Q2","Q3","Q4","NS"))
  )

quad_cols <- c(Q1 = "#D95F02", Q2 = "#7570B3",
               Q3 = "#1B9E77", Q4 = "#80B1D3", NS = "grey80")

# --- 5-A: IL6 / TNF / S100 highlight plot ---
msig_h    <- msigdbr(species = "Homo sapiens", category = "H")
genes_il6 <- msig_h %>%
  dplyr::filter(gs_name == "HALLMARK_IL6_JAK_STAT3_SIGNALING") %>%
  dplyr::pull(gene_symbol) %>% unique()
genes_tnf <- msig_h %>%
  dplyr::filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>%
  dplyr::pull(gene_symbol) %>% unique()

df_sig_il6 <- df_pass %>% dplyr::filter(gene %in% genes_il6) %>%
  dplyr::mutate(signal_legend = "IL6 signaling")
df_sig_tnf <- df_pass %>% dplyr::filter(gene %in% genes_tnf) %>%
  dplyr::mutate(signal_legend = "TNF signaling")
df_sig_s100 <- df_pass %>% dplyr::filter(grepl("^S100", gene)) %>%
  dplyr::mutate(signal_legend = "S100")

df_sig_all <- dplyr::bind_rows(df_sig_il6, df_sig_tnf, df_sig_s100) %>%
  dplyr::mutate(signal_legend = factor(signal_legend,
                                       levels = c("IL6 signaling",
                                                  "TNF signaling", "S100"))) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

signal_cols <- c("IL6 signaling" = "#E41A1C",
                 "TNF signaling" = "#377EB8",
                 "S100"          = "#4DAF4A")

df_lab_main <- df_pass %>%
  dplyr::left_join(
    df_sig_all %>% dplyr::select(gene, signal_legend), by = "gene"
  ) %>%
  dplyr::filter(
    quad %in% c("Q2", "Q4"),
    signal_legend %in% c("IL6 signaling", "TNF signaling"),
    abs(log2FC_good) >= FC_THR,
    abs(log2FC_poor) >= FC_THR
  ) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

df_lab_s100 <- df_sig_s100 %>%
  dplyr::mutate(
    nudge_x = ifelse(log2FC_good < 0, -2.0, 2.0),
    nudge_y = ifelse(log2FC_poor < 0, -2.0, 2.0)
  )

x_max    <- max(abs(df_pass$log2FC_good), na.rm = TRUE)
y_max    <- max(abs(df_pass$log2FC_poor), na.rm = TRUE)
lim_use  <- max(x_max, y_max, 1)
lab_pos  <- 0.92 * lim_use

p_A <- ggplot(df_pass, aes(x = log2FC_good, y = log2FC_poor)) +
  annotate("rect",
           xmin = -FC_THR, xmax = FC_THR,
           ymin = -FC_THR, ymax = FC_THR,
           fill = "grey90", alpha = 0.6, color = NA) +
  geom_point(aes(color = quad_plot), size = 1.6, alpha = 0.75) +
  scale_color_manual(values = quad_cols, guide = "none") +
  geom_point(data = df_sig_all,
             aes(fill = signal_legend),
             shape = 21, size = 2.8, alpha = 0.95,
             color = "grey25", stroke = 0.35) +
  scale_fill_manual(values = signal_cols, name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4, alpha = 1))) +
  # S100 outer circle
  geom_point(data = df_sig_s100,
             aes(x = log2FC_good, y = log2FC_poor),
             shape = 21, size = 4.0, fill = NA,
             color = "#E31A1C", stroke = 1.1, alpha = 1) +
  ggrepel::geom_text_repel(
    data = df_lab_main, aes(label = gene),
    color = "black", size = 3.6,
    box.padding = 0.4, point.padding = 0.25,
    force = 2, max.overlaps = Inf,
    min.segment.length = 0, segment.color = "grey70", seed = 123
  ) +
  { if (nrow(df_lab_s100) > 0)
    ggrepel::geom_text_repel(
      data = df_lab_s100, aes(label = gene),
      color = "#E41A1C",
      nudge_x = df_lab_s100$nudge_x, nudge_y = df_lab_s100$nudge_y,
      box.padding = 0.9, point.padding = 0.7,
      segment.color = "#E41A1C", segment.size = 0.6,
      max.overlaps = Inf, min.segment.length = 0, seed = 123
    ) } +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = c(-FC_THR, FC_THR),
             linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = c(-FC_THR, FC_THR),
             linetype = "dotted", color = "grey60") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey70") +
  coord_cartesian(xlim = c(-lim_use, lim_use),
                  ylim = c(-lim_use, lim_use)) +
  annotate("text", x =  lab_pos, y =  lab_pos,
           label = "Poor d7 up\nGood d7 up",
           color = quad_cols[["Q1"]], fontface = "bold", size = 5, hjust = 1) +
  annotate("text", x = -lab_pos, y =  lab_pos,
           label = "Poor d7 up\nGood d7 down",
           color = quad_cols[["Q2"]], fontface = "bold", size = 5, hjust = 0) +
  annotate("text", x = -lab_pos, y = -lab_pos,
           label = "Poor d7 down\nGood d7 down",
           color = quad_cols[["Q3"]], fontface = "bold", size = 5, hjust = 0) +
  annotate("text", x =  lab_pos, y = -lab_pos,
           label = "Poor d7 down\nGood d7 up",
           color = quad_cols[["Q4"]], fontface = "bold", size = 5, hjust = 1) +
  labs(
    x     = paste0("log2FC (Good d7 vs Good d1)  [adj.p<", PADJ_CUT, "]"),
    y     = paste0("log2FC (Poor d7 vs Poor d1)  [adj.p<", PADJ_CUT, "]"),
    title = ""
  ) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

png("quadrant_DEG_S100_IL6_TNF.png", width = 2700, height = 2400, res = 300)
print(p_A)
dev.off()
message("Saved: quadrant_DEG_S100_IL6_TNF.png")

# --- 5-B: Q2/Q4 concept-labeled quadrant plot ---
df_q2q4_fc <- df_pass %>%
  dplyr::filter(quad %in% c("Q2","Q4"), outside_center,
                abs(log2FC_good) >= FC_THR, abs(log2FC_poor) >= FC_THR) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

cat("[Q2+Q4 FC-pass genes]: ", nrow(df_q2q4_fc), "\n")

# Genes of interest with biological concept labels
genes_A <- c("CD83", "MAFB", "THBD", "CAPG", "MARCKS")   # Myeloid activation
genes_B <- c("FABP5", "S100A12")                           # Inflammatory remodeling
genes_E <- c("LCK", "LAT", "CD2")                         # TCR signaling
genes_F <- c("KLRK1", "GZMA", "GZMH")                     # Cytotoxic effector
genes_G <- c("IL7R")                                       # T cell survival

genes_keep <- unique(c(genes_A, genes_B, genes_E, genes_F, genes_G))

df_focus <- df_q2q4_fc %>%
  dplyr::filter(gene %in% genes_keep) %>%
  dplyr::mutate(
    concept = dplyr::case_when(
      gene %in% genes_A ~ "Myeloid activation",
      gene %in% genes_B ~ "Inflammatory remodeling",
      gene %in% genes_E ~ "TCR signaling",
      gene %in% genes_F ~ "Cytotoxic effector",
      gene %in% genes_G ~ "T cell survival"
    ),
    concept = factor(concept,
                     levels = c("Myeloid activation", "Inflammatory remodeling",
                                "TCR signaling", "Cytotoxic effector",
                                "T cell survival"))
  )

missing_genes <- setdiff(genes_keep, df_focus$gene)
if (length(missing_genes) > 0)
  cat("[Missing genes after filters]:", paste(missing_genes, collapse = ", "), "\n")

concept_cols <- c(
  "Myeloid activation"      = "#4C78A8",
  "Inflammatory remodeling" = "#F58518",
  "TCR signaling"           = "#54A24B",
  "Cytotoxic effector"      = "#B279A2",
  "T cell survival"         = "#3B5B8A"
)

quad_cols_B <- c(Q1 = "#E08A00", Q2 = "#7570B3",
                 Q3 = "#1B9E77", Q4 = "#74ADD1")

df_pass <- df_pass %>%
  dplyr::mutate(
    quad_plot = dplyr::case_when(
      log2FC_good >= 0 & log2FC_poor >= 0 ~ "Q1",
      log2FC_good <  0 & log2FC_poor >= 0 ~ "Q2",
      log2FC_good <  0 & log2FC_poor <  0 ~ "Q3",
      log2FC_good >= 0 & log2FC_poor <  0 ~ "Q4"
    )
  )

lim_use <- 5.5
lab_pos <- lim_use * 0.88

# Nudge values for gene labels
nudge_x_vec <- dplyr::case_when(
  df_focus$gene == "CD83"    ~  0.0, df_focus$gene == "MARCKS"  ~  1.2,
  df_focus$gene == "FABP5"   ~ -1.0, df_focus$gene == "S100A12" ~  0.8,
  df_focus$gene == "THBD"    ~ -1.2, df_focus$gene == "MAFB"    ~ -1.2,
  df_focus$gene == "CAPG"    ~  1.2, df_focus$gene == "LCK"     ~  1.2,
  df_focus$gene == "LAT"     ~  1.2, df_focus$gene == "CD2"     ~  1.2,
  df_focus$gene == "GZMA"    ~  0.8, df_focus$gene == "GZMH"    ~ -0.8,
  df_focus$gene == "IL7R"    ~ -0.8, df_focus$gene == "KLRK1"   ~ -0.5,
  TRUE ~ 0
)
nudge_y_vec <- dplyr::case_when(
  df_focus$gene == "CD83"    ~  0.6, df_focus$gene == "MARCKS"  ~  0.7,
  df_focus$gene == "FABP5"   ~  0.3, df_focus$gene == "S100A12" ~  0.8,
  df_focus$gene == "THBD"    ~  0.2, df_focus$gene == "MAFB"    ~ -0.3,
  df_focus$gene == "CAPG"    ~  0.2, df_focus$gene == "LCK"     ~  0.3,
  df_focus$gene == "LAT"     ~ -0.3, df_focus$gene == "CD2"     ~ -0.6,
  df_focus$gene == "GZMA"    ~  0.3, df_focus$gene == "GZMH"    ~ -0.4,
  df_focus$gene == "IL7R"    ~ -0.5, df_focus$gene == "KLRK1"   ~ -0.7,
  TRUE ~ 0
)

p_main <- ggplot(df_pass, aes(x = log2FC_good, y = log2FC_poor)) +
  annotate("rect",
           xmin = -FC_THR, xmax = FC_THR,
           ymin = -FC_THR, ymax = FC_THR,
           fill = "grey95", alpha = 0.9, color = NA) +
  geom_point(aes(color = quad_plot), size = 1.6, alpha = 0.65) +
  scale_color_manual(values = quad_cols_B, guide = "none") +
  geom_point(data = df_focus, aes(fill = concept),
             shape = 21, size = 3.0, alpha = 0.95,
             color = "black", stroke = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = concept_cols, guide = "none") +
  ggrepel::geom_text_repel(
    data               = df_focus,
    aes(label = gene),
    size               = 3.0,
    fontface           = "plain",
    box.padding        = 0.5,
    point.padding      = 0.1,
    max.overlaps       = Inf,
    min.segment.length = 0,
    segment.color      = "black",
    segment.size       = 0.48,
    segment.alpha      = 0.8,
    force              = 4,
    force_pull         = 0.3,
    seed               = 42,
    show.legend        = FALSE,
    nudge_x            = nudge_x_vec,
    nudge_y            = nudge_y_vec
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = c(-FC_THR, FC_THR),
             linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = c(-FC_THR, FC_THR),
             linetype = "dotted", color = "grey60") +
  coord_cartesian(xlim = c(-lim_use, lim_use), ylim = c(-lim_use, lim_use)) +
  annotate("text", x =  lab_pos, y =  lab_pos,
           label = "Poor_d7 up\nGood_d7 up",
           color = quad_cols_B[["Q1"]], fontface = "plain", size = 3.5,
           hjust = 1, vjust = 1) +
  annotate("text", x = -lab_pos, y =  lab_pos,
           label = "Poor_d7 up\nGood_d7 down",
           color = quad_cols_B[["Q2"]], fontface = "plain", size = 3.5,
           hjust = 0, vjust = 1) +
  annotate("text", x = -lab_pos, y = -lab_pos,
           label = "Poor_d7 down\nGood_d7 down",
           color = quad_cols_B[["Q3"]], fontface = "plain", size = 3.5,
           hjust = 0, vjust = 0) +
  annotate("text", x =  lab_pos, y = -lab_pos,
           label = "Poor_d7 down\nGood_d7 up",
           color = quad_cols_B[["Q4"]], fontface = "plain", size = 3.5,
           hjust = 1, vjust = 0) +
  labs(
    x     = paste0("log2FC (Good_d7 vs Good_d1)  [adj.p<", PADJ_CUT, "]"),
    y     = paste0("log2FC (Poor_d7 vs Poor_d1)  [adj.p<", PADJ_CUT, "]"),
    title = ""
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid   = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text    = element_text(size = 10))

# Separate legends for Q2 and Q4 concepts
make_concept_legend <- function(df_sub, concept_levels) {
  p_leg <- ggplot(df_sub) +
    geom_point(aes(x = 1, y = 1, fill = concept),
               shape = 21, size = 3.5, color = "black", stroke = 0.25) +
    scale_fill_manual(values  = concept_cols,
                      breaks  = concept_levels,
                      drop    = FALSE,
                      name    = "") +
    guides(fill = guide_legend(title    = NULL,
                               byrow    = TRUE,
                               override.aes = list(shape = 21, size = 3.5,
                                                   alpha = 1, color = "black"))) +
    theme_void(base_size = 11) +
    theme(legend.position   = "right",
          legend.text       = element_text(size = 10),
          legend.key.height = grid::unit(0.5, "cm"),
          legend.spacing.y  = grid::unit(0.3, "cm"))
  cowplot::get_legend(p_leg)
}

leg_q2 <- make_concept_legend(
  df_focus %>% dplyr::filter(concept %in% c("Myeloid activation",
                                            "Inflammatory remodeling")),
  c("Myeloid activation", "Inflammatory remodeling")
)
leg_q4 <- make_concept_legend(
  df_focus %>% dplyr::filter(concept %in% c("TCR signaling",
                                            "Cytotoxic effector",
                                            "T cell survival")),
  c("TCR signaling", "Cytotoxic effector", "T cell survival")
)

p_final <- cowplot::ggdraw() +
  cowplot::draw_plot(p_main) +
  cowplot::draw_grob(leg_q2, x = 0.15, y = 0.62, width = 0.28, height = 0.22) +
  cowplot::draw_grob(leg_q4, x = 0.55, y = 0.15, width = 0.34, height = 0.26)

png("quadrant_DEG.png", width = 2000, height = 1900, res = 300)
print(p_final)
dev.off()
message("Saved: quadrant_DEG.png")

# =============================================================================
# SECTION 6: DEG Heatmap (Q2 / Q4 genes across cell types)
# =============================================================================
message("\n>>> Building DEG heatmap")

celltype_col    <- "predicted.celltype.major"
celltype_levels <- c("CD14 Mono", "CD16 Mono", "DC", "B",
                     "NK", "CD4 T", "CD8 T", "other T")

obj[[celltype_col]] <- factor(as.character(obj[[celltype_col]][, 1]),
                              levels = celltype_levels)

# Q2/Q4 gene set from df_pass (already computed in Section 5)
df_q2q4_heat <- df_pass %>%
  dplyr::mutate(score_rank = abs(log2FC_good) + abs(log2FC_poor)) %>%
  dplyr::filter(quad %in% c("Q2", "Q4"), outside_center,
                abs(log2FC_good) >= FC_THR, abs(log2FC_poor) >= FC_THR) %>%
  dplyr::arrange(dplyr::desc(score_rank)) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

genes_q2  <- df_q2q4_heat %>% dplyr::filter(quad == "Q2") %>% dplyr::pull(gene)
genes_q4  <- df_q2q4_heat %>% dplyr::filter(quad == "Q4") %>% dplyr::pull(gene)
genes_use <- unique(c(genes_q2, genes_q4))

cat("[Q2 genes] N =", length(genes_q2), "\n")
cat("[Q4 genes] N =", length(genes_q4), "\n")

# Mean expression matrix (genes x cell type)
mat_mean <- SeuratExtend::CalcStats(
  obj,
  features = genes_use,
  group.by = celltype_col,
  assay    = "RNA",
  slot     = "data",
  method   = "mean"
)

keep_cols <- intersect(celltype_levels, colnames(mat_mean))
mat_mean  <- mat_mean[, keep_cols, drop = FALSE]

gene_order <- genes_use[genes_use %in% rownames(mat_mean)]
mat_mean   <- mat_mean[gene_order, , drop = FALSE]

# Row-wise z-score, capped at +-2
mat_plot <- t(scale(t(mat_mean)))
mat_plot[is.na(mat_plot)] <- 0
mat_plot <- pmax(pmin(mat_plot, 2.0), -2.0)

heat_cols <- colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(101)
heat_bk   <- seq(-2.0, 2.0, length.out = 101)

# ggplot2 heatmap (italic gene names, Q2/Q4 divider)
n_q4      <- sum(genes_q4 %in% rownames(mat_plot))
divider_y <- n_q4 + 0.5

df_heat <- mat_plot %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "celltype", values_to = "zscore") %>%
  dplyr::mutate(
    gene     = factor(gene, levels = rev(rownames(mat_plot))),
    celltype = factor(celltype, levels = celltype_levels),
    quad     = ifelse(gene %in% genes_q2, "Q2", "Q4")
  )

p_heat <- ggplot(df_heat, aes(x = celltype, y = gene, fill = zscore)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colours = c("#2166AC", "#FFFFFF", "#B2182B"),
    limits  = c(-2.0, 2.0),
    oob     = scales::squish,
    name    = "z-score\n(row-wise)"
  ) +
  geom_hline(yintercept = divider_y, color = "black", linewidth = 0.8) +
  scale_x_discrete(position = "bottom") +
  labs(x = NULL, y = NULL, title = "") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y     = element_text(face = "italic", size = 7, hjust = 1),
    axis.text.x     = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    panel.grid      = element_blank(),
    legend.position = "right",
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 9),
    plot.margin     = margin(t = 5, r = 20, b = 5, l = 5)
  )

ggsave("DEG_heatmap_Q2Q4.png",
       p_heat,
       width  = 7,
       height = max(6, 0.12 * nrow(mat_plot) + 2),
       dpi    = 300,
       bg     = "white")

message("Saved: DEG_heatmap_Q2Q4.png")

# =============================================================================
# SECTION 7: GO Enrichment by group pair (using run_go_for_pair wrapper)
# =============================================================================
message("\n>>> GO per group pair (run_go_for_pair)")

res_poor <- run_go_for_pair(
  ident_d7  = "Poor_d7",
  ident_d1  = "Poor_d1",
  out_png   = file.path(getwd(), "GO_BP_Poor_d7_vs_Poor_d1.png"),
  title_txt = "GO(BP): Poor_d7 vs Poor_d1",
  top_n_up  = 8,
  top_n_down = 8
)

res_good <- run_go_for_pair(
  ident_d7  = "Good_d7",
  ident_d1  = "Good_d1",
  out_png   = file.path(getwd(), "GO_BP_Good_d7_vs_Good_d1.png"),
  title_txt = "GO(BP): Good_d7 vs Good_d1",
  top_n_up  = 8,
  top_n_down = 8
)

message("\n>>> All analyses complete.")