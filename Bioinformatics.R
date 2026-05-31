# ============================================================
# Transcriptomic Biological Annotation for Clinical_Intra_Peri Model
# ============================================================
# Project endpoint:
#   Label / Growth_Group: Slow growth vs Fast growth
#   0 / Slow = VDT > 400 days
#   1 / Fast = VDT <= 400 days
#
# Required input files in working directory:
#   1. Clinical_Intra_Peri_Final_Validation_Results.xlsx
#   2. RNA.xlsx
#   3. counts_anno.xlsx
#
# Main analyses:
#   1. Fast vs Slow:
#      - DESeq2 differential expression analysis
#      - GO / KEGG enrichment based on DEGs
#      - Hallmark GSEA
#      - PCA scatter plot without fixed x/y axes
#
#   2. Clinical_Intra_Peri high-risk vs low-risk:
#      - Hallmark GSEA only
#
#   3. Clinical_Intra_Peri continuous LP:
#      - Spearman correlation between gene expression and model LP
#      - Hallmark GSEA
#
#   4. Sankey plot:
#      - Three analysis modules -> significant Hallmark pathways -> leading genes intersecting Fast-vs-Slow DEGs
#
# GSEA significance criterion:
#   qvalue < 0.05
#
# Notes:
#   - This script does NOT perform ssGSEA.
#   - Model high-risk vs low-risk is used for GSEA ranking only; no DEG/GO/KEGG reporting is performed for this comparison.
#   - The high-risk threshold is preferentially extracted from the training cohort Youden threshold in the model results workbook.
# ============================================================


# ============================================================
# 0. Package installation and loading
# ============================================================

cran_packages <- c(
  "readxl", "openxlsx", "dplyr", "tidyr", "stringr", "ggplot2",
  "ggrepel", "pheatmap", "tibble", "purrr", "forcats",
  "ggalluvial", "scales", "patchwork"
)

bioc_packages <- c(
  "DESeq2", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
  "msigdbr", "SummarizedExperiment"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

while (!is.null(dev.list())) {
  dev.off()
}


# ============================================================
# 1. File paths and analysis parameters
# ============================================================

clinical_model_file <- "Clinical_Intra_Peri_Final_Validation_Results.xlsx"
rna_info_file <- "RNA.xlsx"
counts_file <- "counts_anno.xlsx"

final_model <- "Clinical_Intra_Peri"
final_pred_col <- "Clinical_Intra_Peri_Pred"

output_dir <- "Transcriptomic_Annotation_Clinical_Intra_Peri_Hallmark_GSEA_Final"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dir_tables <- file.path(output_dir, "00_Tables")
dir_deg <- file.path(output_dir, "01_Fast_vs_Slow_DEG")
dir_go_kegg <- file.path(output_dir, "02_Fast_vs_Slow_GO_KEGG")
dir_gsea_fast <- file.path(output_dir, "03_Hallmark_GSEA_Fast_vs_Slow")
dir_gsea_model <- file.path(output_dir, "04_Hallmark_GSEA_Model_High_vs_Low")
dir_gsea_lp <- file.path(output_dir, "05_Hallmark_GSEA_Clinical_Intra_Peri_LP")
dir_sankey <- file.path(output_dir, "06_Sankey_Module_Pathway_LeadingGenes")
dir_pca <- file.path(output_dir, "07_PCA_FreeAxis")

for (dd in c(dir_tables, dir_deg, dir_go_kegg, dir_gsea_fast, dir_gsea_model, dir_gsea_lp, dir_sankey, dir_pca)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

deg_log2fc_cutoff <- 1
deg_padj_cutoff <- 0.05

gsea_qvalue_cutoff <- 0.05

top_heatmap_genes <- 50
top_enrich_terms <- 20

top_pca_genes <- 5000

set.seed(1111)


# ============================================================
# 2. General helper functions
# ============================================================

clean_id <- function(x) {
  trimws(as.character(x))
}

clip_prob <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

logit_prob <- function(p) {
  qlogis(clip_prob(p))
}

write_xlsx_safe <- function(x, file) {
  if (is.data.frame(x)) {
    openxlsx::write.xlsx(x, file, overwrite = TRUE)
    return(invisible(TRUE))
  }
  
  if (is.list(x)) {
    old_names <- names(x)
    if (is.null(old_names)) {
      old_names <- paste0("Sheet", seq_along(x))
    }
    
    new_names <- substr(old_names, 1, 31)
    new_names <- make.unique(new_names, sep = "_")
    new_names <- substr(new_names, 1, 31)
    
    if (anyDuplicated(new_names)) {
      new_names <- paste0("Sheet", seq_along(x))
    }
    
    names(x) <- new_names
    openxlsx::write.xlsx(x, file, overwrite = TRUE)
    return(invisible(TRUE))
  }
  
  openxlsx::write.xlsx(x, file, overwrite = TRUE)
  invisible(TRUE)
}

save_pdf_safe <- function(plot_obj, filename, width = 7.2, height = 6.2) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    onefile = TRUE,
    useDingbats = FALSE
  )
  print(plot_obj)
  grDevices::dev.off()
  
  if (!file.exists(filename)) {
    stop("PDF file was not successfully generated: ", filename)
  }
  
  cat("PDF generated: ", normalizePath(filename), "\n")
}

guess_col <- function(df, candidates, required = TRUE, object_name = "data") {
  hit <- candidates[candidates %in% colnames(df)]
  
  if (length(hit) > 0) {
    return(hit[1])
  }
  
  if (required) {
    stop(
      object_name,
      " does not contain any of these columns: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  NA_character_
}

standardize_label <- function(x) {
  x0 <- trimws(as.character(x))
  
  out <- dplyr::case_when(
    x0 %in% c(
      "1", "Fast", "fast", "Fast growth", "fast growth", "快速", "快",
      "VDT<=400", "VDT ≤400 days", "VDT<=400 days", "VDT ≤ 400 days"
    ) ~ "Fast",
    x0 %in% c(
      "0", "Slow", "slow", "Slow growth", "slow growth", "缓慢", "慢",
      "VDT>400", "VDT >400 days", "VDT > 400 days"
    ) ~ "Slow",
    TRUE ~ NA_character_
  )
  
  factor(out, levels = c("Slow", "Fast"))
}

extract_gene_symbol <- function(df) {
  symbol_candidates <- c(
    "SYMBOL", "Symbol", "symbol", "GeneSymbol", "gene_symbol",
    "Gene", "gene", "GeneName", "gene_name", "external_gene_name"
  )
  
  hit <- symbol_candidates[symbol_candidates %in% colnames(df)]
  
  if (length(hit) > 0) {
    return(clean_id(df[[hit[1]]]))
  }
  
  clean_id(df[[1]])
}

collapse_duplicate_genes <- function(count_mat) {
  count_df <- as.data.frame(count_mat)
  count_df$GeneSymbol <- rownames(count_df)
  
  count_df <- count_df %>%
    dplyr::filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::summarise(
      dplyr::across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  out <- as.data.frame(count_df)
  rownames(out) <- out$GeneSymbol
  out$GeneSymbol <- NULL
  
  as.matrix(out)
}

convert_symbol_to_entrez <- function(symbols) {
  symbols <- unique(symbols)
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  
  if (length(symbols) == 0) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  }
  
  suppressMessages(
    bitr_df <- clusterProfiler::bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )
  
  bitr_df %>%
    dplyr::distinct(SYMBOL, ENTREZID, .keep_all = TRUE)
}

make_rank_vector_entrez <- function(res_df, stat_col = "stat") {
  map_df <- convert_symbol_to_entrez(res_df$Gene)
  
  rank_df <- res_df %>%
    dplyr::inner_join(map_df, by = c("Gene" = "SYMBOL")) %>%
    dplyr::filter(!is.na(.data[[stat_col]]), is.finite(.data[[stat_col]])) %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::slice_max(order_by = abs(.data[[stat_col]]), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  gene_list <- rank_df[[stat_col]]
  names(gene_list) <- rank_df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gene_list
}

sanitize_filename <- function(x) {
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9_\\-]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  substr(x, 1, 160)
}

extract_leading_entrez <- function(core_enrichment_string) {
  if (is.na(core_enrichment_string) || core_enrichment_string == "") {
    return(character(0))
  }
  unique(unlist(strsplit(core_enrichment_string, "/")))
}


# ============================================================
# 3. Plot helper functions
# ============================================================

save_volcano <- function(res_df, title, out_pdf, out_png) {
  plot_df <- res_df %>%
    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(padj),
      is.finite(log2FoldChange),
      is.finite(padj)
    ) %>%
    dplyr::mutate(
      padj_plot = pmax(padj, 1e-300),
      Regulation = dplyr::case_when(
        padj < deg_padj_cutoff & log2FoldChange >= deg_log2fc_cutoff ~ "Up",
        padj < deg_padj_cutoff & log2FoldChange <= -deg_log2fc_cutoff ~ "Down",
        TRUE ~ "Not significant"
      ),
      negLog10Padj = -log10(padj_plot)
    )
  
  if (nrow(plot_df) == 0) {
    warning("No valid genes available for volcano plot.")
    return(NULL)
  }
  
  label_df <- plot_df %>%
    dplyr::filter(Regulation != "Not significant") %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 15)
  
  p <- ggplot(plot_df, aes(x = log2FoldChange, y = negLog10Padj, color = Regulation)) +
    geom_point(alpha = 0.75, size = 1.5, na.rm = TRUE) +
    geom_vline(
      xintercept = c(-deg_log2fc_cutoff, deg_log2fc_cutoff),
      linetype = "dashed",
      color = "grey45"
    ) +
    geom_hline(
      yintercept = -log10(deg_padj_cutoff),
      linetype = "dashed",
      color = "grey45"
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = Gene),
      size = 3,
      max.overlaps = 30,
      box.padding = 0.4,
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = c(
        "Up" = "#E64B35FF",
        "Down" = "#4DBBD5FF",
        "Not significant" = "grey70"
      )
    ) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 adjusted P value",
      color = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  save_pdf_safe(p, out_pdf, width = 7.2, height = 6.2)
  
  grDevices::png(out_png, width = 2160, height = 1860, res = 300, type = "cairo")
  print(p)
  grDevices::dev.off()
  
  p
}

save_deg_heatmap <- function(vsd_mat, meta_df, res_df, group_col, title, out_pdf, out_png) {
  sig_genes <- res_df %>%
    dplyr::filter(!is.na(padj), is.finite(padj), padj < deg_padj_cutoff) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = top_heatmap_genes) %>%
    dplyr::pull(Gene)
  
  sig_genes <- intersect(sig_genes, rownames(vsd_mat))
  
  if (length(sig_genes) < 2) {
    warning("Fewer than two significant genes; heatmap will not be generated.")
    return(NULL)
  }
  
  mat <- vsd_mat[sig_genes, , drop = FALSE]
  mat <- t(scale(t(mat)))
  mat[!is.finite(mat)] <- 0
  
  ann_col <- meta_df %>%
    dplyr::select(all_of(group_col)) %>%
    as.data.frame()
  
  rownames(ann_col) <- rownames(meta_df)
  ann_col <- ann_col[colnames(mat), , drop = FALSE]
  
  grDevices::pdf(out_pdf, width = 8.2, height = 8.8, useDingbats = FALSE)
  pheatmap::pheatmap(
    mat,
    annotation_col = ann_col,
    show_colnames = FALSE,
    fontsize_row = 7,
    main = title,
    color = colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(100)
  )
  grDevices::dev.off()
  
  grDevices::png(out_png, width = 2200, height = 2400, res = 300, type = "cairo")
  pheatmap::pheatmap(
    mat,
    annotation_col = ann_col,
    show_colnames = FALSE,
    fontsize_row = 7,
    main = title,
    color = colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(100)
  )
  grDevices::dev.off()
  
  invisible(TRUE)
}

save_enrich_dotplot <- function(enrich_obj, title, out_pdf, out_png, show_n = 20) {
  if (is.null(enrich_obj)) return(NULL)
  if (nrow(as.data.frame(enrich_obj)) == 0) return(NULL)
  
  p <- enrichplot::dotplot(enrich_obj, showCategory = show_n) +
    ggtitle(title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  save_pdf_safe(p, out_pdf, width = 8.5, height = 6.8)
  
  grDevices::png(out_png, width = 2550, height = 2040, res = 300, type = "cairo")
  print(p)
  grDevices::dev.off()
  
  p
}


# ============================================================
# 4. Custom Hallmark GSEA curve functions
# ============================================================

compute_gsea_curve <- function(gene_list, gene_set_entrez, p = 1) {
  gene_list <- sort(gene_list, decreasing = TRUE)
  genes <- names(gene_list)
  N <- length(gene_list)
  
  hits <- genes %in% gene_set_entrez
  Nh <- sum(hits)
  
  if (Nh == 0) {
    stop("No overlap between ranked gene list and gene set.")
  }
  
  weights <- abs(gene_list)^p
  hit_weights <- ifelse(hits, weights, 0)
  norm_hit <- sum(hit_weights)
  
  running_hit <- cumsum(hit_weights / norm_hit)
  running_miss <- cumsum(ifelse(!hits, 1 / (N - Nh), 0))
  running_es <- running_hit - running_miss
  
  data.frame(
    Index = seq_len(N),
    Gene = genes,
    RankMetric = as.numeric(gene_list),
    Hit = hits,
    RunningES = running_es
  )
}

plot_one_gsea_curve <- function(gene_list,
                                pathway_entrez,
                                pathway_name,
                                pathway_id,
                                NES,
                                pvalue,
                                FDR,
                                module_name,
                                A_label,
                                B_label,
                                out_pdf,
                                out_png) {
  curve_df <- compute_gsea_curve(gene_list, pathway_entrez)
  
  if (NES >= 0) {
    peak_idx <- which.max(curve_df$RunningES)
  } else {
    peak_idx <- which.min(curve_df$RunningES)
  }
  
  peak_es <- curve_df$RunningES[peak_idx]
  zero_idx <- which.min(abs(curve_df$RankMetric))
  
  p_text <- ifelse(
    is.na(pvalue),
    "p=NA",
    ifelse(pvalue < 0.001, "p<0.001", paste0("p=", sprintf("%.3f", pvalue)))
  )
  
  fdr_text <- ifelse(
    is.na(FDR),
    "FDR: NA",
    ifelse(FDR < 0.001, "FDR<0.001", paste0("FDR: ", sprintf("%.3f", FDR)))
  )
  
  stat_text <- paste0("NES: ", sprintf("%.2f", NES), ", ", p_text, ", ", fdr_text)
  title_text <- paste0(pathway_name, "(", pathway_id, ")")
  
  hit_df <- curve_df %>% dplyr::filter(Hit)
  
  y_range <- range(curve_df$RunningES, na.rm = TRUE)
  y_span <- diff(y_range)
  if (y_span == 0) y_span <- 1
  
  stat_x <- round(nrow(curve_df) * 0.60)
  stat_y <- y_range[2] - 0.10 * y_span
  
  p1 <- ggplot(curve_df, aes(x = Index, y = RunningES)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.6) +
    geom_line(color = "#00A087FF", linewidth = 1.1) +
    geom_vline(xintercept = peak_idx, linetype = "dashed", color = "#00A087FF", linewidth = 0.7) +
    annotate(
      "text",
      x = peak_idx,
      y = peak_es,
      label = paste0("ES: ", sprintf("%.3f", peak_es)),
      hjust = ifelse(peak_idx < nrow(curve_df) / 2, -0.05, 1.05),
      vjust = ifelse(peak_es >= 0, -0.5, 1.2),
      size = 3.4,
      color = "black"
    ) +
    annotate(
      "text",
      x = stat_x,
      y = stat_y,
      label = stat_text,
      hjust = 0,
      size = 3.6,
      color = "black"
    ) +
    labs(title = title_text, x = NULL, y = "Enrichment score (ES)") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold")
    )
  
  p2 <- ggplot(hit_df, aes(x = Index)) +
    geom_segment(aes(xend = Index, y = 0, yend = 1), color = "black", linewidth = 0.35) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line.y = element_blank(),
      axis.title = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(0, 5, 0, 5)
    )
  
  p3 <- ggplot(curve_df, aes(x = Index, y = 1, fill = RankMetric)) +
    geom_tile(height = 1) +
    annotate("text", x = 1, y = 1.65, label = A_label, hjust = 0, size = 3.6, fontface = "bold") +
    annotate("text", x = nrow(curve_df), y = 1.65, label = B_label, hjust = 1, size = 3.6, fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, guide = "none") +
    scale_y_continuous(limits = c(0.5, 1.8), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(plot.margin = margin(0, 5, 0, 5))
  
  p4 <- ggplot(curve_df, aes(x = Index, y = RankMetric)) +
    geom_area(fill = "grey75", alpha = 0.75) +
    geom_line(color = "grey35", linewidth = 0.45) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.45) +
    geom_vline(xintercept = zero_idx, linetype = "dashed", color = "grey45", linewidth = 0.6) +
    annotate(
      "text",
      x = zero_idx,
      y = max(curve_df$RankMetric, na.rm = TRUE) * 0.85,
      label = paste0("Zero score at ", zero_idx),
      hjust = ifelse(zero_idx < nrow(curve_df) / 2, -0.05, 1.05),
      size = 3.4,
      color = "black"
    ) +
    labs(x = "Rank in ordered gene list", y = "Ranked Metric (signal_to_noise)") +
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold")
    )
  
  final_plot <- p1 / p2 / p3 / p4 +
    patchwork::plot_layout(heights = c(3.2, 0.45, 0.35, 1.7))
  
  save_pdf_safe(final_plot, out_pdf, width = 7.6, height = 7.8)
  
  grDevices::png(out_png, width = 2280, height = 2340, res = 300, type = "cairo")
  print(final_plot)
  grDevices::dev.off()
  
  final_plot
}

plot_significant_hallmark_gsea_curves <- function(gsea_obj,
                                                  gene_list,
                                                  hallmark_term2gene,
                                                  module_name,
                                                  A_label,
                                                  B_label,
                                                  out_dir,
                                                  q_cutoff = 0.05) {
  if (is.null(gsea_obj)) return(data.frame())
  
  gsea_df <- as.data.frame(gsea_obj)
  if (nrow(gsea_df) == 0) return(data.frame())
  
  sig_df <- gsea_df %>%
    dplyr::filter(!is.na(qvalue), qvalue < q_cutoff) %>%
    dplyr::arrange(qvalue)
  
  if (nrow(sig_df) == 0) {
    write_xlsx_safe(
      data.frame(Message = paste0("No significant Hallmark GSEA pathways with qvalue < ", q_cutoff)),
      file.path(out_dir, paste0(module_name, "_No_Significant_Hallmark_GSEA.xlsx"))
    )
    return(sig_df)
  }
  
  write_xlsx_safe(
    sig_df,
    file.path(out_dir, paste0(module_name, "_Significant_Hallmark_GSEA_qvalue_lt_0.05.xlsx"))
  )
  
  for (i in seq_len(nrow(sig_df))) {
    pathway_id <- sig_df$ID[i]
    pathway_name <- sig_df$Description[i]
    
    pathway_entrez <- hallmark_term2gene %>%
      dplyr::filter(gs_name == pathway_id) %>%
      dplyr::pull(entrez_gene) %>%
      as.character() %>%
      unique()
    
    file_base <- sanitize_filename(paste0(module_name, "_", pathway_id))
    
    tryCatch(
      {
        plot_one_gsea_curve(
          gene_list = gene_list,
          pathway_entrez = pathway_entrez,
          pathway_name = pathway_name,
          pathway_id = pathway_id,
          NES = sig_df$NES[i],
          pvalue = sig_df$pvalue[i],
          FDR = sig_df$qvalue[i],
          module_name = module_name,
          A_label = A_label,
          B_label = B_label,
          out_pdf = file.path(out_dir, paste0(file_base, "_GSEA_curve.pdf")),
          out_png = file.path(out_dir, paste0(file_base, "_GSEA_curve.png"))
        )
      },
      error = function(e) {
        warning("GSEA curve plotting failed for ", pathway_id, ": ", e$message)
      }
    )
  }
  
  sig_df
}


# ============================================================
# 5. Read input files and build count matrix
# ============================================================

if (!file.exists(clinical_model_file)) stop("Cannot find file: ", clinical_model_file)
if (!file.exists(rna_info_file)) stop("Cannot find file: ", rna_info_file)
if (!file.exists(counts_file)) stop("Cannot find file: ", counts_file)

model_sheets <- readxl::excel_sheets(clinical_model_file)

score_sheet <- dplyr::case_when(
  "final_score_prediction" %in% model_sheets ~ "final_score_prediction",
  "score_prediction_data" %in% model_sheets ~ "score_prediction_data",
  TRUE ~ model_sheets[1]
)

metrics_sheet <- dplyr::case_when(
  "final_metrics" %in% model_sheets ~ "final_metrics",
  "metrics_7models_train_test" %in% model_sheets ~ "metrics_7models_train_test",
  TRUE ~ NA_character_
)

model_score_df <- readxl::read_excel(clinical_model_file, sheet = score_sheet) %>%
  as.data.frame()

colnames(model_score_df) <- trimws(colnames(model_score_df))

if (!final_pred_col %in% colnames(model_score_df)) {
  pred_candidates <- grep("Clinical.*Intra.*Peri.*Pred|Clinical_Intra_Peri", colnames(model_score_df), value = TRUE)
  if (length(pred_candidates) == 0) {
    stop("No Clinical_Intra_Peri prediction column was found in the model results workbook.")
  }
  final_pred_col <- pred_candidates[1]
}

model_patient_col <- guess_col(
  model_score_df,
  c("PatientID", "patientID", "Patient_Id", "ID", "id"),
  required = TRUE,
  object_name = "Clinical model workbook"
)

model_score_df <- model_score_df %>%
  dplyr::mutate(
    PatientID = clean_id(.data[[model_patient_col]]),
    Clinical_Intra_Peri_Pred_for_RNA = clip_prob(.data[[final_pred_col]]),
    Clinical_Intra_Peri_LP = logit_prob(Clinical_Intra_Peri_Pred_for_RNA)
  )

final_threshold <- NA_real_

if (!is.na(metrics_sheet)) {
  metrics_df_tmp <- readxl::read_excel(clinical_model_file, sheet = metrics_sheet) %>%
    as.data.frame()
  
  colnames(metrics_df_tmp) <- trimws(colnames(metrics_df_tmp))
  
  if (all(c("Model", "Dataset", "Threshold") %in% colnames(metrics_df_tmp))) {
    thr <- metrics_df_tmp %>%
      dplyr::filter(Model == final_model, Dataset %in% c("Train", "Training")) %>%
      dplyr::pull(Threshold)
    
    if (length(thr) > 0) {
      final_threshold <- as.numeric(thr[1])
    }
  }
}

rna_info <- readxl::read_excel(rna_info_file, guess_max = 100000) %>%
  as.data.frame()

colnames(rna_info) <- trimws(colnames(rna_info))

rna_patient_col <- guess_col(
  rna_info,
  c("PatientID", "patientID", "Patient_Id", "ID", "id"),
  required = TRUE,
  object_name = "RNA.xlsx"
)

rna_sample_col <- guess_col(
  rna_info,
  c("SampleID", "sampleID", "Sample_Id", "sample_id", "RNA_ID", "RNAID"),
  required = FALSE,
  object_name = "RNA.xlsx"
)

rna_label_col <- guess_col(
  rna_info,
  c("Label", "label", "GrowthLabel", "Growth_Group", "GrowthGroup", "VDT_group", "VDT_Group"),
  required = TRUE,
  object_name = "RNA.xlsx"
)

rna_info <- rna_info %>%
  dplyr::mutate(
    PatientID = clean_id(.data[[rna_patient_col]]),
    SampleID = if (!is.na(rna_sample_col)) clean_id(.data[[rna_sample_col]]) else clean_id(.data[[rna_patient_col]]),
    Growth_Group = standardize_label(.data[[rna_label_col]])
  )

if (any(is.na(rna_info$Growth_Group))) {
  stop("Some labels in RNA.xlsx cannot be recognized as Fast or Slow. Please check Label coding.")
}

counts_raw <- readxl::read_excel(counts_file, guess_max = 100000) %>%
  as.data.frame()

colnames(counts_raw) <- trimws(colnames(counts_raw))

gene_symbol <- extract_gene_symbol(counts_raw)

sample_ids <- unique(rna_info$SampleID)
patient_ids <- unique(rna_info$PatientID)

count_cols_by_sample <- intersect(colnames(counts_raw), sample_ids)
count_cols_by_patient <- intersect(colnames(counts_raw), patient_ids)

if (length(count_cols_by_sample) > 0) {
  count_cols <- count_cols_by_sample
  match_mode <- "SampleID"
} else if (length(count_cols_by_patient) > 0) {
  count_cols <- count_cols_by_patient
  match_mode <- "PatientID"
} else {
  stop("Sample columns in counts_anno.xlsx cannot be matched with either SampleID or PatientID in RNA.xlsx.")
}

message("Counts columns matched by: ", match_mode)
message("Matched RNA-seq samples before model matching: ", length(count_cols))

count_mat <- counts_raw[, count_cols, drop = FALSE]

for (cc in colnames(count_mat)) {
  count_mat[[cc]] <- suppressWarnings(as.numeric(count_mat[[cc]]))
}

count_mat <- as.matrix(count_mat)
rownames(count_mat) <- gene_symbol
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)
count_mat <- collapse_duplicate_genes(count_mat)

keep_genes <- rowSums(count_mat >= 10) >= max(3, floor(0.10 * ncol(count_mat)))
count_mat <- count_mat[keep_genes, , drop = FALSE]


# ============================================================
# 6. Match RNA metadata with model predictions
# ============================================================

if (match_mode == "SampleID") {
  meta_df <- rna_info %>%
    dplyr::filter(SampleID %in% colnames(count_mat)) %>%
    dplyr::distinct(SampleID, .keep_all = TRUE)
  
  meta_df <- meta_df[match(colnames(count_mat), meta_df$SampleID), , drop = FALSE]
  meta_df$Matched_Count_Column <- meta_df$SampleID
  
} else {
  meta_df <- rna_info %>%
    dplyr::filter(PatientID %in% colnames(count_mat)) %>%
    dplyr::distinct(PatientID, .keep_all = TRUE)
  
  meta_df <- meta_df[match(colnames(count_mat), meta_df$PatientID), , drop = FALSE]
  meta_df$Matched_Count_Column <- meta_df$PatientID
}

if (any(is.na(meta_df$Matched_Count_Column))) {
  stop("Some samples in metadata cannot be matched with count matrix columns.")
}

meta_df <- meta_df %>%
  dplyr::left_join(
    model_score_df %>%
      dplyr::select(
        PatientID,
        Clinical_Intra_Peri_Pred_for_RNA,
        Clinical_Intra_Peri_LP,
        dplyr::any_of(c("Clinical_Score", "Intra_Radscore", "Peri_Radscore"))
      ),
    by = "PatientID"
  )

rownames(meta_df) <- meta_df$Matched_Count_Column
meta_df <- meta_df[colnames(count_mat), , drop = FALSE]

keep_samples <- !is.na(meta_df$Clinical_Intra_Peri_Pred_for_RNA)

if (sum(keep_samples) < 6) {
  stop(
    "Too few RNA-seq samples were matched with Clinical_Intra_Peri predictions. Matched sample count = ",
    sum(keep_samples),
    ". Please check PatientID consistency."
  )
}

meta_df <- meta_df[keep_samples, , drop = FALSE]
matched_cols <- meta_df$Matched_Count_Column

missing_cols <- setdiff(matched_cols, colnames(count_mat))
if (length(missing_cols) > 0) {
  stop("These Matched_Count_Column values are not in count_mat: ", paste(missing_cols, collapse = ", "))
}

count_mat <- count_mat[, matched_cols, drop = FALSE]
rownames(meta_df) <- matched_cols

if (!identical(colnames(count_mat), rownames(meta_df))) {
  stop("count_mat column names are not identical to meta_df row names.")
}

if (!is.finite(final_threshold)) {
  final_threshold <- median(meta_df$Clinical_Intra_Peri_Pred_for_RNA, na.rm = TRUE)
  warning(
    "Training Youden threshold was not found in the model workbook. RNA cohort median predicted probability is used as fallback threshold: ",
    final_threshold
  )
}

meta_df <- meta_df %>%
  dplyr::mutate(
    Model_Risk_Group = ifelse(
      Clinical_Intra_Peri_Pred_for_RNA >= final_threshold,
      "Model-high-risk",
      "Model-low-risk"
    ),
    Model_Risk_Group = factor(Model_Risk_Group, levels = c("Model-low-risk", "Model-high-risk")),
    Growth_Group = factor(Growth_Group, levels = c("Slow", "Fast"))
  )

write_xlsx_safe(
  list(
    RNA_metadata = meta_df %>% tibble::rownames_to_column("Count_Matrix_Column"),
    matrix_info = data.frame(
      Genes = nrow(count_mat),
      Samples = ncol(count_mat),
      Match_Mode = match_mode,
      Model_Risk_Threshold = final_threshold
    ),
    count_columns = data.frame(Count_Column = colnames(count_mat))
  ),
  file.path(dir_tables, "RNA_Metadata_Matched_With_Model_Prediction.xlsx")
)

cat("\nRNA samples by Growth Group:\n")
print(table(meta_df$Growth_Group))

cat("\nRNA samples by Model Risk Group:\n")
print(table(meta_df$Model_Risk_Group))


# ============================================================
# 7. Hallmark gene sets
# ============================================================

hallmark_msig <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")

hallmark_term2gene <- hallmark_msig %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::mutate(entrez_gene = as.character(entrez_gene)) %>%
  dplyr::distinct()

hallmark_symbol_to_entrez <- hallmark_msig %>%
  dplyr::select(gs_name, gene_symbol, entrez_gene) %>%
  dplyr::mutate(entrez_gene = as.character(entrez_gene)) %>%
  dplyr::distinct()


# ============================================================
# 8. Fast vs Slow: DESeq2, GO/KEGG, Hallmark GSEA
# ============================================================

dds_growth <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta_df,
  design = ~ Growth_Group
)

dds_growth <- dds_growth[rowSums(DESeq2::counts(dds_growth)) > 10, ]
dds_growth <- DESeq2::DESeq(dds_growth)

res_growth <- DESeq2::results(
  dds_growth,
  contrast = c("Growth_Group", "Fast", "Slow")
)

res_growth_df <- as.data.frame(res_growth) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(
    Regulation = dplyr::case_when(
      !is.na(padj) & padj < deg_padj_cutoff & log2FoldChange >= deg_log2fc_cutoff ~ "Up_in_Fast",
      !is.na(padj) & padj < deg_padj_cutoff & log2FoldChange <= -deg_log2fc_cutoff ~ "Down_in_Fast",
      TRUE ~ "Not_significant"
    )
  )

deg_growth_sig <- res_growth_df %>%
  dplyr::filter(
    !is.na(padj),
    is.finite(padj),
    padj < deg_padj_cutoff,
    abs(log2FoldChange) >= deg_log2fc_cutoff
  )

write_xlsx_safe(
  list(DEG_all = res_growth_df, DEG_sig = deg_growth_sig),
  file.path(dir_deg, "DESeq2_Fast_vs_Slow_Results.xlsx")
)

save_volcano(
  res_growth_df,
  "Fast growth vs Slow growth",
  file.path(dir_deg, "Volcano_Fast_vs_Slow.pdf"),
  file.path(dir_deg, "Volcano_Fast_vs_Slow.png")
)

vsd_growth <- DESeq2::vst(dds_growth, blind = FALSE)
vsd_mat <- SummarizedExperiment::assay(vsd_growth)

save_deg_heatmap(
  vsd_mat = vsd_mat,
  meta_df = meta_df,
  res_df = res_growth_df,
  group_col = "Growth_Group",
  title = "Top DEGs: Fast vs Slow",
  out_pdf = file.path(dir_deg, "Heatmap_Top_DEGs_Fast_vs_Slow.pdf"),
  out_png = file.path(dir_deg, "Heatmap_Top_DEGs_Fast_vs_Slow.png")
)

deg_symbols <- deg_growth_sig$Gene
deg_map <- convert_symbol_to_entrez(deg_symbols)
deg_entrez <- unique(deg_map$ENTREZID)

ego_bp <- NULL
ego_cc <- NULL
ego_mf <- NULL
ekegg <- NULL

if (length(deg_entrez) >= 5) {
  ego_bp <- clusterProfiler::enrichGO(
    gene = deg_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.25,
    readable = TRUE
  )
  
  ego_cc <- clusterProfiler::enrichGO(
    gene = deg_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.25,
    readable = TRUE
  )
  
  ego_mf <- clusterProfiler::enrichGO(
    gene = deg_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.25,
    readable = TRUE
  )
  
  ekegg <- clusterProfiler::enrichKEGG(
    gene = deg_entrez,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.25
  )
}

write_xlsx_safe(
  list(
    GO_BP = if (!is.null(ego_bp)) as.data.frame(ego_bp) else data.frame(),
    GO_CC = if (!is.null(ego_cc)) as.data.frame(ego_cc) else data.frame(),
    GO_MF = if (!is.null(ego_mf)) as.data.frame(ego_mf) else data.frame(),
    KEGG = if (!is.null(ekegg)) as.data.frame(ekegg) else data.frame(),
    DEG_map = deg_map
  ),
  file.path(dir_go_kegg, "GO_KEGG_Enrichment_Fast_vs_Slow_DEGs.xlsx")
)

save_enrich_dotplot(
  ego_bp,
  "GO Biological Process: Fast vs Slow DEGs",
  file.path(dir_go_kegg, "GO_BP_Dotplot_Fast_vs_Slow.pdf"),
  file.path(dir_go_kegg, "GO_BP_Dotplot_Fast_vs_Slow.png"),
  top_enrich_terms
)

save_enrich_dotplot(
  ekegg,
  "KEGG: Fast vs Slow DEGs",
  file.path(dir_go_kegg, "KEGG_Dotplot_Fast_vs_Slow.pdf"),
  file.path(dir_go_kegg, "KEGG_Dotplot_Fast_vs_Slow.png"),
  top_enrich_terms
)

gene_list_growth <- make_rank_vector_entrez(res_growth_df, stat_col = "stat")

gsea_hallmark_growth <- NULL

if (length(gene_list_growth) >= 50) {
  gsea_hallmark_growth <- clusterProfiler::GSEA(
    geneList = gene_list_growth,
    TERM2GENE = hallmark_term2gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
}

gsea_growth_df <- if (!is.null(gsea_hallmark_growth)) as.data.frame(gsea_hallmark_growth) else data.frame()

sig_gsea_growth <- plot_significant_hallmark_gsea_curves(
  gsea_obj = gsea_hallmark_growth,
  gene_list = gene_list_growth,
  hallmark_term2gene = hallmark_term2gene,
  module_name = "Fast_vs_Slow",
  A_label = "A",
  B_label = "B",
  out_dir = dir_gsea_fast,
  q_cutoff = gsea_qvalue_cutoff
)

write_xlsx_safe(
  list(GSEA_Fast_all = gsea_growth_df, GSEA_Fast_sig = sig_gsea_growth),
  file.path(dir_gsea_fast, "Hallmark_GSEA_Fast_vs_Slow_Results.xlsx")
)


# ============================================================
# 9. Clinical_Intra_Peri high-risk vs low-risk: Hallmark GSEA only
# ============================================================

dds_model <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta_df,
  design = ~ Model_Risk_Group
)

dds_model <- dds_model[rowSums(DESeq2::counts(dds_model)) > 10, ]
dds_model <- DESeq2::DESeq(dds_model)

res_model <- DESeq2::results(
  dds_model,
  contrast = c("Model_Risk_Group", "Model-high-risk", "Model-low-risk")
)

res_model_rank_df <- as.data.frame(res_model) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::arrange(padj)

gene_list_model <- make_rank_vector_entrez(res_model_rank_df, stat_col = "stat")

gsea_hallmark_model <- NULL

if (length(gene_list_model) >= 50) {
  gsea_hallmark_model <- clusterProfiler::GSEA(
    geneList = gene_list_model,
    TERM2GENE = hallmark_term2gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
}

gsea_model_df <- if (!is.null(gsea_hallmark_model)) as.data.frame(gsea_hallmark_model) else data.frame()

sig_gsea_model <- plot_significant_hallmark_gsea_curves(
  gsea_obj = gsea_hallmark_model,
  gene_list = gene_list_model,
  hallmark_term2gene = hallmark_term2gene,
  module_name = "Model_High_vs_Low",
  A_label = "A",
  B_label = "B",
  out_dir = dir_gsea_model,
  q_cutoff = gsea_qvalue_cutoff
)

write_xlsx_safe(
  list(Model_rank = res_model_rank_df, GSEA_Model_all = gsea_model_df, GSEA_Model_sig = sig_gsea_model),
  file.path(dir_gsea_model, "Hallmark_GSEA_Model_High_vs_Low_Results.xlsx")
)


# ============================================================
# 10. Clinical_Intra_Peri continuous LP: Hallmark GSEA
# ============================================================

expr_vst <- SummarizedExperiment::assay(DESeq2::vst(dds_growth, blind = FALSE))
expr_vst <- expr_vst[, rownames(meta_df), drop = FALSE]

lp <- meta_df$Clinical_Intra_Peri_LP
names(lp) <- rownames(meta_df)

cor_gene_df <- lapply(
  rownames(expr_vst),
  function(g) {
    x <- as.numeric(expr_vst[g, ])
    ct <- suppressWarnings(cor.test(x, lp, method = "spearman", exact = FALSE))
    
    data.frame(
      Gene = g,
      Spearman_rho = as.numeric(ct$estimate),
      P_value = ct$p.value,
      stringsAsFactors = FALSE
    )
  }
) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    FDR = p.adjust(P_value, method = "BH"),
    Rank_score = Spearman_rho * -log10(pmax(P_value, 1e-300))
  ) %>%
  dplyr::arrange(desc(Rank_score))

write_xlsx_safe(
  cor_gene_df,
  file.path(dir_gsea_lp, "Gene_Correlation_With_Clinical_Intra_Peri_LP.xlsx")
)

gene_list_lp <- make_rank_vector_entrez(
  cor_gene_df %>% dplyr::rename(stat = Rank_score),
  stat_col = "stat"
)

gsea_hallmark_lp <- NULL

if (length(gene_list_lp) >= 50) {
  gsea_hallmark_lp <- clusterProfiler::GSEA(
    geneList = gene_list_lp,
    TERM2GENE = hallmark_term2gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
}

gsea_lp_df <- if (!is.null(gsea_hallmark_lp)) as.data.frame(gsea_hallmark_lp) else data.frame()

sig_gsea_lp <- plot_significant_hallmark_gsea_curves(
  gsea_obj = gsea_hallmark_lp,
  gene_list = gene_list_lp,
  hallmark_term2gene = hallmark_term2gene,
  module_name = "Clinical_Intra_Peri_LP",
  A_label = "A",
  B_label = "B",
  out_dir = dir_gsea_lp,
  q_cutoff = gsea_qvalue_cutoff
)

write_xlsx_safe(
  list(Gene_Corr_LP = cor_gene_df, GSEA_LP_all = gsea_lp_df, GSEA_LP_sig = sig_gsea_lp),
  file.path(dir_gsea_lp, "Hallmark_GSEA_Clinical_Intra_Peri_LP_Results.xlsx")
)


# ============================================================
# 11. PCA plot with fixed x/y axes
# ============================================================

common_samples <- intersect(colnames(expr_vst), rownames(meta_df))

if (length(common_samples) < 3) {
  stop("Fewer than three samples are available for PCA. Please check sample names.")
}

expr_pca <- expr_vst[, common_samples, drop = FALSE]
meta_pca <- meta_df[common_samples, , drop = FALSE]

gene_var <- apply(expr_pca, 1, var, na.rm = TRUE)
gene_var <- gene_var[is.finite(gene_var)]

top_genes_pca <- names(sort(gene_var, decreasing = TRUE))[1:min(top_pca_genes, length(gene_var))]
expr_pca_top <- expr_pca[top_genes_pca, , drop = FALSE]

pca_res <- prcomp(t(expr_pca_top), center = TRUE, scale. = FALSE)
pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)

pca_df <- as.data.frame(pca_res$x[, 1:2, drop = FALSE]) %>%
  tibble::rownames_to_column("Count_Matrix_Column") %>%
  dplyr::left_join(
    meta_pca %>% tibble::rownames_to_column("Count_Matrix_Column"),
    by = "Count_Matrix_Column"
  ) %>%
  dplyr::mutate(Growth_Group = factor(Growth_Group, levels = c("Slow", "Fast")))

write_xlsx_safe(
  list(
    PCA_coordinates = pca_df,
    PCA_variance = data.frame(
      PC = paste0("PC", seq_along(pca_var)),
      Variance_Explained = pca_var
    )
  ),
  file.path(dir_pca, "PCA_Fast_vs_Slow_FixedAxis_Data.xlsx")
)

p_pca_fixed <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Growth_Group)) +
  geom_point(size = 3.6, alpha = 0.9) +
  scale_color_manual(values = c("Slow" = "#4DBBD5FF", "Fast" = "#E64B35FF")) +
  coord_cartesian(xlim = c(-10, 30), ylim = c(-10, 10)) +
  labs(
    title = "PCA plot of RNA-seq samples",
    x = paste0("PC1 (", sprintf("%.1f", pca_var[1] * 100), "%)"),
    y = paste0("PC2 (", sprintf("%.1f", pca_var[2] * 100), "%)"),
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.position = "right",
    legend.background = element_blank()
  )

save_pdf_safe(
  p_pca_fixed,
  file.path(dir_pca, "PCA_Fast_vs_Slow_FixedAxis.pdf"),
  width = 7.2,
  height = 6.2
)


# ============================================================
# 12. Sankey plot: three modules -> significant Hallmark pathways -> leading genes intersecting DEGs
# ============================================================

build_module_pathway_gene_edges <- function(sig_gsea_df,
                                            module_label,
                                            deg_symbols,
                                            hallmark_symbol_to_entrez) {
  if (is.null(sig_gsea_df) || nrow(sig_gsea_df) == 0) {
    return(data.frame())
  }
  
  edges <- list()
  
  for (i in seq_len(nrow(sig_gsea_df))) {
    pathway_id <- sig_gsea_df$ID[i]
    pathway_desc <- sig_gsea_df$Description[i]
    pathway_nes <- sig_gsea_df$NES[i]
    pathway_qvalue <- sig_gsea_df$qvalue[i]
    
    leading_entrez <- extract_leading_entrez(sig_gsea_df$core_enrichment[i])
    
    leading_symbols <- hallmark_symbol_to_entrez %>%
      dplyr::filter(gs_name == pathway_id, as.character(entrez_gene) %in% as.character(leading_entrez)) %>%
      dplyr::pull(gene_symbol) %>%
      unique()
    
    intersect_genes <- intersect(leading_symbols, deg_symbols)
    
    if (length(intersect_genes) > 0) {
      edges[[length(edges) + 1]] <- data.frame(
        Module = module_label,
        Pathway_ID = pathway_id,
        Pathway = pathway_desc,
        NES = pathway_nes,
        qvalue = pathway_qvalue,
        LeadingGene_DEG_Intersection = intersect_genes,
        stringsAsFactors = FALSE
      )
    }
  }
  
  dplyr::bind_rows(edges)
}

deg_intersection_symbols <- unique(deg_growth_sig$Gene)
deg_intersection_symbols <- deg_intersection_symbols[!is.na(deg_intersection_symbols) & deg_intersection_symbols != ""]

sankey_df <- dplyr::bind_rows(
  build_module_pathway_gene_edges(
    sig_gsea_df = sig_gsea_growth,
    module_label = "Fast vs Slow",
    deg_symbols = deg_intersection_symbols,
    hallmark_symbol_to_entrez = hallmark_symbol_to_entrez
  ),
  build_module_pathway_gene_edges(
    sig_gsea_df = sig_gsea_model,
    module_label = "Clinical_Intra_Peri high-risk vs low-risk",
    deg_symbols = deg_intersection_symbols,
    hallmark_symbol_to_entrez = hallmark_symbol_to_entrez
  ),
  build_module_pathway_gene_edges(
    sig_gsea_df = sig_gsea_lp,
    module_label = "Clinical_Intra_Peri continuous LP",
    deg_symbols = deg_intersection_symbols,
    hallmark_symbol_to_entrez = hallmark_symbol_to_entrez
  )
)

if (nrow(sankey_df) > 0) {
  sankey_df <- sankey_df %>%
    dplyr::left_join(
      res_growth_df %>% dplyr::select(Gene, log2FoldChange, padj, Regulation),
      by = c("LeadingGene_DEG_Intersection" = "Gene")
    ) %>%
    dplyr::mutate(
      Pathway_clean = stringr::str_replace_all(Pathway, "HALLMARK_", ""),
      Pathway_clean = stringr::str_replace_all(Pathway_clean, "_", " "),
      Pathway_clean = stringr::str_to_title(Pathway_clean),
      Gene_clean = LeadingGene_DEG_Intersection,
      DEG_direction = dplyr::case_when(
        Regulation == "Up_in_Fast" ~ "Up in Fast",
        Regulation == "Down_in_Fast" ~ "Down in Fast",
        TRUE ~ "DEG"
      )
    )
  
  write_xlsx_safe(
    sankey_df,
    file.path(dir_sankey, "Sankey_ThreeModules_Pathways_LeadingGenes_DEG_Intersection.xlsx")
  )
  
  p_sankey_module <- ggplot(
    sankey_df,
    aes(axis1 = Module, axis2 = Pathway_clean, axis3 = Gene_clean, y = 1)
  ) +
    ggalluvial::geom_alluvium(
      aes(fill = Module),
      width = 1 / 12,
      alpha = 0.80,
      color = "grey70",
      linewidth = 0.15
    ) +
    ggalluvial::geom_stratum(
      width = 1 / 12,
      fill = "grey96",
      color = "black",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 3,
      color = "black"
    ) +
    scale_x_discrete(
      limits = c("Analysis module", "Significant Hallmark pathway", "Leading gene ∩ DEG"),
      expand = c(0.08, 0.08)
    ) +
    scale_fill_manual(
      values = c(
        "Fast vs Slow" = "#E64B35FF",
        "Clinical_Intra_Peri high-risk vs low-risk" = "#4DBBD5FF",
        "Clinical_Intra_Peri continuous LP" = "#00A087FF"
      ),
      drop = FALSE
    ) +
    labs(
      title = "Significant Hallmark pathways and leading genes intersecting Fast-vs-Slow DEGs",
      x = NULL,
      y = "Gene count",
      fill = "Analysis module"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  save_pdf_safe(
    p_sankey_module,
    file.path(dir_sankey, "Sankey_ThreeModules_Pathways_LeadingGenes_DEG_Intersection.pdf"),
    width = 13.5,
    height = 7.8
  )
  
  grDevices::png(
    file.path(dir_sankey, "Sankey_ThreeModules_Pathways_LeadingGenes_DEG_Intersection.png"),
    width = 3900,
    height = 2250,
    res = 300,
    type = "cairo"
  )
  print(p_sankey_module)
  grDevices::dev.off()
  
} else {
  write_xlsx_safe(
    data.frame(
      Message = "No overlap among significant Hallmark GSEA leading genes and Fast-vs-Slow DEGs.",
      Criterion = paste0("Hallmark GSEA qvalue < ", gsea_qvalue_cutoff)
    ),
    file.path(dir_sankey, "Sankey_No_Overlap_Message.xlsx")
  )
  
  warning("No overlap among significant Hallmark leading genes and Fast-vs-Slow DEGs. Sankey plot was not generated.")
}


# ============================================================
# 13. Integrated summary workbook and RData
# ============================================================

meta_export_df <- meta_df

if ("Count_Matrix_Column" %in% colnames(meta_export_df)) {
  meta_export_df <- meta_export_df %>% dplyr::select(-Count_Matrix_Column)
}

meta_export_df <- meta_export_df %>% tibble::rownames_to_column("Count_Matrix_Column")

summary_list <- list(
  RNA_metadata = meta_export_df,
  DEG_all = res_growth_df,
  DEG_sig = deg_growth_sig,
  GO_BP = if (!is.null(ego_bp)) as.data.frame(ego_bp) else data.frame(),
  GO_CC = if (!is.null(ego_cc)) as.data.frame(ego_cc) else data.frame(),
  GO_MF = if (!is.null(ego_mf)) as.data.frame(ego_mf) else data.frame(),
  KEGG = if (!is.null(ekegg)) as.data.frame(ekegg) else data.frame(),
  GSEA_Fast_all = gsea_growth_df,
  GSEA_Fast_sig = sig_gsea_growth,
  Model_rank = res_model_rank_df,
  GSEA_Model_all = gsea_model_df,
  GSEA_Model_sig = sig_gsea_model,
  Gene_Corr_LP = cor_gene_df,
  GSEA_LP_all = gsea_lp_df,
  GSEA_LP_sig = sig_gsea_lp,
  PCA_coordinates = pca_df,
  Sankey_edges = if (exists("sankey_df")) sankey_df else data.frame()
)

write_xlsx_safe(
  summary_list,
  file.path(output_dir, "Transcriptomic_Annotation_Clinical_Intra_Peri_Hallmark_GSEA_All_Results.xlsx")
)

save(
  count_mat,
  meta_df,
  dds_growth,
  dds_model,
  vsd_mat,
  res_growth_df,
  deg_growth_sig,
  ego_bp,
  ego_cc,
  ego_mf,
  ekegg,
  gsea_hallmark_growth,
  gsea_hallmark_model,
  gsea_hallmark_lp,
  sig_gsea_growth,
  sig_gsea_model,
  sig_gsea_lp,
  cor_gene_df,
  pca_df,
  sankey_df,
  file = file.path(output_dir, "Transcriptomic_Annotation_Clinical_Intra_Peri_Hallmark_GSEA_All_Objects.RData")
)


# ============================================================
# 14. Console summary
# ============================================================

cat("\n================ Transcriptomic biological annotation completed ================\n\n")

cat("Output directory:\n")
cat(normalizePath(output_dir), "\n\n")

cat("Sample matching mode: ", match_mode, "\n")
cat("RNA samples: ", ncol(count_mat), "\n")
cat("Genes retained: ", nrow(count_mat), "\n")
cat("Clinical_Intra_Peri high-risk threshold: ", final_threshold, "\n\n")

cat("Fast vs Slow groups:\n")
print(table(meta_df$Growth_Group))

cat("\nClinical_Intra_Peri high-risk vs low-risk groups:\n")
print(table(meta_df$Model_Risk_Group))

cat("\nFast vs Slow DEG counts:\n")
if (nrow(deg_growth_sig) > 0) {
  print(table(deg_growth_sig$Regulation))
} else {
  cat("0\n")
}

cat("\nSignificant Hallmark GSEA pathways, qvalue < 0.05:\n")
cat("Fast vs Slow: ", nrow(sig_gsea_growth), "\n")
cat("Model high-risk vs low-risk: ", nrow(sig_gsea_model), "\n")
cat("Clinical_Intra_Peri continuous LP: ", nrow(sig_gsea_lp), "\n\n")

cat("Main outputs:\n")
cat(file.path(output_dir, "Transcriptomic_Annotation_Clinical_Intra_Peri_Hallmark_GSEA_All_Results.xlsx"), "\n")
cat(file.path(dir_deg, "Volcano_Fast_vs_Slow.pdf"), "\n")
cat(file.path(dir_deg, "Heatmap_Top_DEGs_Fast_vs_Slow.pdf"), "\n")
cat(file.path(dir_go_kegg, "GO_BP_Dotplot_Fast_vs_Slow.pdf"), "\n")
cat(file.path(dir_go_kegg, "KEGG_Dotplot_Fast_vs_Slow.pdf"), "\n")
cat(file.path(dir_pca, "PCA_Fast_vs_Slow_FreeAxis.pdf"), "\n")
cat(file.path(dir_gsea_fast, "Fast_vs_Slow_*_GSEA_curve.pdf"), "\n")
cat(file.path(dir_gsea_model, "Model_High_vs_Low_*_GSEA_curve.pdf"), "\n")
cat(file.path(dir_gsea_lp, "Clinical_Intra_Peri_LP_*_GSEA_curve.pdf"), "\n")
cat(file.path(dir_sankey, "Sankey_ThreeModules_Pathways_LeadingGenes_DEG_Intersection.pdf"), "\n")

cat("\nInterpretation guide:\n")
cat("1. Fast vs Slow: transcriptomic basis of true rapid growth.\n")
cat("2. Clinical_Intra_Peri high-risk vs low-risk: Hallmark biological annotation of model-defined high-risk phenotype.\n")
cat("3. Clinical_Intra_Peri continuous LP: Hallmark pathways associated with continuous model risk.\n")
cat("4. Sankey: three modules, significant Hallmark pathways, and leading genes intersecting Fast-vs-Slow DEGs.\n")

cat("\n================ Done ================\n")
