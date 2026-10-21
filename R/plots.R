#' Build the scaled DEG matrix used by the heatmap
#'
#' Subsets each study's normalised counts to the selected genes, z-scores every
#' gene within each study, and combines the studies with
#' [combine_study_matrices()]. A gene with constant counts within a study gets
#' a z-score of 0 in that study instead of `NaN`.
#'
#' @param norm_list A list of normalised count matrices (genes x samples), one
#'   per study.
#' @param genes Character vector of genes to keep (e.g. the meta-analysis DEGs).
#' @param samples Optional character vector giving the column order.
#' @return A numeric matrix with `genes` as rows and samples as columns.
#' @export
#' @examples
#' a <- matrix(c(1, 2, 3, 4, 5, 6), 2, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
#' scale_deg_counts(list(A = a), genes = c("g2", "g1"))
scale_deg_counts <- function(norm_list, genes, samples = NULL) {
  scaled <- lapply(norm_list, function(m) {
    x <- m[genes, , drop = FALSE]
    z <- (x - rowMeans(x)) / apply(x, 1, stats::sd)
    z[!is.finite(z)] <- 0
    z
  })
  combine_study_matrices(scaled, samples)
}

#' Number of row clusters for the DEG heatmap
#'
#' Returns the number of gene (row) clusters to cut the heatmap dendrogram
#' into. The user's choice `k` is used when it is a valid whole number, and it
#' is capped by the number of genes shown, so it depends on the number of
#' genes and not on the number of samples.
#'
#' @param n_genes Number of genes (rows) in the heatmap.
#' @param k Requested number of row clusters, e.g. from the "Number of gene
#'   clusters" input. Missing, non-numeric or values below 1 fall back to the
#'   default of 2 (up- and down-regulated genes).
#' @return An integer between 1 and `n_genes` (at least 1).
#' @export
#' @examples
#' heatmap_row_clusters(100)
#' heatmap_row_clusters(100, k = 5)
#' heatmap_row_clusters(3, k = 10)
#' heatmap_row_clusters(100, k = NA)
heatmap_row_clusters <- function(n_genes, k = 2L) {
  k <- suppressWarnings(as.integer(k))
  if (length(k) != 1 || is.na(k) || k < 1L) {
    k <- 2L
  }
  as.integer(max(1L, min(k, n_genes)))
}

#' Scaled DEG matrix of a TMA result
#'
#' Normalised counts of the meta-analysis DEGs, z-scored within each study
#' (see [scale_deg_counts()]); the data behind the heatmap.
#'
#' @param result A `tma_result` from [run_tma_pipeline()].
#' @param corrected If `TRUE`, use the batch-corrected log2 values (see
#'   [batch_diagnostics]) instead of the normalised counts. By default they
#'   are used when the result was created with `batch_plots = TRUE`.
#' @param genes Genes to show (default all meta-analysis DEGs); see
#'   [heatmap_genes()] to select them by cutoffs.
#' @return A numeric matrix (meta DEGs x samples).
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' dim(heatmap_matrix(res))
#' }
heatmap_matrix <- function(result, corrected = use_corrected(result), genes = result$meta_degs) {
  values <- if (corrected) result$batch$corrected_log else result$norm_counts
  samples_by_study <- split(result$design$sample, result$design$study)[result$studies]
  norm_list <- lapply(samples_by_study, function(s) values[, s, drop = FALSE])
  scale_deg_counts(norm_list, genes, result$design$sample)
}

#' Genes shown in the heatmap
#'
#' Selects meta-analysis DEGs for display: those with an adjusted
#' meta-analysis p-value at or below `padj_cutoff` and an absolute mean log2
#' fold change at or above `lfc_cutoff`, optionally only the `top_n` with the
#' smallest adjusted p-values. This only changes what the heatmap shows; the
#' DEGs themselves are defined by the analysis.
#'
#' @param result A `tma_result` from [run_tma_pipeline()].
#' @param padj_cutoff Largest adjusted meta-analysis p-value to show.
#' @param lfc_cutoff Smallest absolute mean log2 fold change to show.
#' @param top_n Show at most this many genes; `NULL` or `NA` shows all.
#' @return A character vector of genes, in the order of `result$meta_degs`.
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' length(heatmap_genes(res))
#' length(heatmap_genes(res, lfc_cutoff = 2, top_n = 50))
#' }
heatmap_genes <- function(result, padj_cutoff = result$alpha, lfc_cutoff = 0, top_n = NULL) {
  degs <- result$meta_degs
  padj <- result$meta[degs, "adjpval"]
  mean_fc <- abs(mean_fold_change(result$fold_changes[degs, , drop = FALSE]))
  keep <- padj <= padj_cutoff & mean_fc >= lfc_cutoff
  keep[is.na(keep)] <- FALSE
  shown <- degs[keep]
  if (length(top_n) == 1 && !is.na(top_n) && top_n < length(shown)) {
    best <- order(padj[keep])[seq_len(max(0, top_n))]
    shown <- shown[sort(best)]
  }
  shown
}

#' Plots of a TMA result
#'
#' The figures shown in the results tabs. Each takes a `tma_result` (or data
#' derived from it) and returns a plot object, except `plot_pvalue_hist()` and
#' `plot_heatmap_static()`, which draw on the current graphics device.
#'
#' @param result A `tma_result` from [run_tma_pipeline()].
#' @param mat A scaled DEG matrix from [heatmap_matrix()].
#' @param k_row Number of gene (row) clusters, see [heatmap_row_clusters()].
#' @param k_col Number of sample (column) clusters.
#' @param palette Heatmap colours, one of [heatmap_palettes()].
#' @param cluster_rows,cluster_cols Cluster (and reorder) the genes and the
#'   samples. Without clustering, `k_row` and `k_col` are ignored.
#' @param method Agglomeration method for the clustering, passed to
#'   [stats::hclust()].
#' @param show_rownames Show the gene names.
#' @param p_cutoff,fc_cutoff Adjusted p-value and absolute log2 fold-change
#'   thresholds drawn on the volcano plot.
#' @param labels Which genes to label on the volcano plot: `"significant"`
#'   (genes past both thresholds), `"top"` (the `top_n` genes with the
#'   smallest adjusted p-values among those) or `"none"`.
#' @param top_n Number of genes to label with `labels = "top"`.
#' @param y_max Largest -log10 adjusted p-value shown on the volcano plot;
#'   genes above it are drawn at the top edge (and counted in the caption).
#'   `NULL` shows the whole range.
#' @param design Experimental design, used to annotate the static heatmap.
#' @param sets Named list of gene sets for the Venn and UpSet plots.
#' @param volcano A data frame from [volcano_data()].
#' @param go_table A table from [run_go_enrichment()].
#' @param kegg_table A table from [run_kegg_enrichment()].
#' @param gene Gene name.
#' @param nrows Number of facet rows of the counts plot.
#' @return A ggplot, plotly, or ComplexHeatmap object, or `NULL` invisibly for
#'   the functions that draw directly.
#' @name tma_plots
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' plot_pca(res)
#' plot_volcano(volcano_data(res))
#' plot_gene_counts(res, res$meta_degs[1])
#' }
NULL

#' @rdname tma_plots
#' @export
plot_pvalue_hist <- function(result) {
  method <- names(meta_methods())[meta_methods() == result$meta_method]
  graphics::hist(result$meta$rawpval, main = paste(method, "Method p-Values Distribution"),
                 breaks = 100, col = "grey", xlab = "Raw p-values (meta-analysis)")
  invisible(NULL)
}

#' @rdname tma_plots
#' @export
plot_pca <- function(result) {
  design <- result$design
  rownames(design) <- design$sample
  values <- if (use_corrected(result)) result$batch$corrected_log else log2(result$norm_counts + 1)
  se <- SummarizedExperiment::SummarizedExperiment(values, colData = design)
  pca <- suppressMessages(DESeq2::plotPCA(DESeq2::DESeqTransform(se), intgroup = "condition",
                                          returnData = TRUE))
  percent <- round(100 * attr(pca, "percentVar"))
  pca$study <- design$study[match(pca$name, design$sample)]
  # no fixed aspect ratio: DESeq2::plotPCA() scales the axes by the explained
  # variance, which squashes the plot when PC1 dominates
  ggplot2::ggplot(pca, ggplot2::aes(x = PC1, y = PC2, colour = condition, shape = study)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = paste0("PC1: ", percent[1], "% variance"),
                  y = paste0("PC2: ", percent[2], "% variance"))
}

# Whether the PCA and heatmap use the batch-corrected values
use_corrected <- function(result) {
  !is.null(result$batch) && isTRUE(result$batch$plots)
}

#' @rdname tma_plots
#' @export
plot_heatmap <- function(mat, k_row, k_col, palette = "yellow_black_purple",
                         cluster_rows = TRUE, cluster_cols = TRUE, method = "complete",
                         show_rownames = FALSE) {
  dendrogram <- if (cluster_rows && cluster_cols) "both" else if (cluster_rows) "row" else
    if (cluster_cols) "column" else "none"
  suppressWarnings(heatmaply::heatmaply(
    mat,
    k_row = if (cluster_rows) k_row else NA,
    k_col = if (cluster_cols) k_col else NA,
    Rowv = cluster_rows, Colv = cluster_cols, dendrogram = dendrogram,
    hclust_method = method,
    colors = heatmap_palette(palette),
    showticklabels = c(TRUE, show_rownames)
  ))
}

#' @rdname tma_plots
#' @export
plot_heatmap_static <- function(mat, design, k_row, k_col, palette = "yellow_black_purple",
                                cluster_rows = TRUE, cluster_cols = TRUE, method = "complete",
                                show_rownames = FALSE) {
  annot <- design[, c("condition", "study")]
  rownames(annot) <- design$sample
  pheatmap::pheatmap(mat,
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols,
                     clustering_method = method,
                     cutree_rows = if (cluster_rows) k_row else NA,
                     cutree_cols = if (cluster_cols) k_col else NA,
                     show_rownames = show_rownames,
                     color = heatmap_palette(palette),
                     annotation_col = annot)
  invisible(NULL)
}

#' Heatmap colour palettes
#'
#' @return A named character vector of palette codes, named by their labels.
#'   `heatmap_palette()` returns the colours of one palette.
#' @param name A palette code from `heatmap_palettes()`.
#' @param n Number of colours.
#' @export
#' @examples
#' heatmap_palettes()
#' heatmap_palette("blue_white_red", 5)
heatmap_palettes <- function() {
  c("Yellow - black - purple" = "yellow_black_purple",
    "Blue - white - red" = "blue_white_red",
    "Viridis" = "viridis")
}

#' @rdname heatmap_palettes
#' @export
heatmap_palette <- function(name = "yellow_black_purple", n = 150) {
  name <- match.arg(name, heatmap_palettes())
  switch(name,
         yellow_black_purple = grDevices::colorRampPalette(c("yellow", "black", "purple"))(n),
         blue_white_red = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(n),
         viridis = grDevices::hcl.colors(n, "viridis"))
}

#' Gene sets for the Venn and UpSet plots
#'
#' @param result A `tma_result` from [run_tma_pipeline()].
#' @return A named list: the direction-consistent DEGs of each study and the
#'   meta-analysis DEGs (`Meta`).
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' lengths(deg_sets(res))
#' }
deg_sets <- function(result) {
  c(result$study_degs, list(Meta = result$meta_degs))
}

#' @rdname tma_plots
#' @export
plot_venn <- function(sets) {
  ggvenn::ggvenn(sets)
}

#' @rdname tma_plots
#' @export
plot_upset <- function(sets) {
  m <- ComplexHeatmap::make_comb_mat(sets)
  ComplexHeatmap::UpSet(m,
                        top_annotation = ComplexHeatmap::upset_top_annotation(m, add_numbers = TRUE),
                        right_annotation = ComplexHeatmap::upset_right_annotation(m, add_numbers = TRUE))
}

#' Data for the volcano plot
#'
#' @param result A `tma_result` from [run_tma_pipeline()].
#' An adjusted p-value of exactly 0 (possible only when a study's own p-value
#' underflowed) cannot be drawn on a -log10 scale. Such values are set to the
#' smallest non-zero p-value in the data and flagged in `capped`.
#'
#' @return A data frame with `genes`, `fc` (mean log2 fold change across
#'   studies), `pval` (adjusted meta-analysis p-value) and `capped` (`TRUE`
#'   where a p-value of 0 was replaced).
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' head(volcano_data(res))
#' }
volcano_data <- function(result) {
  pval <- result$meta$adjpval
  capped <- !is.na(pval) & pval == 0
  if (any(capped)) pval[capped] <- min(pval[!is.na(pval) & pval > 0])
  data.frame(genes = result$genes,
             fc = mean_fold_change(result$fold_changes),
             pval = pval,
             capped = capped)
}

#' @rdname tma_plots
#' @export
plot_volcano <- function(volcano, p_cutoff = 0.05, fc_cutoff = 1,
                         labels = c("significant", "top", "none"), top_n = 10, y_max = NULL) {
  labels <- match.arg(labels)
  passing <- volcano$pval <= p_cutoff & abs(volcano$fc) >= fc_cutoff
  passing[is.na(passing)] <- FALSE
  select <- switch(labels,
                   significant = NULL,
                   top = volcano$genes[passing][utils::head(order(volcano$pval[passing]), top_n)],
                   none = character(0))
  # genes above the y-axis maximum are drawn at the top edge
  capped <- 0
  if (length(y_max) == 1 && !is.na(y_max) && y_max > 0) {
    above <- !is.na(volcano$pval) & -log10(volcano$pval) > y_max
    capped <- sum(above)
    volcano$pval[above] <- 10^-y_max
  }
  plot <- suppressWarnings(EnhancedVolcano::EnhancedVolcano(
    volcano,
    lab = volcano$genes,
    selectLab = select,
    x = "fc",
    y = "pval",
    title = NULL,
    subtitle = NULL,
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    border = "full",
    caption = if (capped > 0) sprintf("%d gene(s) above the y-axis maximum are drawn at the top", capped) else NULL
  ))
  # label the threshold line so its position is clear on any y scale
  x_right <- max(abs(volcano$fc), na.rm = TRUE)
  plot + ggplot2::annotate("text", x = x_right, y = -log10(p_cutoff), hjust = 1, vjust = -0.4,
                           size = 3.5, label = paste("adj. p =", format(p_cutoff, digits = 3)))
}

#' Which plot shows the overlap of DEG lists
#'
#' Venn diagrams are readable for at most 4 sets (here: the studies plus the
#' meta-analysis); with more sets, only the UpSet plot is shown.
#'
#' @param sets Named list of gene sets, see [deg_sets()].
#' @return `"venn"` (both plots are offered, the Venn diagram first) or
#'   `"upset"` (only the UpSet plot).
#' @export
#' @examples
#' overlap_plot_type(list(A = "g1", B = "g2", Meta = "g1"))
overlap_plot_type <- function(sets) {
  if (length(sets) <= 4) "venn" else "upset"
}

#' @rdname tma_plots
#' @export
plot_go_lollipop <- function(go_table) {
  for_plot <- data.frame(full_term = paste(go_table$term.id, go_table$term.label),
                         fold_enrichment = go_table$fold_enrichment,
                         fdr = go_table$fdr,
                         genes = go_table$number_in_list,
                         plus_minus = go_table$plus_minus,
                         category = go_table$category)
  for_plot <- utils::head(for_plot, 30)
  ggplot2::ggplot(data = for_plot, ggplot2::aes(x = stats::reorder(full_term, fold_enrichment),
                                                y = fold_enrichment)) +
    ggplot2::geom_segment(ggplot2::aes(x = stats::reorder(full_term, fold_enrichment),
                                       xend = stats::reorder(full_term, fold_enrichment),
                                       y = 0,
                                       yend = fold_enrichment,
                                       colour = fdr),
                          linewidth = 1.5) +
    ggplot2::scale_color_gradient2(low = "blue", high = "red", mid = "green",
                                   midpoint = mean(for_plot$fdr),
                                   name = bquote(-log[10](FDR))) +
    ggplot2::geom_point(ggplot2::aes(colour = fdr, size = genes)) +
    ggplot2::labs(x = "GO terms", y = "Fold enrichment", size = "No. of genes") +
    ggplot2::facet_grid(rows = ggplot2::vars(category), scales = "free", space = "free") +
    ggplot2::coord_flip() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(color = "grey20", size = 12),
                   axis.text.y = ggplot2::element_text(color = "grey20", size = 12))
}

#' @rdname tma_plots
#' @export
plot_kegg_lollipop <- function(kegg_table) {
  for_plot <- utils::head(kegg_table, 30)
  for_plot$log_fdr <- log10(for_plot$qvalue) * -1
  ggplot2::ggplot(data = for_plot, ggplot2::aes(x = stats::reorder(Description, FoldEnrichment),
                                                y = FoldEnrichment)) +
    ggplot2::geom_segment(ggplot2::aes(x = stats::reorder(Description, FoldEnrichment),
                                       xend = stats::reorder(Description, FoldEnrichment),
                                       y = 0,
                                       yend = FoldEnrichment,
                                       colour = log_fdr),
                          linewidth = 1.5) +
    ggplot2::scale_color_gradient2(low = "blue", high = "red", mid = "green",
                                   midpoint = mean(for_plot$log_fdr),
                                   name = bquote(-log[10](FDR))) +
    ggplot2::geom_point(ggplot2::aes(colour = log_fdr, size = Count)) +
    ggplot2::labs(x = "KEGG pathways", y = "Fold enrichment", size = "No. of genes") +
    ggplot2::facet_grid(rows = ggplot2::vars(subcategory), scales = "free", space = "free") +
    ggplot2::coord_flip() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(color = "grey20", size = 12),
                   axis.text.y = ggplot2::element_text(color = "grey20", size = 12))
}

#' @rdname tma_plots
#' @export
plot_gene_counts <- function(result, gene, nrows = length(result$studies) %/% 2) {
  design <- result$design
  # normalised counts + 0.5, as DESeq2::plotCounts() returns them
  counts_df <- data.frame(study = design$study,
                          condition = design$condition,
                          count = result$norm_counts[gene, design$sample] + 0.5)
  ggplot2::ggplot(counts_df, ggplot2::aes(x = condition, y = count)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(w = 0.1, h = 0)) +
    ggplot2::labs(title = paste("Counts of", gene, "across studies"),
                  x = "Condition", y = "Normalized counts") +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~study, nrow = max(1, nrows))
}
