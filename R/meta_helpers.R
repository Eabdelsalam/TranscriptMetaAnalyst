#' Combine per-study matrices into one gene x sample matrix
#'
#' Binds the per-study matrices column-wise after aligning their rows on gene
#' names, so every value stays attached to its gene and sample.
#'
#' @param mat_list A list of numeric matrices (genes x samples) with gene names
#'   as row names and sample names as column names.
#' @param samples Optional character vector of sample names giving the column
#'   order of the result.
#' @return A numeric matrix with genes as rows and samples as columns.
#' @export
#' @examples
#' a <- matrix(1:4, 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
#' b <- matrix(5:6, 2, dimnames = list(c("g2", "g1"), "s3"))
#' combine_study_matrices(list(A = a, B = b))
combine_study_matrices <- function(mat_list, samples = NULL) {
  genes <- rownames(mat_list[[1]])
  if (is.null(genes)) {
    stop("Each study matrix needs gene names as row names.", call. = FALSE)
  }
  aligned <- lapply(mat_list, function(m) {
    missing_genes <- setdiff(genes, rownames(m))
    if (length(missing_genes) > 0) {
      stop("Study matrices do not share the same genes.", call. = FALSE)
    }
    m[genes, , drop = FALSE]
  })
  combined <- do.call(cbind, unname(aligned))
  if (!is.null(samples)) {
    combined <- combined[, samples, drop = FALSE]
  }
  combined
}

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
#' The number of row clusters must not exceed the number of genes shown, so it
#' depends on the number of genes and not on the number of samples.
#'
#' @param n_genes Number of genes (rows) in the heatmap.
#' @param k Preferred number of row clusters (default 2: up- and
#'   down-regulated genes).
#' @return An integer between 1 and `k`.
#' @export
#' @examples
#' heatmap_row_clusters(100)
#' heatmap_row_clusters(1)
heatmap_row_clusters <- function(n_genes, k = 2L) {
  as.integer(max(1L, min(k, n_genes)))
}

#' Mean log2 fold change of each gene across studies
#'
#' @param fc_list A list of numeric vectors of log2 fold changes, one per
#'   study, all in the same gene order.
#' @return A numeric vector with one mean log2 fold change per gene.
#' @export
#' @examples
#' mean_fold_change(list(A = c(1, -2), B = c(3, NA)))
mean_fold_change <- function(fc_list) {
  rowMeans(do.call(cbind, unname(fc_list)), na.rm = TRUE)
}
