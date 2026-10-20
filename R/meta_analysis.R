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

#' Mean log2 fold change of each gene across studies
#'
#' @param fc_list A list of numeric vectors of log2 fold changes, one per
#'   study, all in the same gene order, or a matrix with one column per study.
#' @return A numeric vector with one mean log2 fold change per gene.
#' @export
#' @examples
#' mean_fold_change(list(A = c(1, -2), B = c(3, NA)))
mean_fold_change <- function(fc_list) {
  fc <- if (is.matrix(fc_list)) fc_list else do.call(cbind, unname(fc_list))
  rowMeans(fc, na.rm = TRUE)
}

#' Meta-significant, direction-consistent DEGs
#'
#' Selects the genes that are significant in the meta-analysis and whose
#' fold changes have the same sign in every study. These are the genes used for
#' the heatmap, the Venn/UpSet plots and the GO/KEGG enrichment.
#'
#' @param fc_select A data frame with gene names as row names, a
#'   `DE.meta` column (1 = meta-significant) and a `signFC` column
#'   (1 or -1 = consistent direction, 0 = conflicting).
#' @return A character vector of gene names.
#' @export
#' @examples
#' fc_select <- data.frame(DE.meta = c(1, 1, 0), signFC = c(1, 0, -1),
#'                         row.names = c("g1", "g2", "g3"))
#' select_meta_degs(fc_select)
select_meta_degs <- function(fc_select) {
  keep <- fc_select$DE.meta == 1 & abs(fc_select$signFC) == 1
  rownames(fc_select)[which(keep)]
}

#' Meta-analysis of per-study p-values
#'
#' Combines the raw p-values of the studies with Fisher's method
#' ([metaRNASeq::fishercomb()]) or the inverse-normal method (as in
#' [metaRNASeq::invnorm()]: each study weighted by the square root of its
#' number of samples), and adjusts the combined p-values with
#' Benjamini-Hochberg. As in DESeq2's independent filtering, a gene whose
#' adjusted p-value is `NA` in a study does not contribute that study's
#' p-value.
#'
#' The combined p-values are computed from the upper tail of the chi-squared
#' (Fisher) or normal (inverse normal) distribution, and the inverse-normal
#' z-scores from the upper tail of the normal distribution, so p-values far
#' below 1e-16 stay accurate instead of being rounded to 0 (or, for the
#' z-scores, to infinity).
#'
#' For the inverse-normal method, p-values of exactly 1 are moved just below 1,
#' so that a study with p = 1 cannot turn the combined statistic into `NaN`.
#'
#' @param pvalues Matrix of raw p-values (genes x studies).
#' @param padj Matrix of adjusted p-values with the same layout.
#' @param method `"fisher"` or `"invnorm"`, see [meta_methods()].
#' @param nrep Number of samples in each study (used by `"invnorm"`).
#' @param alpha Significance level for the meta-analysis.
#' @return A data frame with one row per gene and columns `TestStatistic`,
#'   `rawpval` and `adjpval`.
#' @export
#' @examples
#' p <- matrix(c(0.001, 0.2, 0.5, 0.004, 0.3, 0.9), ncol = 2)
#' run_meta(p, p)
#' run_meta(p, p, method = "invnorm", nrep = c(4, 6))
run_meta <- function(pvalues, padj, method = "fisher", nrep = NULL, alpha = 0.05) {
  method <- match.arg(method, meta_methods())
  pvalues[is.na(padj)] <- NA
  if (method == "invnorm") {
    if (length(nrep) != ncol(pvalues)) {
      stop("`nrep` needs one sample count per study.", call. = FALSE)
    }
    pvalues <- pmin(pvalues, 1 - .Machine$double.eps)
  }
  if (method == "fisher") {
    indpval <- lapply(seq_len(ncol(pvalues)), function(i) pvalues[, i])
    stat <- metaRNASeq::fishercomb(indpval, BHth = alpha)$TestStatistic
    # metaRNASeq returns 1 - pchisq(), which rounds p-values below ~1e-16 to
    # exactly 0; take the upper tail directly to keep them accurate
    rawpval <- stats::pchisq(stat, df = 2 * rowSums(!is.na(pvalues)), lower.tail = FALSE)
  } else {
    stat <- inverse_normal_statistic(pvalues, nrep)
    rawpval <- stats::pnorm(stat, lower.tail = FALSE)
  }
  data.frame(TestStatistic = stat, rawpval = rawpval,
             adjpval = stats::p.adjust(rawpval, method = "BH"))
}

# Inverse-normal combination as in metaRNASeq::invnorm(): each study's
# z-score weighted by sqrt(n_s / sum(n)), over the studies with a p-value.
# The z-scores use qnorm(p, lower.tail = FALSE) instead of qnorm(1 - p), which
# is Inf for p below ~1e-16.
inverse_normal_statistic <- function(pvalues, nrep) {
  z <- stats::qnorm(pvalues, lower.tail = FALSE)
  n <- matrix(nrep, nrow = nrow(z), ncol = ncol(z), byrow = TRUE)
  n[is.na(z)] <- 0
  weighted <- sqrt(n / rowSums(n)) * z
  stat <- rowSums(weighted, na.rm = TRUE)
  stat[rowSums(!is.na(z)) == 0] <- NA
  stat
}

#' Meta-analysis methods
#'
#' @return A named character vector: the method codes used by [run_meta()],
#'   named by their labels.
#' @export
#' @examples
#' meta_methods()
meta_methods <- function() {
  c("Fisher" = "fisher", "Inverse normal" = "invnorm")
}

#' Direction of each gene's fold change across studies
#'
#' @param fold_changes Matrix of log2 fold changes (genes x studies).
#' @return A numeric vector: `1` or `-1` when the gene changes in the same
#'   direction in every study, `0` when the directions conflict, and `NA`
#'   when a fold change is missing.
#' @export
#' @examples
#' consistent_direction(matrix(c(1, -2, 0.5, 3, -1, -0.5), ncol = 2))
consistent_direction <- function(fold_changes) {
  sum_signs <- rowSums(sign(fold_changes))
  ifelse(abs(sum_signs) == ncol(fold_changes), sign(sum_signs), 0)
}

#' Significant DEGs of one study
#'
#' The per-study DEG list offered for download: genes with an adjusted p-value
#' below `alpha` and an absolute log2 fold change above `lfc`.
#'
#' @param result A `tma_result` from [run_tma_pipeline()].
#' @param study Study name.
#' @param lfc Absolute log2 fold-change threshold.
#' @return The rows of the study's DE table that pass both thresholds.
#' @export
#' @examples
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' head(study_significant(res, "A"))
study_significant <- function(result, study, lfc = 1) {
  tab <- result$de[[study]]
  tab[which(tab$padj < result$alpha & abs(tab$log2_fc) > lfc), , drop = FALSE]
}

#' @export
print.tma_result <- function(x, ...) {
  cat("TMA meta-analysis result\n")
  cat(" ", length(x$studies), "studies,", length(x$genes), "genes,",
      nrow(x$design), "samples\n")
  cat("  DE methods:", paste(unique(x$de_methods), collapse = ", "),
      " meta-analysis:", names(meta_methods())[meta_methods() == x$meta_method], "\n")
  cat(" ", length(x$meta_degs), "meta-analysis DEGs (alpha =", x$alpha, ")\n")
  invisible(x)
}
