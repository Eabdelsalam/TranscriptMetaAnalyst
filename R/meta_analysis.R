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
#' Selects the genes that are significant in the Fisher meta-analysis and whose
#' fold changes have the same sign in every study. These are the genes used for
#' the heatmap, the Venn/UpSet plots and the GO/KEGG enrichment.
#'
#' @param fc_select A data frame with gene names as row names, a
#'   `DE.fishercomb` column (1 = meta-significant) and a `signFC` column
#'   (1 or -1 = consistent direction, 0 = conflicting).
#' @return A character vector of gene names.
#' @export
#' @examples
#' fc_select <- data.frame(DE.fishercomb = c(1, 1, 0), signFC = c(1, 0, -1),
#'                         row.names = c("g1", "g2", "g3"))
#' select_meta_degs(fc_select)
select_meta_degs <- function(fc_select) {
  keep <- fc_select$DE.fishercomb == 1 & abs(fc_select$signFC) == 1
  rownames(fc_select)[which(keep)]
}

#' Fisher meta-analysis of per-study p-values
#'
#' Combines the raw p-values of the studies with Fisher's method
#' ([metaRNASeq::fishercomb()]) and adjusts them with Benjamini-Hochberg. As in
#' DESeq2's independent filtering, a gene whose adjusted p-value is `NA` in a
#' study does not contribute that study's p-value.
#'
#' @param pvalues Matrix of raw p-values (genes x studies).
#' @param padj Matrix of adjusted p-values with the same layout.
#' @param alpha Significance level for the meta-analysis.
#' @return A data frame with one row per gene and columns `TestStatistic`,
#'   `rawpval` and `adjpval`.
#' @export
#' @examples
#' p <- matrix(c(0.001, 0.2, 0.5, 0.004, 0.3, 0.9), ncol = 2)
#' run_fisher_meta(p, p)
run_fisher_meta <- function(pvalues, padj, alpha = 0.05) {
  pvalues[is.na(padj)] <- NA
  fish <- metaRNASeq::fishercomb(lapply(seq_len(ncol(pvalues)), function(i) pvalues[, i]),
                                 BHth = alpha)
  data.frame(TestStatistic = fish$TestStatistic, rawpval = fish$rawpval, adjpval = fish$adjpval)
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

#' Run the TMA RNA-seq meta-analysis
#'
#' Removes genes with low counts in every study ([filter_low_counts()]), runs
#' DESeq2 in every study, combines the p-values with Fisher's method, removes
#' genes whose direction conflicts between studies, and collects everything
#' the results tabs need in one `tma_result` object.
#'
#' @param counts Counts data frame with a `Geneid` column and one column per
#'   sample.
#' @param design Experimental design data frame with `sample`, `study` and
#'   `condition` columns.
#' @param alpha Significance level for the per-study and meta-analysis tests.
#' @param filter If `TRUE`, remove genes with low counts in every study
#'   before the analysis (see [filter_low_counts()]).
#' @param BPPARAM BiocParallel back-end for the per-study analyses, see
#'   [tma_bpparam()].
#' @return An object of class `tma_result`, a list with
#'   \describe{
#'     \item{design}{the experimental design.}
#'     \item{studies, genes}{study names and the names of the analysed genes.}
#'     \item{de}{per-study DE tables (see [run_deseq2()]).}
#'     \item{norm_counts}{normalised counts of all samples (genes x samples).}
#'     \item{fold_changes}{log2 fold changes (genes x studies).}
#'     \item{meta}{meta-analysis table: `TestStatistic`, `rawpval`, `adjpval`
#'       and `significant` ("yes"/"no"), one row per gene.}
#'     \item{fc_select}{per-study DE flags, the Fisher flag `DE.fishercomb`,
#'       fold changes and the direction `signFC` of each gene.}
#'     \item{meta_degs}{meta-significant, direction-consistent DEGs.}
#'     \item{study_degs}{direction-consistent DEGs of each study.}
#'     \item{idd_irr}{integration-driven discoveries and lost DEGs
#'       ([metaRNASeq::IDD.IRR()]).}
#'     \item{alpha}{the significance level used.}
#'   }
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' length(res$meta_degs)
#' }
run_tma_pipeline <- function(counts, design, alpha = 0.05, filter = TRUE,
                             BPPARAM = tma_bpparam(length(unique(design$study)))) {
  studies <- split_studies(counts, design)
  if (filter) studies <- filter_low_counts(studies)
  de <- run_de_studies(studies, alpha = alpha, BPPARAM = BPPARAM)
  genes <- rownames(studies[[1]]$counts)
  study_ids <- names(studies)
  tables <- lapply(de, `[[`, "table")
  de_column <- function(col) {
    m <- do.call(cbind, lapply(tables, `[[`, col))
    dimnames(m) <- list(genes, study_ids)
    m
  }
  pvalues <- de_column("pvalue")
  padj <- de_column("padj")
  fold_changes <- de_column("log2_fc")

  meta <- run_fisher_meta(pvalues, padj, alpha = alpha)
  sign_fc <- consistent_direction(fold_changes)
  consistent <- abs(sign_fc) == 1 & !is.na(sign_fc)
  study_flags <- ifelse(padj <= alpha, 1, 0)

  fc_select <- data.frame(study_flags,
                          DE.fishercomb = ifelse(meta$adjpval <= alpha, 1, 0),
                          fold_changes,
                          signFC = sign_fc,
                          row.names = genes)
  meta_degs <- select_meta_degs(fc_select)
  study_degs <- lapply(stats::setNames(study_ids, study_ids), function(s) {
    genes[which(consistent & study_flags[, s] == 1)]
  })
  meta$significant <- ifelse(meta$adjpval <= alpha & consistent, "yes", "no")
  rownames(meta) <- genes

  structure(list(
    design = design,
    studies = study_ids,
    genes = genes,
    de = tables,
    norm_counts = combine_study_matrices(lapply(de, `[[`, "norm_counts"), design$sample),
    fold_changes = fold_changes,
    meta = meta,
    fc_select = fc_select,
    meta_degs = meta_degs,
    study_degs = study_degs,
    idd_irr = metaRNASeq::IDD.IRR(meta_degs, unname(study_degs)),
    alpha = alpha
  ), class = "tma_result")
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
#' \donttest{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' head(study_significant(res, "A"))
#' }
study_significant <- function(result, study, lfc = 1) {
  tab <- result$de[[study]]
  tab[which(tab$padj < result$alpha & abs(tab$log2_fc) > lfc), , drop = FALSE]
}

#' @export
print.tma_result <- function(x, ...) {
  cat("TMA meta-analysis result\n")
  cat(" ", length(x$studies), "studies,", length(x$genes), "genes,",
      nrow(x$design), "samples\n")
  cat(" ", length(x$meta_degs), "meta-analysis DEGs (alpha =", x$alpha, ")\n")
  invisible(x)
}
