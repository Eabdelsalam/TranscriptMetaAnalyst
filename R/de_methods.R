#' Differential expression of one study with DESeq2
#'
#' Fits `~ condition` with DESeq2 and returns a standardised result table and
#' the normalised counts. The full `DESeqDataSet` is not kept, so the result is
#' small enough to hold for many studies.
#'
#' @param counts Count matrix (genes x samples) with gene names as row names.
#' @param design Data frame with a `condition` column and sample names as row
#'   names, in the same order as the columns of `counts`.
#' @param alpha Significance level used by DESeq2's independent filtering.
#' @return A list with
#'   \describe{
#'     \item{table}{data frame with `gene`, `base_mean`, `log2_fc`, `lfc_se`,
#'       `stat`, `pvalue` and `padj`, one row per gene in the order of `counts`.}
#'     \item{norm_counts}{matrix of size-factor normalised counts.}
#'   }
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' study <- split_studies(demo$counts, demo$design)$A
#' de <- run_deseq2(study$counts, study$design)
#' head(de$table)
#' }
run_deseq2 <- function(counts, design, alpha = 0.05) {
  dds <- suppressMessages(suppressWarnings(
    DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = design, design = ~condition)
  ))
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, alpha = alpha)
  list(
    table = data.frame(
      gene = rownames(res),
      base_mean = res$baseMean,
      log2_fc = res$log2FoldChange,
      lfc_se = res$lfcSE,
      stat = res$stat,
      pvalue = res$pvalue,
      padj = res$padj
    ),
    norm_counts = DESeq2::counts(dds, normalized = TRUE)
  )
}

#' Differential expression of every study
#'
#' Runs [run_deseq2()] on every study. The studies are independent, so they
#' are fitted in parallel with BiocParallel.
#'
#' @param studies Output of [split_studies()].
#' @param alpha Significance level passed to [run_deseq2()].
#' @param BPPARAM A BiocParallel back-end; by default [tma_bpparam()].
#' @return A named list with one [run_deseq2()] result per study.
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' de <- run_de_studies(split_studies(demo$counts, demo$design))
#' names(de)
#' }
run_de_studies <- function(studies, alpha = 0.05, BPPARAM = tma_bpparam(length(studies))) {
  BiocParallel::bplapply(studies, function(s) run_deseq2(s$counts, s$design, alpha = alpha),
                         BPPARAM = BPPARAM)
}

#' Parallel back-end for the per-study analyses
#'
#' Uses one worker per study, up to the number of cores minus one, as forked
#' processes (`MulticoreParam`). On Windows, or when only one worker is
#' available, the studies are run one after the other (`SerialParam`). Set
#' `options(tma.workers = n)` to choose the number of workers.
#'
#' @param n_jobs Number of independent jobs (studies).
#' @return A `BiocParallelParam` object.
#' @export
#' @examples
#' tma_bpparam(3)
tma_bpparam <- function(n_jobs) {
  cores <- parallel::detectCores()
  if (is.na(cores)) cores <- 1L
  workers <- getOption("tma.workers", max(1L, cores - 1L))
  workers <- min(n_jobs, workers)
  if (.Platform$OS.type == "windows" || workers < 2) {
    BiocParallel::SerialParam()
  } else {
    BiocParallel::MulticoreParam(workers = workers)
  }
}
