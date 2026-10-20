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

#' Differential expression of one study with edgeR
#'
#' Fits `~ condition` with the edgeR quasi-likelihood pipeline (TMM
#' normalisation, `glmQLFit()`, `glmQLFTest()`) and returns the same
#' standardised table as [run_deseq2()].
#'
#' @inheritParams run_deseq2
#' @return A list with `table` and `norm_counts`, as for [run_deseq2()]. `stat`
#'   is the quasi-likelihood F statistic and `lfc_se` is `NA` (edgeR does not
#'   report it).
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' study <- split_studies(demo$counts, demo$design)$A
#' head(run_edger(study$counts, study$design)$table)
#' }
run_edger <- function(counts, design, alpha = 0.05) {
  condition <- factor(design$condition)
  model <- stats::model.matrix(~condition)
  y <- edgeR::DGEList(counts = counts, group = condition)
  y <- edgeR::normLibSizes(y)
  y <- edgeR::estimateDisp(y, model)
  fit <- edgeR::glmQLFit(y, model)
  tab <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
  norm_counts <- size_factor_normalise(counts)
  list(
    table = data.frame(
      gene = rownames(counts),
      base_mean = rowMeans(norm_counts),
      log2_fc = tab$logFC,
      lfc_se = NA_real_,
      stat = tab$F,
      pvalue = tab$PValue,
      padj = tab$FDR
    ),
    norm_counts = norm_counts
  )
}

#' Differential expression of one study with limma-voom
#'
#' Fits `~ condition` with limma after the voom transformation (TMM
#' normalisation, `voom()`, `lmFit()`, `eBayes()`) and returns the same
#' standardised table as [run_deseq2()].
#'
#' @inheritParams run_deseq2
#' @return A list with `table` and `norm_counts`, as for [run_deseq2()]. `stat`
#'   is the moderated t statistic.
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' study <- split_studies(demo$counts, demo$design)$A
#' head(run_limma_voom(study$counts, study$design)$table)
#' }
run_limma_voom <- function(counts, design, alpha = 0.05) {
  condition <- factor(design$condition)
  model <- stats::model.matrix(~condition)
  y <- edgeR::normLibSizes(edgeR::DGEList(counts = counts, group = condition))
  fit <- limma::eBayes(limma::lmFit(limma::voom(y, model), model))
  tab <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  norm_counts <- size_factor_normalise(counts)
  list(
    table = data.frame(
      gene = rownames(counts),
      base_mean = rowMeans(norm_counts),
      log2_fc = tab$logFC,
      lfc_se = unname(fit$stdev.unscaled[, 2] * sqrt(fit$s2.post)),
      stat = tab$t,
      pvalue = tab$P.Value,
      padj = tab$adj.P.Val
    ),
    norm_counts = norm_counts
  )
}

# Median-of-ratios normalised counts (as DESeq2), used for the plots whatever
# the DE method, so the figures are comparable between studies and methods
size_factor_normalise <- function(counts) {
  sweep(counts, 2, DESeq2::estimateSizeFactorsForMatrix(counts), "/")
}

#' Differential expression methods
#'
#' The DE methods TMA offers, and a dispatcher that runs one of them on a
#' study.
#'
#' @return `de_methods()` returns the method names. `run_de()` returns the
#'   result of [run_deseq2()], [run_edger()] or [run_limma_voom()].
#' @name de_methods
#' @examples
#' de_methods()
NULL

#' @rdname de_methods
#' @export
de_methods <- function() {
  c("DESeq2", "edgeR", "limma-voom")
}

#' @rdname de_methods
#' @inheritParams run_deseq2
#' @param method One of [de_methods()].
#' @export
run_de <- function(counts, design, method = "DESeq2", alpha = 0.05) {
  method <- match.arg(method, de_methods())
  fun <- switch(method,
                "DESeq2" = run_deseq2,
                "edgeR" = run_edger,
                "limma-voom" = run_limma_voom)
  fun(counts, design, alpha = alpha)
}

#' Differential expression of every study
#'
#' Runs the chosen DE method on every study. The studies are independent, so
#' they are fitted in parallel with BiocParallel.
#'
#' @param studies Output of [split_studies()].
#' @param methods DE method for all studies (one of [de_methods()]), or a
#'   named character vector with one method per study.
#' @param alpha Significance level passed to the DE method.
#' @param BPPARAM A BiocParallel back-end; by default [tma_bpparam()].
#' @return A named list with one DE result (see [run_deseq2()]) per study.
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' de <- run_de_studies(split_studies(demo$counts, demo$design))
#' names(de)
#' }
run_de_studies <- function(studies, methods = "DESeq2", alpha = 0.05,
                           BPPARAM = tma_bpparam(length(studies))) {
  methods <- study_methods(methods, names(studies))
  jobs <- lapply(names(studies), function(id) c(studies[[id]], method = methods[[id]]))
  res <- BiocParallel::bplapply(jobs, function(s) {
    run_de(s$counts, s$design, method = s$method, alpha = alpha)
  }, BPPARAM = BPPARAM)
  stats::setNames(res, names(studies))
}

# Expands a single DE method, or checks a per-study vector of methods
study_methods <- function(methods, study_ids) {
  bad <- setdiff(methods, de_methods())
  if (length(bad) > 0) {
    stop("Unknown DE method: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (length(methods) == 1 && is.null(names(methods))) {
    return(stats::setNames(rep(methods, length(study_ids)), study_ids))
  }
  missing_ids <- setdiff(study_ids, names(methods))
  if (length(missing_ids) > 0) {
    stop("No DE method given for study: ", paste(missing_ids, collapse = ", "), call. = FALSE)
  }
  methods[study_ids]
}

#' Parallel back-end for the per-study analyses
#'
#' Uses one worker per study, up to [BiocParallel::multicoreWorkers()] (the
#' number of cores minus two, within any limit set for the machine), as forked
#' processes (`MulticoreParam`). On Windows, or when only one worker is
#' available, the studies are run one after the other (`SerialParam`). Set
#' `options(tma.workers = n)` to choose the number of workers.
#'
#' @param n_jobs Number of independent jobs (studies).
#' @return A `BiocParallelParam` object.
#' @export
#' @examples
#' tma_bpparam(2)
tma_bpparam <- function(n_jobs) {
  workers <- min(n_jobs, getOption("tma.workers", BiocParallel::multicoreWorkers()))
  if (.Platform$OS.type == "windows" || workers < 2) {
    BiocParallel::SerialParam()
  } else {
    BiocParallel::MulticoreParam(workers = workers)
  }
}
