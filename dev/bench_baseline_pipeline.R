# Baseline ("before") pipeline for the benchmark in dev/benchmark.R.
#
# A non-reactive copy of the computation in inst/TMA/server.R as of the
# fix/critical-bugs branch (before the Phase 2 refactor): DESeq2 per study run
# one after the other, full DESeqDataSet objects kept, Fisher meta-analysis,
# direction filter, and the matrices for the PCA, heatmap and volcano plot.
# It is kept frozen so the "after" numbers can always be compared against it.

baseline_pipeline <- function(counts, design) {
  studies <- unique(design$study)
  genenames <- counts$Geneid

  studies_cts <- lapply(studies, function(x) {
    m <- counts[, colnames(counts) %in% design$sample[design$study == x]]
    rownames(m) <- counts$Geneid
    m
  })
  names(studies_cts) <- studies
  studies_design <- lapply(studies, function(x) {
    d <- design[design$study == x, ]
    rownames(d) <- d$sample
    d
  })
  names(studies_design) <- studies

  studies_deseq <- lapply(studies, function(x) {
    dds <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(
      countData = studies_cts[[x]], colData = studies_design[[x]], design = ~condition
    ))
    DESeq2::DESeq(dds, quiet = TRUE)
  })
  names(studies_deseq) <- studies
  studies_results <- lapply(studies_deseq, function(x) DESeq2::results(x, alpha = 0.05))

  rawpval <- lapply(studies_results, function(x) x$pvalue)
  adjpval <- lapply(studies_results, function(x) x$padj)
  fold_changes <- lapply(studies_results, function(x) x$log2FoldChange)

  de_flags <- as.data.frame(mapply(adjpval, FUN = function(x) ifelse(x <= 0.05, 1, 0)))
  colnames(de_flags) <- studies
  rownames(de_flags) <- genenames
  for (i in studies) {
    rawpval[[i]][is.na(adjpval[[i]])] <- NA
  }
  fishcomb <- metaRNASeq::fishercomb(rawpval, BHth = 0.05)

  signs <- mapply(fold_changes, FUN = sign)
  sum_signs <- apply(signs, 1, sum)
  common_sign <- ifelse(abs(sum_signs) == ncol(signs), sign(sum_signs), 0)
  fc_select <- data.frame(de_flags,
                          DE.fishercomb = ifelse(fishcomb$adjpval <= 0.05, 1, 0),
                          do.call(cbind, fold_changes),
                          signFC = common_sign)
  meta_degs <- rownames(fc_select)[which(fc_select$DE.fishercomb == 1 & abs(fc_select$signFC) == 1)]

  norm_counts <- lapply(studies_deseq, function(x) DESeq2::counts(x, normalized = TRUE))
  list(
    dds = studies_deseq,
    results = studies_results,
    fishcomb = fishcomb,
    fc_select = fc_select,
    meta_degs = meta_degs,
    pca_matrix = combine_study_matrices(norm_counts, design$sample),
    heatmap_matrix = scale_deg_counts(norm_counts, meta_degs, design$sample),
    mean_fc = mean_fold_change(fold_changes)
  )
}
