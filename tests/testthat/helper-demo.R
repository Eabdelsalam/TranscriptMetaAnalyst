# Runs the Demo meta-analysis the same way as inst/TMA/server.R (once per test run)
demo_cache <- new.env()

demo_pipeline <- function() {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("metaRNASeq")
  counts_path <- system.file("extdata", "example_rnaseq_counts.txt", package = "TranscriptMetaAnalyst")
  design_path <- system.file("extdata", "example_rnaseq_design.txt", package = "TranscriptMetaAnalyst")
  skip_if(counts_path == "" || design_path == "", "demo data files not available")

  if (!is.null(demo_cache$result)) {
    return(demo_cache$result)
  }

  counts <- utils::read.delim(counts_path)
  design <- utils::read.delim(design_path)
  studies <- unique(design$study)
  genenames <- counts$Geneid

  dds <- lapply(studies, function(x) {
    study_design <- design[design$study == x, ]
    rownames(study_design) <- study_design$sample
    study_counts <- counts[, study_design$sample]
    rownames(study_counts) <- genenames
    dds <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(
      countData = study_counts, colData = study_design, design = ~ condition
    ))
    DESeq2::DESeq(dds, quiet = TRUE)
  })
  names(dds) <- studies
  res <- lapply(dds, DESeq2::results, alpha = 0.05)

  adjpval <- lapply(res, function(x) x$padj)
  rawpval <- lapply(res, function(x) x$pvalue)
  fold_changes <- lapply(res, function(x) x$log2FoldChange)
  for (i in studies) {
    rawpval[[i]][is.na(adjpval[[i]])] <- NA
  }
  fishcomb <- metaRNASeq::fishercomb(rawpval, BHth = 0.05)

  de_flags <- as.data.frame(lapply(adjpval, function(x) ifelse(x <= 0.05, 1, 0)))
  signs <- sapply(fold_changes, sign)
  sum_signs <- rowSums(signs)
  sign_fc <- ifelse(abs(sum_signs) == ncol(signs), sign(sum_signs), 0)
  fc_select <- data.frame(de_flags,
                          DE.fishercomb = ifelse(fishcomb$adjpval <= 0.05, 1, 0),
                          do.call(cbind, fold_changes),
                          signFC = sign_fc,
                          row.names = genenames)
  norm_counts <- lapply(dds, DESeq2::counts, normalized = TRUE)

  demo_cache$result <- list(
    design = design, studies = studies, genenames = genenames, dds = dds,
    fold_changes = fold_changes, fishcomb = fishcomb, fc_select = fc_select,
    norm_counts = norm_counts, meta_degs = select_meta_degs(fc_select)
  )
  demo_cache$result
}

# Reads fixtures/demo-reference.txt: the recorded DEG list and the package
# versions it was produced with
read_demo_reference <- function() {
  lines <- readLines(test_path("fixtures", "demo-reference.txt"))
  meta <- grep("^# [A-Za-z0-9_]+: ", lines, value = TRUE)
  fields <- stats::setNames(sub("^# [^:]+: ", "", meta), sub("^# ([^:]+): .*", "\\1", meta))
  list(
    degs = grep("^#", lines, value = TRUE, invert = TRUE),
    deseq2 = fields[["DESeq2"]],
    metarnaseq = fields[["metaRNASeq"]],
    consistent_genes = as.integer(fields[["consistent_genes"]])
  )
}
