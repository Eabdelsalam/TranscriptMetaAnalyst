# Runs the meta-analysis on the bundled demo data (once per test run)
demo_cache <- new.env()

demo_result <- function() {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("metaRNASeq")
  if (is.null(demo_cache$result)) {
    demo <- tma_demo_data()
    demo_cache$result <- run_tma_pipeline(demo$counts, demo$design)
  }
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
