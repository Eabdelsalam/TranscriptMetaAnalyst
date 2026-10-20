#' Path to a bundled example file
#'
#' Returns the path of the example RNA-seq counts or experimental design file
#' shipped with the package (used by the Demo tab and as download templates).
#'
#' @param which Either `"counts"` or `"design"`.
#' @return The full path to the file.
#' @export
#' @examples
#' tma_example_file("design")
tma_example_file <- function(which = c("counts", "design")) {
  which <- match.arg(which)
  system.file("extdata", paste0("example_rnaseq_", which, ".txt"),
              package = "TranscriptMetaAnalyst", mustWork = TRUE)
}

#' Read a tab-delimited TMA input file
#'
#' Reads a counts or experimental design file the same way the app does.
#'
#' @param path Path to a tab-delimited text file with a header line.
#' @return A data frame.
#' @export
#' @examples
#' design <- read_tma_table(tma_example_file("design"))
#' head(design)
read_tma_table <- function(path) {
  utils::read.delim(path, header = TRUE, fill = TRUE, quote = "")
}

#' Load the bundled demo dataset
#'
#' @return A list with `counts` (data frame with a `Geneid` column and one
#'   column per sample) and `design` (data frame with `sample`, `study` and
#'   `condition` columns).
#' @export
#' @examples
#' demo <- tma_demo_data()
#' table(demo$design$study)
tma_demo_data <- function() {
  list(
    counts = read_tma_table(tma_example_file("counts")),
    design = read_tma_table(tma_example_file("design"))
  )
}
