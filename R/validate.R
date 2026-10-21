#' Validate TMA input files
#'
#' Check the experimental design file, the counts file, and whether the two
#' match. Each function returns the problems it found as a character vector,
#' so the app can show them to the user and the checks can be tested without
#' Shiny.
#'
#' @param design Experimental design data frame with `sample`, `study` and
#'   `condition` columns.
#' @param counts Counts data frame with a `Geneid` column and one column per
#'   sample.
#' @return A character vector of error messages; empty when the input is valid.
#' @name validate_inputs
#' @examples
#' demo <- tma_demo_data()
#' validate_design(demo$design)
#' validate_counts(demo$counts)
#' validate_match(demo$counts, demo$design)
#' validate_design(data.frame(sample = "S1"))
NULL

#' @rdname validate_inputs
#' @export
validate_design <- function(design) {
  required <- c("sample", "study", "condition")
  missing_cols <- setdiff(required, colnames(design))
  if (length(missing_cols) > 0) {
    return(sprintf("The column '%s' (case sensitive) should exist in the Experimental Design file!",
                   missing_cols[1]))
  }
  if (length(unique(design$condition)) != 2) {
    return("The Experimental Design file should have two unique conditions only!")
  }
  if (anyDuplicated(design$sample) > 0) {
    return("The 'sample' column in the Experimental Design file has duplicated sample names, please correct and re-upload it!")
  }
  character(0)
}

#' @rdname validate_inputs
#' @export
validate_counts <- function(counts) {
  if (!"Geneid" %in% colnames(counts)) {
    return("The column 'Geneid' (case sensitive) should exist in the counts file!")
  }
  if (anyDuplicated(counts$Geneid) > 0) {
    return("The column 'Geneid' has duplicated gene names, please correct and re-upload it!")
  }
  character(0)
}

#' @rdname validate_inputs
#' @export
validate_match <- function(counts, design) {
  if (!all(design$sample %in% colnames(counts))) {
    return("All the samples included in the experimental design 'sample' column should have a corresponding counts column with the same sample name in the counts file!")
  }
  character(0)
}
