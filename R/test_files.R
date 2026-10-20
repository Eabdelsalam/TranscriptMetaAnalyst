#' Input checks that show an alert (superseded)
#'
#' Superseded by [validate_design()], [validate_counts()] and
#' [validate_match()]; kept until the app uses the results module.
#'
#' @param exp_design,counts Input data frames.
#' @return Called for the alert it shows.
#' @name legacy_input_checks
#' @keywords internal
NULL

#' @rdname legacy_input_checks
#' @export
test_exp_design <- function(exp_design) {
  col_names <- colnames(exp_design)
  
  # test that essential columns existing
  if (!"sample" %in% col_names) {
    shinyalert(text = "The column 'sample'(case sensitive) should exist in the Experimental Design file!")
  }
  
  else if (!"study" %in% col_names) {
    shinyalert(text = "The column 'study'(case sensitive) should exist in the Experimental Design file!")
  }
  
  else if (!"condition" %in% col_names) {
    shinyalert(text = "The column 'condition'(case sensitive) should exist in the Experimental Design file!")
  }
  
  # test the presence of two conditions
  else if (!length(unique(exp_design$condition)) == 2) {
    shinyalert(text = "The Experimental Design file should have two unique conditions only!")
  }
  
  # test if sample names has duplication
  else if (sum(duplicated(exp_design$sample)) > 0) {
    shinyalert(text = "The sample 'column' in the Experimental Design file duplicated sample names, please correct and re-upload it!")
  }
}

#' @rdname legacy_input_checks
#' @export
test_counts <- function(counts) {
  col_names <- colnames(counts)
  
  # test that essential columns existing
  if (!"Geneid" %in% col_names) {
    shinyalert(text = "The column 'Geneid'(case sensitive) should exist in the counts file!")
  }
  
  # test if gene names has duplication
  if (sum(duplicated(counts$Geneid)) > 0) {
    shinyalert(text = "The column 'Geneid' has duplicated gene names, please correct and re-upload it!")
  }
}

#' @rdname legacy_input_checks
#' @export
test_match <- function(counts, exp_design) {
  if (!all(exp_design$sample %in% colnames(counts))) {
    shinyalert(text = "All the samples included in the experimental design 'sample' column should have a correponding counts column with the same sample name in the counts file, and vice versa!")
  }
}