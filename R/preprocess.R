#' Split the combined counts table into one dataset per study
#'
#' @param counts Counts data frame with a `Geneid` column and one column per
#'   sample (extra annotation columns are ignored).
#' @param design Experimental design data frame with `sample`, `study` and
#'   `condition` columns.
#' @return A named list with one element per study (in the order the studies
#'   first appear in `design`). Each element is a list with `counts`, a count
#'   matrix (genes x samples, columns in design order), and `design`, the
#'   study's rows of `design` with sample names as row names.
#' @export
#' @examples
#' demo <- tma_demo_data()
#' studies <- split_studies(demo$counts, demo$design)
#' names(studies)
#' dim(studies$A$counts)
split_studies <- function(counts, design) {
  study_ids <- unique(design$study)
  studies <- lapply(study_ids, function(id) {
    study_design <- design[design$study == id, , drop = FALSE]
    rownames(study_design) <- study_design$sample
    study_counts <- as.matrix(counts[, study_design$sample, drop = FALSE])
    rownames(study_counts) <- counts$Geneid
    list(counts = study_counts, design = study_design)
  })
  stats::setNames(studies, study_ids)
}
