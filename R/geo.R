#' Folder where downloaded GEO series are cached
#'
#' @return The path of the cache folder (created if needed), inside the user's
#'   R cache directory ([tools::R_user_dir()]).
#' @export
#' @examples
#' \dontrun{
#' tma_geo_cache()
#' }
tma_geo_cache <- function() {
  dir <- file.path(tools::R_user_dir("TranscriptMetaAnalyst", which = "cache"), "geo")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir
}

#' Download a GEO series
#'
#' Downloads the series matrix of a GEO series (and its platform annotation)
#' with [GEOquery::getGEO()]. Files are cached, so a series is downloaded only
#' once. Needs an internet connection for the first download.
#'
#' @param gse A GEO series accession such as `"GSE82120"`.
#' @param destdir Cache folder.
#' @return A named list of `ExpressionSet` objects, one per platform (named by
#'   the platform accession, e.g. `"GPL15867"`).
#' @export
#' @examples
#' \dontrun{
#' sets <- fetch_geo("GSE82120")
#' names(sets)
#' }
fetch_geo <- function(gse, destdir = tma_geo_cache()) {
  gse <- toupper(trimws(gse))
  if (!grepl("^GSE[0-9]+$", gse)) {
    stop(enrichment_error(sprintf("'%s' is not a GEO series accession (it should look like GSE12345).", gse)))
  }
  sets <- tryCatch(
    suppressMessages(GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = TRUE, destdir = destdir)),
    error = function(e) {
      stop(enrichment_error(sprintf("%s could not be downloaded from GEO (%s). Check the accession and your internet connection.",
                                    gse, conditionMessage(e))))
    }
  )
  stats::setNames(sets, vapply(sets, Biobase::annotation, character(1)))
}

#' Inspect a GEO series before importing it
#'
#' Helpers that suggest how to import a series: which sample annotation column
#' holds the condition, which feature annotation column holds the gene, and
#' whether the series is a two-colour experiment.
#'
#' @param eset An `ExpressionSet` from [fetch_geo()].
#' @return `geo_condition_columns()`: names of the sample columns (title,
#'   source and characteristics) with at least two different values, the most
#'   informative first. `geo_gene_columns()`: names of the feature columns, the
#'   most likely gene column (symbol, ORF, locus tag or gene) first.
#'   `geo_is_two_colour()`: `TRUE` for two-channel arrays (also when the
#'   series lists single channels labelled with both Cy3 and Cy5), whose
#'   channels should not be used as independent samples.
#' @name geo_inspect
#' @examples
#' # a small series standing in for a downloaded one
#' samples <- paste0("GSM", 1:4)
#' es <- Biobase::ExpressionSet(
#'   matrix(c(5, 6, 8, 9, 7, 7, 8, 8), 2, byrow = TRUE,
#'          dimnames = list(c("p1", "p2"), samples)),
#'   phenoData = Biobase::AnnotatedDataFrame(data.frame(
#'     title = paste0("sample", 1:4), channel_count = "1", label_ch1 = "biotin",
#'     characteristics_ch1 = rep(c("co2: high", "co2: low"), each = 2), row.names = samples)),
#'   featureData = Biobase::AnnotatedDataFrame(data.frame(
#'     ID = c("p1", "p2"), ORF = c("sll0001", "sll0002"), row.names = c("p1", "p2"))))
#' geo_condition_columns(es)
#' geo_gene_columns(es)
#' geo_is_two_colour(es)
#'
#' \dontrun{
#' es <- fetch_geo("GSE82120")[[1]]
#' geo_condition_columns(es)
#' geo_gene_columns(es)
#' }
NULL

#' @rdname geo_inspect
#' @export
geo_condition_columns <- function(eset) {
  p <- Biobase::pData(eset)
  cols <- grep("^(characteristics|source_name|title|description)", names(p), value = TRUE)
  n_values <- vapply(cols, function(col) length(unique(p[[col]])), integer(1))
  cols <- cols[n_values >= 2]
  n_values <- n_values[n_values >= 2]
  # characteristics with few distinct values are the best candidates; titles are unique per sample
  cols[order(n_values == nrow(p), !grepl("^characteristics", cols), n_values)]
}

#' @rdname geo_inspect
#' @export
geo_gene_columns <- function(eset) {
  cols <- names(Biobase::fData(eset))
  score <- ifelse(grepl("symbol", cols, ignore.case = TRUE), 1,
                  ifelse(grepl("^orf$|locus|systematic", cols, ignore.case = TRUE), 2,
                         ifelse(grepl("gene", cols, ignore.case = TRUE), 3, 9)))
  cols[order(score)]
}

#' @rdname geo_inspect
#' @export
geo_is_two_colour <- function(eset) {
  p <- Biobase::pData(eset)
  channels <- suppressWarnings(as.integer(p$channel_count))
  # single channels extracted from two-colour arrays carry both dye labels
  dyes <- unique(tolower(as.character(p$label_ch1)))
  any(channels > 1, na.rm = TRUE) || "label_ch2" %in% names(p) ||
    (any(grepl("cy3", dyes)) && any(grepl("cy5", dyes)))
}

#' Turn a GEO series into a study for the meta-analysis
#'
#' Builds, from one platform of a GEO series and the user's choices, the same
#' inputs as an uploaded expression file and design: an expression matrix
#' (probe IDs, genes, samples) and the design rows of the study.
#'
#' @param eset An `ExpressionSet` from [fetch_geo()].
#' @param study Study name (usually the GSE accession).
#' @param condition_column Sample annotation column that holds the condition.
#' @param labels Named character vector mapping the values of
#'   `condition_column` to condition labels. Values mapped to `""` or `NA`, or
#'   not listed, are dropped.
#' @param gene_column Feature annotation column with the gene of each probe, or
#'   `NULL` to keep the probe IDs as genes.
#' @return A list with `matrix` (data frame with `ID`, `Symbol` and one column
#'   per kept sample) and `design` (data frame with `sample`, `study` and
#'   `condition`).
#' @export
#' @examples
#' # a small series standing in for a downloaded one
#' samples <- paste0("GSM", 1:4)
#' es <- Biobase::ExpressionSet(
#'   matrix(c(5, 6, 8, 9, 7, 7, 8, 8), 2, byrow = TRUE,
#'          dimnames = list(c("p1", "p2"), samples)),
#'   phenoData = Biobase::AnnotatedDataFrame(data.frame(
#'     title = paste0("sample", 1:4), channel_count = "1", label_ch1 = "biotin",
#'     characteristics_ch1 = rep(c("co2: high", "co2: low"), each = 2), row.names = samples)),
#'   featureData = Biobase::AnnotatedDataFrame(data.frame(
#'     ID = c("p1", "p2"), ORF = c("sll0001", "sll0002"), row.names = c("p1", "p2"))))
#' study <- geo_study(es, "GSE0", "characteristics_ch1",
#'                    labels = c("co2: high" = "HC", "co2: low" = "LC"), gene_column = "ORF")
#' study$design
#'
#' \dontrun{
#' es <- fetch_geo("GSE63352")[[1]]
#' study <- geo_study(es, "GSE63352", "title",
#'                    labels = c("PCC6803_WT_HC_rep1" = "HC", "PCC6803_WT_HC_rep2" = "HC"),
#'                    gene_column = "ORF")
#' }
geo_study <- function(eset, study, condition_column, labels, gene_column = NULL) {
  p <- Biobase::pData(eset)
  if (!condition_column %in% names(p)) {
    stop("The sample annotation has no column '", condition_column, "'.", call. = FALSE)
  }
  value <- as.character(p[[condition_column]])
  condition <- unname(labels[value])
  keep <- !is.na(condition) & condition != ""
  if (!any(keep)) {
    stop("No sample of ", study, " has a condition label; give a label to at least one value.", call. = FALSE)
  }
  expr <- Biobase::exprs(eset)[, keep, drop = FALSE]
  genes <- if (is.null(gene_column) || gene_column == "") {
    rep(NA_character_, nrow(expr))
  } else {
    as.character(Biobase::fData(eset)[[gene_column]])
  }
  list(
    matrix = data.frame(ID = rownames(expr), Symbol = genes, expr, check.names = FALSE),
    design = data.frame(sample = colnames(expr), study = study, condition = condition[keep])
  )
}
