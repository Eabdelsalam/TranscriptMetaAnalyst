#' Load the bundled microarray demo dataset
#'
#' Normalised (log2) expression of wild-type *Synechocystis* sp. PCC 6803 grown
#' at high CO2 (`HC`) and 3 h or 24 h after a shift to low CO2 (`LC_3h`,
#' `LC_24h`), from four GEO series on two platforms: GSE82120, GSE63352 and
#' GSE68250 (GPL15867) and GSE31672 (GPL6546). To keep the package small, each
#' file holds at most 3 probes per gene (locus tag) and a few control probes.
#'
#' @return A list with `matrices` (a named list with one data frame per study,
#'   as read by [read_array_matrix()]) and `design` (data frame with `sample`,
#'   `study` and `condition`).
#' @export
#' @examples
#' demo <- tma_microarray_demo_data()
#' names(demo$matrices)
#' table(demo$design$study, demo$design$condition)
tma_microarray_demo_data <- function() {
  dir <- system.file("extdata", "microarray_demo", package = "TranscriptMetaAnalyst", mustWork = TRUE)
  files <- list.files(dir, pattern = "^GSE.*\\.txt\\.gz$", full.names = TRUE)
  matrices <- lapply(files, read_array_matrix)
  names(matrices) <- sub("\\.txt\\.gz$", "", basename(files))
  list(matrices = matrices, design = read_tma_table(file.path(dir, "design.txt")))
}

#' Read a microarray expression file
#'
#' Reads one study's normalised expression matrix: a tab-delimited file whose
#' first column holds the probe (or gene) IDs and whose other columns are the
#' samples. An optional column named `Symbol` (also accepted: `Gene`,
#' `Gene symbol`, `Gene_symbol`; any case) gives the gene each probe measures;
#' probes are then collapsed to genes (see [collapse_probes()]). Files may be
#' gzip-compressed.
#'
#' @param path Path to the file.
#' @return A data frame with the columns `ID`, `Symbol` (`NA` if the file has
#'   no symbol column) and one numeric column per sample.
#' @export
#' @examples
#' f <- tempfile(fileext = ".txt")
#' writeLines(c("ID\tSymbol\tS1\tS2", "p1\tgeneA\t5.1\t6.2", "p2\tgeneB\t7.3\t7.0"), f)
#' read_array_matrix(f)
read_array_matrix <- function(path) {
  x <- utils::read.delim(path, check.names = FALSE, quote = "", comment.char = "",
                         stringsAsFactors = FALSE)
  symbol_col <- which(tolower(gsub("[ ._]", "", names(x))) %in% c("symbol", "gene", "genesymbol"))
  symbol <- if (length(symbol_col) > 0) as.character(x[[symbol_col[1]]]) else rep(NA_character_, nrow(x))
  samples <- x[, -c(1, symbol_col), drop = FALSE]
  data.frame(ID = as.character(x[[1]]), Symbol = symbol, samples, check.names = FALSE)
}

#' Detect and apply a log2 transformation
#'
#' Uses the rule of GEO2R: values are considered not log-transformed when the
#' 99th percentile is above 100, or when the range is wider than 50 and the
#' 25th percentile is positive. Such data are log2-transformed (values <= 0
#' become `NA`).
#'
#' @param values Numeric matrix of expression values.
#' @return `needs_log2()` returns `TRUE` or `FALSE`. `log2_if_needed()` returns
#'   a list with `values` (log2 scale) and `transformed` (whether log2 was
#'   applied).
#' @export
#' @examples
#' needs_log2(matrix(c(10, 2000, 50000, 30), 2))
#' needs_log2(matrix(c(7.2, 10.5, 13.1, 8.4), 2))
needs_log2 <- function(values) {
  q <- stats::quantile(values, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
  (q[5] > 100) || (q[6] - q[1] > 50 && q[2] > 0)
}

#' @rdname needs_log2
#' @export
log2_if_needed <- function(values) {
  if (!needs_log2(values)) return(list(values = values, transformed = FALSE))
  values[!is.na(values) & values <= 0] <- NA
  list(values = log2(values), transformed = TRUE)
}

#' Collapse probes to genes
#'
#' Drops probes without a gene and probes that map to several genes (e.g.
#' `"A /// B"`), then keeps, for every gene, the probe with the highest mean
#' expression (D7). Vectorised: the probes are ordered by mean expression and
#' the first probe of each gene is kept.
#'
#' @param values Numeric matrix (probes x samples).
#' @param genes Character vector with the gene of each probe.
#' @return A list with `values` (genes x samples, gene names as row names) and
#'   `report` (numbers of `probes`, probes dropped for `no_gene` and
#'   `multiple_genes`, and resulting `genes`).
#' @export
#' @examples
#' m <- matrix(c(1, 5, 3, 9, 2, 6, 4, 8), 4, dimnames = list(paste0("p", 1:4), c("s1", "s2")))
#' collapse_probes(m, c("geneA", "geneA", "", "geneB /// geneC"))
collapse_probes <- function(values, genes) {
  genes <- trimws(as.character(genes))
  no_gene <- is.na(genes) | genes == "" | genes == "---"
  multiple <- !no_gene & grepl("///|//|,|;", genes)
  keep <- !no_gene & !multiple
  values <- values[keep, , drop = FALSE]
  genes <- genes[keep]
  best_first <- order(rowMeans(values, na.rm = TRUE), decreasing = TRUE)
  best_first <- best_first[!duplicated(genes[best_first])]
  best_first <- best_first[order(genes[best_first])]
  collapsed <- values[best_first, , drop = FALSE]
  rownames(collapsed) <- genes[best_first]
  list(values = collapsed,
       report = list(probes = length(keep), no_gene = sum(no_gene), multiple_genes = sum(multiple),
                     genes = nrow(collapsed)))
}

#' Prepare microarray studies for the meta-analysis
#'
#' Turns the uploaded (or imported) expression matrices into one gene-level
#' log2 matrix per study: each matrix is matched to its study through the
#' design, log2-transformed if needed ([log2_if_needed()]), collapsed from
#' probes to genes when it has a symbol column ([collapse_probes()]), and the
#' genes are intersected across studies (D7).
#'
#' @param matrices A list of data frames from [read_array_matrix()], one per
#'   study.
#' @param design Experimental design with `sample`, `study` and `condition`
#'   (and optionally `batch`). Every design sample must be a column of exactly
#'   one matrix, and each matrix must hold the samples of a single study.
#' @return A list with `studies` (one element per study, in design order, with
#'   `values`, a genes x samples log2 matrix of the common genes, and `design`)
#'   and `report` (a data frame with one row per study: `study`, `samples`,
#'   `probes`, `no_gene`, `multiple_genes`, `genes`, `log2_applied`) with the
#'   number of `common_genes` as an attribute, and `notes` (messages for the
#'   user).
#' @export
#' @examples
#' demo <- tma_microarray_demo_data()
#' prepared <- prepare_array_studies(demo$matrices, demo$design)
#' prepared$report
#' attr(prepared$report, "common_genes")
prepare_array_studies <- function(matrices, design) {
  design <- as.data.frame(design)
  study_ids <- unique(design$study)
  studies <- list()
  rows <- list()
  notes <- character(0)
  for (i in seq_along(matrices)) {
    m <- matrices[[i]]
    label <- if (!is.null(names(matrices)) && nzchar(names(matrices)[i])) names(matrices)[i] else paste("file", i)
    in_design <- intersect(names(m), design$sample)
    if (length(in_design) == 0) {
      stop("The expression file '", label, "' has no column that matches a sample of the design file.", call. = FALSE)
    }
    study <- unique(design$study[design$sample %in% in_design])
    if (length(study) > 1) {
      stop("The expression file '", label, "' contains samples of several studies (",
           paste(study, collapse = ", "), "); use one file per study.", call. = FALSE)
    }
    study_design <- design[design$study == study, , drop = FALSE]
    missing <- setdiff(study_design$sample, names(m))
    if (length(missing) > 0) {
      stop("Samples of study '", study, "' are missing from the expression file '", label, "': ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    values <- as.matrix(m[, study_design$sample, drop = FALSE])
    storage.mode(values) <- "double"
    rownames(values) <- m$ID
    logged <- log2_if_needed(values)
    if (logged$transformed) {
      notes <- c(notes, sprintf("Study %s: the values looked like linear intensities, so they were log2-transformed.", study))
    }
    has_symbol <- any(!is.na(m$Symbol) & m$Symbol != "")
    collapsed <- if (has_symbol) {
      collapse_probes(logged$values, m$Symbol)
    } else {
      keep <- !duplicated(m$ID)
      list(values = logged$values[keep, , drop = FALSE],
           report = list(probes = nrow(m), no_gene = 0L, multiple_genes = 0L, genes = sum(keep)))
    }
    rownames(study_design) <- study_design$sample
    studies[[study]] <- list(values = collapsed$values, design = study_design)
    rows[[study]] <- data.frame(study = study, samples = ncol(values), probes = collapsed$report$probes,
                                no_gene = collapsed$report$no_gene,
                                multiple_genes = collapsed$report$multiple_genes,
                                genes = collapsed$report$genes, log2_applied = logged$transformed)
  }
  missing_studies <- setdiff(study_ids, names(studies))
  if (length(missing_studies) > 0) {
    stop("No expression file was given for study ", paste(missing_studies, collapse = ", "), ".", call. = FALSE)
  }
  studies <- studies[study_ids]
  common <- Reduce(intersect, lapply(studies, function(s) rownames(s$values)))
  if (length(common) == 0) {
    stop("The studies have no gene in common. Check that the files use the same gene IDs (or symbols).", call. = FALSE)
  }
  common <- sort(common)
  studies <- lapply(studies, function(s) {
    s$values <- s$values[common, , drop = FALSE]
    s
  })
  report <- do.call(rbind, unname(rows[study_ids]))
  attr(report, "common_genes") <- length(common)
  notes <- c(notes, sprintf("%d genes are measured in all %d studies and are used in the meta-analysis.",
                            length(common), length(studies)))
  list(studies = studies, report = report, notes = notes)
}
