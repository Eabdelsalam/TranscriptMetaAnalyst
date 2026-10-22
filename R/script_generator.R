#' Reproducible R script of an analysis
#'
#' Writes an R script that repeats an analysis done in the app with the same
#' package functions: it reads the input files, runs [run_tma_analysis()] (or
#' [run_array_analysis()]) with the settings used in the app, and saves, with
#' [save_tma_results()], the tables and figures of every comparison with its
#' plot display options. Enrichment services that were run in the app are run
#' again (they need an internet connection); otherwise a commented-out
#' enrichment section is added.
#'
#' For uploaded data the file paths in the script are **placeholders**: the
#' app does not know where the files are on the user's computer. The script
#' starts with a banner that says so, and the lines to edit are marked. For the
#' demo, the script uses the demo data bundled with the package; in a
#' reproducible analysis bundle (`data = "bundle"`, see [write_tma_bundle()]),
#' it reads the files in `data/` with relative paths.
#'
#' @param settings A list with the analysis settings, as returned by the
#'   settings module: `de_method` (one method or a named vector with one per
#'   study), `meta_method`, `batch_method`, `batch_plots` and `reference`
#'   (`NULL` for the default reference).
#' @param data `"upload"` (placeholder paths), `"demo"` (bundled demo data),
#'   `"bundle"` (the files in `data/` of an analysis bundle) or, for
#'   microarrays, `"geo"` (studies imported from GEO; the script downloads them
#'   again).
#' @param counts_name,design_name Names of the uploaded files, quoted in the
#'   script to help the user find them.
#' @param platform `"rnaseq"` or `"array"` (microarray tab).
#' @param array_files Names of the microarray expression files (uploaded, or in
#'   `data/` of a bundle; for a bundle, named by study).
#' @param geo For `data = "geo"`: a list with one element per study, each a
#'   list with `gse`, `platform`, `condition_column`, `labels` (named
#'   character vector) and `gene_column`, as chosen in the app.
#' @param enrichment Enrichment choices: either one list with `id_type`, `go`,
#'   `kegg` and `string` (organisms) that fills in the commented-out enrichment
#'   section (missing values are left as placeholders), or a list named by
#'   comparison, each with `id_type`, `go`, `kegg`, `string` and `run` (the
#'   services run in the app, which the script runs again, see
#'   [run_tma_enrichment()]).
#' @param plot_options Plot display options: `NULL` (defaults), one list as
#'   from [tma_plot_options()] for all comparisons, or a list of them named by
#'   comparison.
#' @param output_dir Folder, relative to the working directory, where the
#'   script writes its results.
#' @return `tma_script()` returns the script as a character vector of lines;
#'   `write_tma_script()` writes it to `file` and returns `file` invisibly.
#' @export
#' @examples
#' settings <- list(de_method = "DESeq2", meta_method = "fisher",
#'                  batch_method = "none", batch_plots = FALSE, reference = NULL)
#' cat(tma_script(settings, data = "demo"), sep = "\n")
#' cat(tma_script(settings, data = "demo", plot_options = tma_plot_options(heatmap_top = 50)),
#'     sep = "\n")
tma_script <- function(settings, data = c("upload", "demo", "geo", "bundle"), counts_name = NULL,
                       design_name = NULL, enrichment = NULL, output_dir = "tma_results",
                       platform = c("rnaseq", "array"), array_files = NULL, geo = NULL,
                       plot_options = NULL) {
  data <- match.arg(data)
  platform <- match.arg(platform)
  if (data == "geo" && platform != "array") stop("GEO import is only used for microarray data.", call. = FALSE)
  r <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")
  setting <- function(name, default) if (is.null(settings[[name]])) default else settings[[name]]
  version <- tryCatch(as.character(utils::packageVersion("TranscriptMetaAnalyst")),
                      error = function(e) "unknown")

  banner <- if (data == "upload") {
    c("# ======================================================================",
      "# !! The input file paths below are PLACEHOLDERS.                     !!",
      "# !! Replace them with the paths of your files before running this    !!",
      "# !! script.                                                          !!",
      "# ======================================================================")
  }
  input <- if (platform == "array") {
    array_input(data, array_files, design_name, geo, r)
  } else if (data == "upload") {
    c("# Input files (EDIT THESE TWO PATHS)",
      sprintf("counts <- read_tma_table(\"PATH/TO/%s\")  # <- your counts file%s",
              if (is.null(counts_name)) "counts_file.txt" else counts_name,
              if (is.null(counts_name)) "" else " (uploaded in the app as this name)"),
      sprintf("design <- read_tma_table(\"PATH/TO/%s\")  # <- your design file%s",
              if (is.null(design_name)) "design_file.txt" else design_name,
              if (is.null(design_name)) "" else " (uploaded in the app as this name)"))
  } else if (data == "bundle") {
    c("# Input: the data of the original run, in data/",
      "counts <- read_tma_table(\"data/counts.txt\")",
      "design <- read_tma_table(\"data/design.txt\")")
  } else {
    c("# Input: the demo data bundled with TranscriptMetaAnalyst",
      "demo <- tma_demo_data()",
      "counts <- demo$counts",
      "design <- demo$design")
  }

  lines <- c(
    "# Reproducible analysis script generated by TranscriptMetaAnalyst (TMA)",
    sprintf("# TMA version %s, %s", version, format(Sys.Date())),
    "#",
    "# Repeats the analysis done in the app with the same settings and saves the",
    "# tables and figures of every comparison with its plot display options.",
    "# Change the plot options below to customise the figures; they never change",
    "# the DEGs.",
    banner,
    "",
    "library(TranscriptMetaAnalyst)",
    "",
    input,
    "",
    "# Analysis settings used in the app",
    if (platform == "array") "analysis <- run_array_analysis(" else "analysis <- run_tma_analysis(",
    if (platform == "array") "  matrices, design," else "  counts, design,",
    sprintf("  reference = %s,", r(setting("reference", NULL))),
    if (platform == "rnaseq") sprintf("  de_method = %s,", r(setting("de_method", "DESeq2"))),
    sprintf("  meta_method = %s,", r(setting("meta_method", "fisher"))),
    sprintf("  batch_method = %s,", r(setting("batch_method", "none"))),
    sprintf("  batch_plots = %s", r(isTRUE(setting("batch_plots", FALSE)))),
    ")",
    "print(analysis)",
    "",
    plot_options_code(plot_options, r),
    "",
    enrichment_code(enrichment, r),
    "",
    "# Save the tables and figures of every comparison",
    sprintf("output_dir <- %s", r(output_dir)),
    "for (contrast in names(analysis$results)) {",
    "  res <- analysis$results[[contrast]]",
    "  options <- if (is.null(plot_options[[contrast]])) default_options else plot_options[[contrast]]",
    "  enriched <- if (length(enrichment[[contrast]]$run) > 0) run_tma_enrichment(res, enrichment[[contrast]])",
    "  save_tma_results(res, output_dir, paste0(gsub(\"[^A-Za-z0-9]+\", \"_\", contrast), \"_\"),",
    "                   options, enriched)",
    "}",
    "",
    "# Software versions, for the methods section",
    "sessionInfo()"
  )
  # one script line per element (some sections are built as multi-line strings)
  unlist(lapply(lines, function(line) if (identical(line, "")) "" else strsplit(line, "\n", fixed = TRUE)[[1]]))
}

#' @rdname tma_script
#' @param file Path of the script file to write.
#' @param ... Passed to `tma_script()`.
#' @export
write_tma_script <- function(file, ...) {
  writeLines(tma_script(...), file)
  invisible(file)
}

# Plot options section: a tma_plot_options() call per comparison, with every
# option written out
plot_options_code <- function(plot_options, r) {
  call <- function(o, indent) {
    o <- complete_plot_options(o)
    args <- vapply(names(o), function(n) sprintf("%s = %s", n, r(o[[n]])), character(1))
    pad <- strrep(" ", indent)
    paste0("tma_plot_options(\n", pad, paste(args, collapse = paste0(",\n", pad)), ")")
  }
  per_contrast <- !is.null(plot_options) && !is.null(names(plot_options)) &&
    !any(names(plot_options) %in% names(tma_plot_options()))
  c("# Plot display options of each comparison (heatmap and volcano plot content,",
    "# figure format and size; width/height NA = default size)",
    sprintf("default_options <- %s", call(if (per_contrast) NULL else plot_options, 2)),
    if (per_contrast) {
      c("plot_options <- list(",
        paste(vapply(names(plot_options), function(n) {
          sprintf("  %s = %s", r(n), call(plot_options[[n]], 4))
        }, character(1)), collapse = ",\n"),
        ")")
    } else {
      "plot_options <- list()"
    })
}

# Enrichment section: the services run in the app, per comparison; otherwise a
# commented-out template filled with the choices made
enrichment_code <- function(enrichment, r) {
  keys <- c("id_type", "go", "kegg", "string", "run")
  per_contrast <- !is.null(enrichment) && !is.null(names(enrichment)) && !any(names(enrichment) %in% keys)
  ran <- per_contrast && any(vapply(enrichment, function(e) length(e$run) > 0, logical(1)))
  if (ran) {
    entries <- vapply(names(enrichment), function(n) {
      e <- enrichment[[n]][intersect(keys, names(enrichment[[n]]))]
      sprintf("  %s = %s", r(n), r(e))
    }, character(1))
    return(c("# Functional enrichment run in the app (needs an internet connection;",
             "# a service that cannot be reached is skipped with a message)",
             "enrichment <- list(",
             paste(entries, collapse = ",\n"),
             ")"))
  }
  single <- if (per_contrast) enrichment[[1]] else enrichment
  enrich <- function(name, placeholder) {
    value <- single[[name]]
    if (is.null(value) || identical(value, "")) placeholder else r(value)
  }
  c("# Optional: functional enrichment of the meta-analysis DEGs against all",
    "# tested genes (needs an internet connection). Remove the '#' and replace",
    "# any <...> placeholder to run it.",
    "enrichment <- list()",
    sprintf("# id_type <- %s  # one of: %s", enrich("id_type", "\"kegg\""),
            paste(sprintf("\"%s\"", id_types()), collapse = ", ")),
    "# for (contrast in names(analysis$results)) {",
    "#   res <- analysis$results[[contrast]]",
    sprintf("#   go <- run_go_enrichment(res$meta_degs, organism = %s, background = res$genes, id_type = id_type)",
            enrich("go", "<PANTHER taxon ID>")),
    sprintf("#   kegg <- run_kegg_enrichment(res$meta_degs, organism = %s, background = res$genes, id_type = id_type)",
            enrich("kegg", "\"<KEGG organism code>\"")),
    sprintf("#   string <- run_string_enrichment(res$meta_degs, species = %s, background = res$genes, id_type = id_type)",
            enrich("string", "<STRING taxon ID>")),
    "#   print(string)",
    "# }")
}

# Input section of a microarray script
array_input <- function(data, array_files, design_name, geo, r) {
  if (data == "demo") {
    return(c("# Input: the microarray demo data bundled with TranscriptMetaAnalyst",
             "demo <- tma_microarray_demo_data()",
             "matrices <- demo$matrices",
             "design <- demo$design"))
  }
  if (data == "bundle") {
    files <- vapply(names(array_files), function(study) {
      sprintf("  %s = read_array_matrix(%s)", r(study), r(file.path("data", array_files[[study]])))
    }, character(1))
    return(c("# Input: the expression data (one file per study) and design of the original",
             "# run, in data/",
             "matrices <- list(",
             paste(files, collapse = ",\n"),
             ")",
             "design <- read_tma_table(\"data/design.txt\")"))
  }
  if (data == "geo") {
    studies <- vapply(geo, function(g) {
      sprintf("  %s = geo_study(fetch_geo(%s)[[%s]], %s,\n    condition_column = %s,\n    labels = %s,\n    gene_column = %s)",
              g$gse, r(g$gse), r(g$platform), r(g$gse), r(g$condition_column), r(g$labels),
              r(if (is.null(g$gene_column) || g$gene_column == "") NULL else g$gene_column))
    }, character(1))
    return(c("# Input: the GEO series imported in the app (downloaded and cached;",
             "# needs an internet connection the first time)",
             "studies <- list(",
             paste(studies, collapse = ",\n"),
             ")",
             "matrices <- lapply(studies, `[[`, \"matrix\")",
             "design <- do.call(rbind, unname(lapply(studies, `[[`, \"design\")))"))
  }
  files <- if (length(array_files) == 0) "expression_file.txt" else array_files
  c("# Input files (EDIT THESE PATHS): one expression file per study, and the design file",
    "matrices <- list(",
    paste0(sprintf("  read_array_matrix(\"PATH/TO/%s\")", files), c(rep(",", length(files) - 1), ""),
           "  # <- your expression file"),
    ")",
    sprintf("design <- read_tma_table(\"PATH/TO/%s\")  # <- your design file",
            if (is.null(design_name)) "design_file.txt" else design_name))
}
