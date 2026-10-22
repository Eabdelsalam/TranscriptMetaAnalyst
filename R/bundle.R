#' Plot and display options
#'
#' All display options of a results tab: what the heatmap and volcano plot
#' show and the format and size of the saved figures. They change only the
#' figures, never the DEGs. The app keeps one set per comparison, and
#' [save_tma_results()] and the reproducible analysis bundle use them.
#'
#' @param heatmap_k Number of gene clusters in the heatmap.
#' @param heatmap_top Show only the top N DEGs (smallest adjusted p-values);
#'   `NA` shows all.
#' @param heatmap_padj,heatmap_lfc Show only DEGs with an adjusted p-value of at
#'   most `heatmap_padj` (`NA`: the significance level of the analysis) and an
#'   absolute mean log2 fold change of at least `heatmap_lfc`.
#' @param heatmap_palette One of [heatmap_palettes()].
#' @param heatmap_method Clustering method (`"complete"`, `"average"`,
#'   `"ward.D2"` or `"single"`).
#' @param heatmap_cluster_rows,heatmap_cluster_cols Cluster genes and samples.
#' @param heatmap_names Show gene names.
#' @param volcano_p,volcano_fc Adjusted p-value and |log2 fold change| cutoffs
#'   drawn on the volcano plot.
#' @param volcano_labels Gene labels: `"significant"`, `"top"` or `"none"`.
#' @param volcano_top N for `volcano_labels = "top"`.
#' @param volcano_ymax Y-axis maximum (-log10 adjusted p-value); `NA` shows all.
#' @param format File format of the figures, one of [plot_formats()].
#' @param width,height Figure size in cm; `NA` uses each figure's default size.
#' @return A named list.
#' @export
#' @examples
#' tma_plot_options(heatmap_top = 50, volcano_p = 0.01, format = "png")
tma_plot_options <- function(heatmap_k = 2, heatmap_top = NA, heatmap_padj = NA, heatmap_lfc = 0,
                             heatmap_palette = "yellow_black_purple", heatmap_method = "complete",
                             heatmap_cluster_rows = TRUE, heatmap_cluster_cols = TRUE,
                             heatmap_names = FALSE, volcano_p = 0.05, volcano_fc = 1,
                             volcano_labels = "significant", volcano_top = 10, volcano_ymax = NA,
                             format = "pdf", width = NA, height = NA) {
  list(heatmap_k = heatmap_k, heatmap_top = heatmap_top, heatmap_padj = heatmap_padj,
       heatmap_lfc = heatmap_lfc, heatmap_palette = heatmap_palette, heatmap_method = heatmap_method,
       heatmap_cluster_rows = heatmap_cluster_rows, heatmap_cluster_cols = heatmap_cluster_cols,
       heatmap_names = heatmap_names, volcano_p = volcano_p, volcano_fc = volcano_fc,
       volcano_labels = volcano_labels, volcano_top = volcano_top, volcano_ymax = volcano_ymax,
       format = unname(match.arg(format, plot_formats())), width = width, height = height)
}

# Complete a (possibly partial) list of plot options with the defaults
complete_plot_options <- function(options) {
  utils::modifyList(tma_plot_options(), as.list(options))
}

# Default figure sizes (cm) of a result, as used by the app's downloads
plot_sizes <- function(result) {
  n <- length(result$studies)
  rows <- max(1, n %/% 2)
  list(pdist = c(15.24, 11.43), pca = c(20, 16), heatmap = c(17.78, 25.4), venn = c(15, 15),
       upset = c(max(6, n * 2) * 2.54, max(4, n) * 2.54), volcano = c(30, 25),
       batch = c(25, 15), enrichment = c(30, 20), counts = c(ceiling(n / rows) * 4 * 2.54, rows * 4 * 2.54))
}

# The heatmap as drawn with the given options: the matrix and how to draw it
heatmap_display <- function(result, options) {
  o <- complete_plot_options(options)
  genes <- heatmap_genes(result, padj_cutoff = number_or(o$heatmap_padj, result$alpha),
                         lfc_cutoff = number_or(o$heatmap_lfc, 0), top_n = number_or(o$heatmap_top, NA))
  mat <- heatmap_matrix(result, genes = genes)
  list(mat = mat, k_row = heatmap_row_clusters(nrow(mat), o$heatmap_k), k_col = length(result$studies),
       palette = o$heatmap_palette, cluster_rows = !isFALSE(o$heatmap_cluster_rows),
       cluster_cols = !isFALSE(o$heatmap_cluster_cols), method = o$heatmap_method,
       show_rownames = isTRUE(o$heatmap_names))
}

# The volcano plot as drawn with the given options
volcano_display <- function(result, options, volcano = volcano_data(result)) {
  o <- complete_plot_options(options)
  plot_volcano(volcano, p_cutoff = number_or(o$volcano_p, 0.05), fc_cutoff = number_or(o$volcano_fc, 1),
               labels = o$volcano_labels, top_n = number_or(o$volcano_top, 10),
               y_max = number_or(o$volcano_ymax, NA))
}

#' Save the tables and figures of a comparison
#'
#' Writes the files the app offers for download for one comparison: the
#' meta-analysis table, the DEGs of each study, the p-value histogram, PCA,
#' heatmap (and its data), Venn diagram (up to 4 sets), UpSet plot, volcano plot
#' (and its data), the batch diagnostics when the design had a batch column,
#' and the tables and figures of any enrichment results given.
#'
#' @param result A `tma_result` (one comparison of [run_tma_analysis()] or
#'   [run_array_analysis()]).
#' @param dir Output folder (created if needed).
#' @param prefix Start of the file names, e.g. `"Treated_vs_Cont_"`.
#' @param options Plot options, see [tma_plot_options()] (a partial list is
#'   completed with the defaults).
#' @param enrichment `NULL`, or a named list of `tma_enrichment` results (names
#'   `go`, `kegg` and/or `string`), e.g. from [run_tma_enrichment()].
#' @return The paths of the files written, invisibly.
#' @export
#' @examples
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' out <- tempfile()
#' files <- save_tma_results(res, out, "Treated_vs_Cont_",
#'                           tma_plot_options(heatmap_top = 50, format = "png"))
#' basename(files)
save_tma_results <- function(result, dir, prefix = "", options = tma_plot_options(), enrichment = NULL) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  o <- complete_plot_options(options)
  sizes <- plot_sizes(result)
  written <- character(0)
  path <- function(name) file.path(dir, paste0(prefix, name))
  csv <- function(x, name, ...) {
    file <- path(name)
    utils::write.csv(x, file, quote = FALSE, ...)
    written <<- c(written, file)
  }
  figure <- function(name, plot, size) {
    file <- path(paste0(name, ".", o$format))
    save_plot(file, plot, format = o$format, width = o$width, height = o$height,
              default_width = size[1], default_height = size[2])
    written <<- c(written, file)
  }

  csv(result$meta, "meta_analysis_results.csv")
  for (study in result$studies) {
    csv(study_significant(result, study), paste0(study, "_DEGs.csv"), row.names = FALSE)
  }
  figure("pdist", function() plot_pvalue_hist(result), sizes$pdist)
  figure("pca", plot_pca(result), sizes$pca)
  heat <- heatmap_display(result, o)
  if (nrow(heat$mat) >= 2) {
    figure("heatmap", function() {
      plot_heatmap_static(heat$mat, result$design, heat$k_row, heat$k_col, palette = heat$palette,
                          cluster_rows = heat$cluster_rows, cluster_cols = heat$cluster_cols,
                          method = heat$method, show_rownames = heat$show_rownames)
    }, sizes$heatmap)
    csv(heat$mat, "heatmap_data.csv")
  }
  sets <- deg_sets(result)
  if (overlap_plot_type(sets) == "venn") figure("venndiagram", plot_venn(sets), sizes$venn)
  figure("upset", plot_upset(sets), sizes$upset)
  volcano <- volcano_data(result)
  figure("volcano_plot", volcano_display(result, o, volcano), sizes$volcano)
  csv(volcano, "volcano_data.csv", row.names = FALSE)
  if (!is.null(result$batch)) {
    figure("batch_pca", plot_batch_pca(result), sizes$batch)
    figure("batch_rle", plot_rle(result), sizes$batch)
  }
  lollipops <- list(go = plot_go_lollipop, kegg = plot_kegg_lollipop, string = plot_string_lollipop)
  for (key in intersect(names(enrichment), names(lollipops))) {
    enr <- enrichment[[key]]
    if (is.null(enr)) next
    file <- path(paste0(key, "_enrichment.txt"))
    write_enrichment(enr, file)
    written <- c(written, file)
    if (nrow(enr$table) > 0) figure(paste0(key, "_plot"), lollipops[[key]](enr$table), sizes$enrichment)
  }
  invisible(written)
}

#' Run the enrichment services chosen for a comparison
#'
#' Runs, for the meta-analysis DEGs of a comparison and against all tested
#' genes, each enrichment service listed in `params$run`. A service that fails
#' (e.g. without an internet connection) is reported with a message and left
#' out, so a script continues.
#'
#' @param result A `tma_result`.
#' @param params A list with `id_type`, `go` (PANTHER taxon ID), `kegg` (KEGG
#'   organism code), `string` (STRING taxon ID) and `run` (the services to run,
#'   any of `"go"`, `"kegg"` and `"string"`).
#' @return A named list of `tma_enrichment` results (possibly empty).
#' @export
#' @examples
#' \dontrun{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' run_tma_enrichment(res, list(id_type = "kegg", string = 1148, run = "string"))
#' }
run_tma_enrichment <- function(result, params) {
  out <- list()
  if (length(result$meta_degs) == 0) return(out)
  calls <- list(
    go = function() run_go_enrichment(result$meta_degs, organism = params$go, background = result$genes,
                                      id_type = params$id_type),
    kegg = function() run_kegg_enrichment(result$meta_degs, organism = params$kegg, background = result$genes,
                                          id_type = params$id_type),
    string = function() run_string_enrichment(result$meta_degs, species = params$string,
                                              background = result$genes, id_type = params$id_type))
  names <- c(go = "PANTHER", kegg = "KEGG", string = "STRING")
  for (key in intersect(params$run, names(calls))) {
    out[[key]] <- tryCatch(calls[[key]](), error = function(e) {
      message(names[[key]], " enrichment was not run: ", enrichment_message(e, names[[key]]))
      NULL
    })
  }
  out[!vapply(out, is.null, logical(1))]
}

#' Write a ready-to-run analysis bundle
#'
#' Writes everything needed to repeat an analysis on another computer: the
#' input data as used, a script that repeats the analysis with every setting of
#' the run (including the plot display options and enrichment choices), the
#' results of the run for comparison, a setup script that installs the same TMA
#' release and checks the package versions, the settings in machine-readable
#' form and the R session information. This is the "Download reproducible
#' analysis (ZIP)" of the app.
#'
#' The bundle contains `README.md`, `setup.R`, `analysis.R`, `settings.json`,
#' `session_info.txt`, `data/` (the input) and `results/` (the output of the
#' run). `analysis.R` reads `data/` with relative paths and writes the same
#' files to `rerun/`, so the two folders can be compared.
#'
#' @param file Path of the bundle: a `.zip` file, or a folder.
#' @param analysis The analysis, from [run_tma_analysis()] or
#'   [run_array_analysis()].
#' @param inputs The input data as used. RNA-seq: a list with `counts` and
#'   `design` (data frames). Microarrays: a list with `matrices` (a named list
#'   of data frames as read by [read_array_matrix()]; names become the file
#'   names) and `design`.
#' @param settings The analysis settings of the run, as for [tma_script()].
#' @param platform `"rnaseq"` or `"array"`.
#' @param source Where the data came from: `"upload"`, `"demo"` or `"geo"`
#'   (microarrays imported from GEO).
#' @param plot_options Plot options of every comparison: a list named by
#'   comparison, each as [tma_plot_options()] (missing comparisons get the
#'   defaults).
#' @param enrichment Enrichment choices of every comparison: a list named by
#'   comparison, each a list with `params` (as for [run_tma_enrichment()]) and,
#'   optionally, `results` (the `tma_enrichment` results of the run, saved in
#'   `results/`).
#' @param geo For `source = "geo"`: the import choices, as for [tma_script()].
#'   The downloaded series files found in the GEO cache are copied to
#'   `data/geo/`.
#' @param input_names Optional names of the uploaded files, for the README.
#' @param run_time Run time of the analysis in seconds, for the README.
#' @param session The R session information of the run.
#' @param geo_cache Folder of the downloaded GEO series.
#' @return `file`, invisibly.
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' analysis <- run_tma_analysis(demo$counts, demo$design, batch_method = "none")
#' settings <- list(de_method = "DESeq2", meta_method = "fisher", batch_method = "none",
#'                  batch_plots = FALSE, reference = "Cont")
#' zip_file <- file.path(tempdir(), "tma_analysis.zip")
#' write_tma_bundle(zip_file, analysis, demo, settings, source = "demo")
#' head(zip::zip_list(zip_file)$filename)
#' }
write_tma_bundle <- function(file, analysis, inputs, settings, platform = c("rnaseq", "array"),
                             source = c("upload", "demo", "geo"), plot_options = NULL, enrichment = NULL,
                             geo = NULL, input_names = NULL, run_time = NA,
                             session = utils::sessionInfo(), geo_cache = tma_geo_cache()) {
  platform <- match.arg(platform)
  source <- match.arg(source)
  zipped <- grepl("\\.zip$", file, ignore.case = TRUE)
  dir <- if (zipped) file.path(tempfile("tma-bundle"), "tma_analysis") else file
  dir.create(file.path(dir, "data"), showWarnings = FALSE, recursive = TRUE)
  contrasts <- names(analysis$results)
  # the reference the analysis used, also when it was guessed
  settings$reference <- analysis$reference
  options <- stats::setNames(lapply(contrasts, function(con) complete_plot_options(plot_options[[con]])), contrasts)
  enrich_params <- stats::setNames(lapply(contrasts, function(con) {
    p <- enrichment[[con]]$params
    if (is.null(p)) list(run = character(0)) else p
  }), contrasts)

  # data, exactly as used
  write_tsv <- function(x, name) {
    utils::write.table(x, file.path(dir, "data", name), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (platform == "rnaseq") {
    write_tsv(inputs$counts, "counts.txt")
    array_files <- NULL
  } else {
    # one file per study, named after the study or the uploaded file (without
    # its extension)
    stems <- make.names(tools::file_path_sans_ext(names(inputs$matrices), compression = TRUE), unique = TRUE)
    array_files <- stats::setNames(paste0(stems, ".txt"), names(inputs$matrices))
    for (study in names(inputs$matrices)) write_tsv(inputs$matrices[[study]], array_files[[study]])
  }
  write_tsv(inputs$design, "design.txt")
  if (source == "geo") write_geo_provenance(file.path(dir, "data", "geo"), geo, geo_cache)

  # results of the run
  for (con in contrasts) {
    save_tma_results(analysis$results[[con]], file.path(dir, "results"), contrast_prefix(con),
                     options[[con]], enrichment[[con]]$results)
  }

  # scripts, settings and session
  writeLines(tma_script(settings, data = "bundle", platform = platform, array_files = array_files,
                        enrichment = enrich_params, plot_options = options, output_dir = "rerun"),
             file.path(dir, "analysis.R"))
  writeLines(bundle_setup(session), file.path(dir, "setup.R"))
  writeLines(utils::capture.output(print(session)), file.path(dir, "session_info.txt"))
  json <- list(
    tma_version = tma_version(), tma_commit = tma_commit(), created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    platform = platform, source = source, input_files = as.list(input_names),
    data_files = list.files(file.path(dir, "data"), recursive = TRUE),
    analysis = c(settings[intersect(names(settings), c("de_method", "meta_method", "batch_method", "batch_plots"))],
                 list(reference = analysis$reference, alpha = analysis$results[[1]]$alpha)),
    geo = geo,
    comparisons = stats::setNames(lapply(contrasts, function(con) {
      list(plot_options = options[[con]], enrichment = enrich_params[[con]])
    }), contrasts))
  jsonlite::write_json(json, file.path(dir, "settings.json"), auto_unbox = TRUE, pretty = TRUE,
                       null = "null", na = "null", digits = NA)
  writeLines(bundle_readme(platform, source, contrasts, input_names, run_time),
             file.path(dir, "README.md"))

  if (zipped) {
    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
    # an absolute path: zip() writes from inside the temporary folder
    target <- file.path(normalizePath(dirname(file)), basename(file))
    zip::zip(target, files = basename(dir), root = dirname(dir))
    unlink(dirname(dir), recursive = TRUE)
  }
  invisible(file)
}

# File-name prefix of a comparison, e.g. "Treated_vs_Cont_"
contrast_prefix <- function(contrast) paste0(gsub("[^A-Za-z0-9]+", "_", contrast), "_")

# Version of the installed TMA, and the commit it was installed from (GitHub
# installs record it), or NA
tma_version <- function() {
  tryCatch(as.character(utils::packageVersion("TranscriptMetaAnalyst")), error = function(e) NA_character_)
}
tma_commit <- function() {
  d <- tryCatch(utils::packageDescription("TranscriptMetaAnalyst"), error = function(e) list())
  sha <- c(d$RemoteSha, d$GithubSHA1)
  if (length(sha) == 0) NA_character_ else sha[[1]]
}

# The downloaded GEO series (from the cache) and the sample-to-condition
# mapping used, for the record
write_geo_provenance <- function(dir, geo, cache) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  rows <- lapply(geo, function(g) {
    data.frame(gse = g$gse, platform = g$platform, condition_column = g$condition_column,
               value = names(g$labels), condition = unname(g$labels),
               gene_column = if (is.null(g$gene_column)) "" else g$gene_column)
  })
  utils::write.csv(do.call(rbind, rows), file.path(dir, "sample_conditions.csv"), row.names = FALSE)
  for (g in geo) {
    downloaded <- list.files(cache, pattern = paste0("^(", g$gse, "|", g$platform, ")[._-]"), full.names = TRUE)
    file.copy(downloaded, dir)
  }
}

# setup.R: R version check, installation of the same TMA release, and a
# comparison of package versions with the recorded session
bundle_setup <- function(session) {
  pkgs <- c(session$otherPkgs, session$loadedOnly)
  versions <- vapply(pkgs, function(p) p$Version, character(1))
  versions <- versions[order(names(versions))]
  r <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")
  version <- tma_version()
  c("# Installs the TranscriptMetaAnalyst release used for this analysis and checks",
    "# that the installed package versions match those of the original run",
    "# (session_info.txt). Run once, before analysis.R:  Rscript setup.R",
    "",
    sprintf("recorded_r <- %s", r(paste(session$R.version$major, session$R.version$minor, sep = "."))),
    "if (getRversion() < \"4.5.0\") stop(\"TranscriptMetaAnalyst needs R 4.5.0 or later; this is R \", getRversion(), \".\")",
    "if (getRversion() != recorded_r) {",
    "  message(\"Note: the analysis was run with R \", recorded_r, \"; this is R \", getRversion(), \".\")",
    "}",
    "",
    "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")",
    sprintf("tma_version <- %s", r(version)),
    "installed <- tryCatch(as.character(packageVersion(\"TranscriptMetaAnalyst\")), error = function(e) NA)",
    "if (!identical(installed, tma_version)) {",
    "  BiocManager::install(paste0(\"Eabdelsalam/TranscriptMetaAnalyst@v\", tma_version), update = FALSE, ask = FALSE)",
    "}",
    "",
    "# package versions of the original run",
    sprintf("recorded <- %s", r(versions)),
    "now <- vapply(names(recorded), function(p) {",
    "  tryCatch(as.character(packageVersion(p)), error = function(e) \"not installed\")",
    "}, character(1))",
    "differ <- now != recorded",
    "if (any(differ)) {",
    "  message(sum(differ), \" package(s) differ from the original run (results may differ slightly):\")",
    "  print(data.frame(package = names(recorded)[differ], recorded = recorded[differ],",
    "                   installed = now[differ], row.names = NULL))",
    "} else {",
    "  message(\"All \", length(recorded), \" packages match the original run.\")",
    "}")
}

# README.md of a bundle
bundle_readme <- function(platform, source, contrasts, input_names, run_time) {
  version <- tma_version()
  commit <- tma_commit()
  data_desc <- switch(source,
    demo = "the demo data bundled with TranscriptMetaAnalyst",
    geo = "microarray series imported from GEO (the extracted expression tables; the downloaded series files and the sample-to-condition mapping are in `data/geo/`)",
    upload = "the uploaded files")
  minutes <- if (is.na(run_time)) "a few minutes" else
    sprintf("about %s (the analysis took %.0f s in the app; saving the figures adds a little)",
            if (run_time < 60) "1 minute" else sprintf("%.0f minutes", ceiling(run_time / 60)), run_time)
  c("# Reproducible TranscriptMetaAnalyst analysis",
    "",
    sprintf("This folder repeats %s meta-analysis made with TranscriptMetaAnalyst (TMA) %s%s, on %s, with every setting of the run.",
            if (platform == "rnaseq") "an RNA-Seq" else "a microarray", version,
            if (is.na(commit)) "" else sprintf(" (commit %s)", substr(commit, 1, 7)), data_desc),
    if (length(input_names) > 0) c("", paste0("Original file names: ", paste(sprintf("`%s`", unlist(input_names)), collapse = ", "), ".")),
    "",
    "## Run it",
    "",
    "From this folder:",
    "",
    "```",
    "Rscript setup.R      # once: installs this TMA release and checks the package versions",
    "Rscript analysis.R   # repeats the analysis and writes the results to rerun/",
    "```",
    "",
    sprintf("Expected run time: %s.", minutes),
    "",
    "## Contents",
    "",
    "| File | Content |",
    "|---|---|",
    "| `analysis.R` | The analysis with all settings of the run: analysis methods, plot display options and enrichment choices. Reads `data/`, writes `rerun/`. |",
    "| `setup.R` | Installs the TMA release used and lists any package whose version differs from the original run. |",
    "| `settings.json` | Every analysis and display setting, machine-readable. |",
    "| `data/` | The input data as used in the analysis. |",
    "| `results/` | The tables and figures of the original run. |",
    "| `session_info.txt` | R and package versions of the original run. |",
    "",
    sprintf("Comparisons: %s.", paste(contrasts, collapse = ", ")),
    "",
    "## Checking the re-run",
    "",
    "`rerun/` gets the same files as `results/`. The tables (`*_meta_analysis_results.csv`, `*_DEGs.csv`, `*_heatmap_data.csv`, `*_volcano_data.csv`) should be identical when the package versions match; figures can differ in small details between systems. Enrichment results come from online services and need an internet connection; they can change when a service updates its database.")
}
