settings_demo <- list(de_method = "DESeq2", meta_method = "fisher", batch_method = "none",
                      batch_plots = FALSE, reference = "Cont")

# Runs analysis.R of an unpacked bundle with Rscript, as a user would (the
# setup step is skipped). When the package is loaded from source (devtools),
# the script loads the same source instead of an installed version.
run_bundle_script <- function(dir) {
  script <- "analysis.R"
  if (requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("TranscriptMetaAnalyst")) {
    lines <- readLines(file.path(dir, script))
    lines[lines == "library(TranscriptMetaAnalyst)"] <-
      sprintf("pkgload::load_all(%s, quiet = TRUE)", deparse(pkgload::pkg_path(testthat::test_path())))
    script <- "analysis_dev.R"
    writeLines(lines, file.path(dir, script))
  }
  out <- withr::with_dir(dir, system2(file.path(R.home("bin"), "Rscript"), script,
                                      stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) stop(paste(utils::tail(out, 20), collapse = "\n"))
  invisible(out)
}

unzip_bundle <- function(zip_file) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  zip::unzip(zip_file, exdir = dir)
  file.path(dir, "tma_analysis")
}

test_that("tma_plot_options has the app's defaults and checks the format", {
  o <- tma_plot_options()
  expect_equal(o$heatmap_k, 2)
  expect_equal(o$volcano_p, 0.05)
  expect_identical(o$format, "pdf")
  expect_null(names(o$format))
  expect_error(tma_plot_options(format = "gif"))
})

test_that("save_tma_results writes the app's files with the chosen options", {
  res <- demo_result()
  dir <- withr::local_tempdir()
  files <- save_tma_results(res, dir, "T_", tma_plot_options(heatmap_top = 20, format = "png"))
  expect_true(all(file.exists(files)))
  expect_setequal(basename(files), paste0("T_", c(
    "meta_analysis_results.csv", "A_DEGs.csv", "B_DEGs.csv", "C_DEGs.csv", "pdist.png", "pca.png",
    "heatmap.png", "heatmap_data.csv", "venndiagram.png", "upset.png", "volcano_plot.png",
    "volcano_data.csv")))
  expect_equal(nrow(utils::read.csv(file.path(dir, "T_heatmap_data.csv"))), 20)
  meta <- utils::read.csv(file.path(dir, "T_meta_analysis_results.csv"), row.names = 1)
  expect_equal(meta$adjpval, res$meta$adjpval)
})

test_that("save_tma_results writes enrichment tables and plots", {
  res <- demo_result()
  dir <- withr::local_tempdir()
  enr <- structure(list(service = "STRING",
                        table = data.frame(category = "Process", term = "GO:1", description = "photosynthesis",
                                           number_of_genes = 5, fdr = 0.001),
                        mapping = mapping_report(res$meta_degs, res$meta_degs, "STRING"),
                        version = "12.0", background = "all tested genes", organism = "Synechocystis",
                        id_type = "kegg"), class = "tma_enrichment")
  local_mocked_bindings(plot_string_lollipop = function(table, ...) ggplot2::ggplot())
  files <- save_tma_results(res, dir, "", tma_plot_options(format = "png"), list(string = enr))
  expect_true(all(c("string_enrichment.txt", "string_plot.png") %in% basename(files)))
  expect_match(readLines(file.path(dir, "string_enrichment.txt"))[1], "STRING")
})

test_that("the script carries the plot options and re-runs the enrichment done in the app", {
  script <- tma_script(settings_demo, data = "bundle",
                       plot_options = list(`Treated vs Cont` = tma_plot_options(heatmap_top = 20, volcano_p = 0.01)),
                       enrichment = list(`Treated vs Cont` = list(id_type = "kegg", string = 1148, run = "string")))
  text <- paste(script, collapse = "\n")
  expect_match(text, 'read_tma_table("data/counts.txt")', fixed = TRUE)
  expect_match(text, "heatmap_top = 20")
  expect_match(text, "volcano_p = 0.01")
  expect_match(text, 'run = "string"')
  expect_match(text, "run_tma_enrichment(res", fixed = TRUE)
  expect_no_error(parse(text = script))
  # without enrichment in the app, the enrichment section stays commented out
  plain <- tma_script(settings_demo, data = "demo",
                      enrichment = list(`Treated vs Cont` = list(id_type = "kegg", kegg = "syn")))
  expect_match(plain, '#   kegg <- run_kegg_enrichment(res$meta_degs, organism = "syn"', all = FALSE, fixed = TRUE)
})

test_that("an RNA-Seq bundle re-runs with Rscript and reproduces the app's tables", {
  skip_on_cran()
  demo <- tma_demo_data()
  analysis <- run_tma_analysis(demo$counts, demo$design, batch_method = "none")
  opts <- list(`Treated vs Cont` = tma_plot_options(heatmap_top = 25, volcano_p = 0.01, volcano_labels = "none",
                                                    format = "png", width = 12))
  zip_file <- file.path(withr::local_tempdir(), "tma_analysis.zip")
  write_tma_bundle(zip_file, analysis, demo, settings_demo, source = "demo", plot_options = opts, run_time = 5)
  dir <- unzip_bundle(zip_file)
  expect_true(all(file.exists(file.path(dir, c("README.md", "setup.R", "analysis.R", "settings.json",
                                               "session_info.txt", "data/counts.txt", "data/design.txt")))))
  expect_equal(read_tma_table(file.path(dir, "data", "counts.txt")), demo$counts)

  # the plot options recorded in settings.json reach analysis.R
  json <- jsonlite::read_json(file.path(dir, "settings.json"))
  recorded <- json$comparisons$`Treated vs Cont`$plot_options
  expect_equal(recorded$heatmap_top, 25)
  expect_equal(recorded$volcano_labels, "none")
  script <- paste(readLines(file.path(dir, "analysis.R")), collapse = "\n")
  for (option in c("heatmap_top = 25", "volcano_p = 0.01", 'volcano_labels = "none"', 'format = "png"', "width = 12")) {
    expect_match(script, option, fixed = TRUE)
  }
  expect_no_error(parse(text = readLines(file.path(dir, "setup.R"))))

  run_bundle_script(dir)
  original <- list.files(file.path(dir, "results"))
  expect_setequal(list.files(file.path(dir, "rerun")), original)
  for (table in grep("\\.csv$", original, value = TRUE)) {
    expect_identical(readLines(file.path(dir, "rerun", table)), readLines(file.path(dir, "results", table)))
  }
  expect_equal(nrow(utils::read.csv(file.path(dir, "rerun", "Treated_vs_Cont_heatmap_data.csv"))), 25)
})

test_that("a microarray bundle re-runs with Rscript and reproduces the app's tables", {
  skip_on_cran()
  demo <- tma_microarray_demo_data()
  analysis <- run_array_analysis(demo$matrices, demo$design, reference = "HC")
  settings <- list(meta_method = "fisher", batch_method = "none", batch_plots = FALSE, reference = "HC")
  zip_file <- file.path(withr::local_tempdir(), "bundle.zip")
  write_tma_bundle(zip_file, analysis, demo, settings, platform = "array", source = "demo")
  dir <- unzip_bundle(zip_file)
  expect_setequal(list.files(file.path(dir, "data")), c(paste0(names(demo$matrices), ".txt"), "design.txt"))
  run_bundle_script(dir)
  tables <- grep("meta_analysis_results\\.csv$", list.files(file.path(dir, "results")), value = TRUE)
  expect_length(tables, 2)
  for (table in tables) {
    expect_identical(readLines(file.path(dir, "rerun", table)), readLines(file.path(dir, "results", table)))
  }
})

test_that("a GEO bundle keeps the downloaded series and the sample-to-condition mapping", {
  res <- demo_result()
  analysis <- fake_analysis(res)
  cache <- withr::local_tempdir()
  writeLines("series", file.path(cache, "GSE1_series_matrix.txt.gz"))
  writeLines("platform", file.path(cache, "GPL9.soft.gz"))
  writeLines("other", file.path(cache, "GSE2_series_matrix.txt.gz"))
  geo <- list(list(gse = "GSE1", platform = "GPL9", condition_column = "title",
                   labels = c(s1 = "Cont", s2 = "Treated"), gene_column = "ORF"))
  matrices <- list(GSE1 = data.frame(ID = c("p1", "p2"), Symbol = c("a", "b"), s1 = 1:2, s2 = 3:4))
  dir <- file.path(withr::local_tempdir(), "bundle")
  write_tma_bundle(dir, analysis, list(matrices = matrices, design = res$design),
                   list(meta_method = "fisher"), platform = "array", source = "geo", geo = geo, geo_cache = cache)
  expect_setequal(list.files(file.path(dir, "data", "geo")),
                  c("GSE1_series_matrix.txt.gz", "GPL9.soft.gz", "sample_conditions.csv"))
  mapping <- utils::read.csv(file.path(dir, "data", "geo", "sample_conditions.csv"))
  expect_equal(mapping$condition, c("Cont", "Treated"))
  expect_true(file.exists(file.path(dir, "data", "GSE1.txt")))
  expect_match(readLines(file.path(dir, "analysis.R")), 'read_array_matrix("data/GSE1.txt")', all = FALSE, fixed = TRUE)
})

test_that("the RNA-Seq tab offers the run as a reproducible analysis ZIP", {
  res <- demo_result()
  local_mocked_bindings(run_tma_analysis = function(counts, design, ...) fake_analysis(res))
  mock_organism_lists()
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(rnaseq_demo = 1)
    for (i in 1:200) {
      later::run_now(0.05)
      session$flushReact()
      if (rnaseq_task$status() != "running") break
    }
    session$setInputs(`rnaseq_c1-heatmap_top` = 30, `rnaseq_c1-volcano_p` = 0.001)
    zip_file <- output$rnaseq_bundle
    expect_match(basename(zip_file), "\\.zip$")
    listed <- zip::zip_list(zip_file)$filename
    expect_true(all(paste0("tma_analysis/", c("README.md", "setup.R", "analysis.R", "settings.json",
                                              "session_info.txt", "data/counts.txt", "data/design.txt")) %in% listed))
    expect_true(any(grepl("^tma_analysis/results/.*meta_analysis_results\\.csv$", listed)))
    dir <- unzip_bundle(zip_file)
    script <- paste(readLines(file.path(dir, "analysis.R")), collapse = "\n")
    expect_match(script, "heatmap_top = 30")
    expect_match(script, "volcano_p = 0.001")
    # the demo's preselected enrichment choices are recorded
    json <- jsonlite::read_json(file.path(dir, "settings.json"))
    expect_equal(json$source, "demo")
  })
})
