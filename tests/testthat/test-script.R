default_settings <- list(de_method = "DESeq2", meta_method = "fisher", batch_method = "none",
                         batch_plots = FALSE, reference = NULL)

# Runs a generated script in a temporary folder; the package is already
# loaded, so the library() line is skipped
run_script <- function(lines) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  withr::local_dir(dir, .local_envir = parent.frame())
  env <- new.env(parent = globalenv())
  eval(parse(text = lines[lines != "library(TranscriptMetaAnalyst)"]), envir = env)
  env
}

test_that("the generated scripts parse", {
  for (data in c("upload", "demo")) {
    expect_no_error(parse(text = tma_script(default_settings, data = data)))
  }
  mixed <- modifyList(default_settings, list(de_method = c(A = "edgeR", B = "DESeq2", C = "limma-voom"),
                                             reference = "Cont"))
  expect_no_error(parse(text = tma_script(mixed)))
})

test_that("the upload script flags its placeholder paths", {
  script <- tma_script(default_settings, counts_name = "my_counts.txt", design_name = "my_design.txt")
  expect_match(script, "file paths below are PLACEHOLDERS", all = FALSE)
  expect_match(script, 'read_tma_table\\("PATH/TO/my_counts.txt"\\)', all = FALSE)
  expect_match(script, 'read_tma_table\\("PATH/TO/my_design.txt"\\)', all = FALSE)
  expect_false(any(grepl("PLACEHOLDERS", tma_script(default_settings, data = "demo"))))
})

test_that("the script records the settings used in the app", {
  settings <- list(de_method = c(A = "edgeR", B = "DESeq2"), meta_method = "invnorm",
                   batch_method = "combat_seq", batch_plots = TRUE, reference = "Cont")
  script <- paste(tma_script(settings, enrichment = list(id_type = "symbol", string = 9606)),
                  collapse = "\n")
  expect_match(script, 'reference = "Cont"')
  expect_match(script, 'de_method = c\\(A = "edgeR", B = "DESeq2"\\)')
  expect_match(script, 'meta_method = "invnorm"')
  expect_match(script, 'batch_method = "combat_seq"')
  expect_match(script, "batch_plots = TRUE")
  expect_match(script, '# id_type <- "symbol"')
  expect_match(script, "species = 9606")
  expect_match(script, "organism = <PANTHER taxon ID>")
})

test_that("the demo script runs and reproduces the app's result", {
  settings <- modifyList(default_settings, list(meta_method = "invnorm"))
  env <- run_script(tma_script(settings, data = "demo"))
  res <- env$analysis$results[["Treated vs Cont"]]
  demo <- tma_demo_data()
  expected <- run_tma_pipeline(demo$counts, demo$design, meta_method = "invnorm", batch_method = "none")
  expect_equal(res$meta_degs, expected$meta_degs)
  files <- list.files("tma_results")
  expect_setequal(files, paste0("Treated_vs_Cont_",
                                c("meta_analysis_results.csv", "pca.pdf", "volcano.pdf",
                                  "heatmap.pdf", "venn.pdf", "upset.pdf")))
})

test_that("the upload script runs once the placeholder paths are filled in", {
  script <- tma_script(default_settings, counts_name = "counts.txt", design_name = "design.txt")
  script <- sub("PATH/TO/counts.txt", tma_example_file("counts"), script, fixed = TRUE)
  script <- sub("PATH/TO/design.txt", tma_example_file("design"), script, fixed = TRUE)
  env <- run_script(script)
  expect_equal(env$analysis$results[[1]]$meta_degs, demo_result()$meta_degs)
})

test_that("the app offers the scripts with the settings of the run", {
  res <- demo_result()
  local_mocked_bindings(run_tma_analysis = function(counts, design, reference, de_method,
                                                     meta_method, batch_method, batch_plots, ...) {
    fake_analysis(res)
  })
  mock_organism_lists()
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(rnaseq_counts = list(datapath = tma_example_file("counts"), name = "my_counts.txt"),
                      rnaseq_design = list(datapath = tma_example_file("design"), name = "my_design.txt"),
                      `rnaseq_settings-meta_method` = "invnorm",
                      rnaseq_analyze = 1)
    for (i in 1:200) {
      later::run_now(0.05)
      session$flushReact()
      if (rnaseq_task$status() != "running") break
    }
    session$setInputs(`rnaseq_c1-id_type` = "kegg", `rnaseq_c1-organism_kegg` = "syn")
    script <- readLines(output$rnaseq_script)
    expect_match(script, "PATH/TO/my_counts.txt", all = FALSE)
    expect_match(script, 'meta_method = "invnorm"', all = FALSE)
    expect_match(script, 'reference = "Cont"', all = FALSE)
    expect_match(script, 'organism = "syn"', all = FALSE)
    expect_no_error(parse(text = script))

    # the demo on the same tab gives a script with the bundled demo data
    session$setInputs(`rnaseq_settings-de_method` = "edgeR", rnaseq_demo = 1)
    for (i in 1:200) {
      later::run_now(0.05)
      session$flushReact()
      if (rnaseq_task$status() != "running") break
    }
    demo_script <- readLines(output$rnaseq_script)
    expect_match(demo_script, "demo <- tma_demo_data()", all = FALSE, fixed = TRUE)
    expect_match(demo_script, 'de_method = "edgeR"', all = FALSE)
    expect_false(any(grepl("PATH/TO", demo_script)))
  })
})
