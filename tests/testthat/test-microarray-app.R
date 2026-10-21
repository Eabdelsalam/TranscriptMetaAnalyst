wait_for_task <- function(task, session) {
  for (i in 1:400) {
    later::run_now(0.05)
    session$flushReact()
    if (task$status() != "running") break
  }
  # one more tick so the outputs see the finished task
  later::run_now(0.1)
  session$flushReact()
}

test_that("the demo button runs the microarray demo and fills the results", {
  mock_organism_lists()
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(array_demo = 1)
    wait_for_task(array_task, session)
    expect_equal(array_task$status(), "success")
    expect_equal(array_analysis()$contrasts, c("LC_24h vs HC", "LC_3h vs HC"))
    prep <- as.character(output$array_preparation$html)
    expect_match(prep, "3358 genes are measured in all 3 studies")
    expect_match(output$array_preparation_table, "GSE82120")
    results <- as.character(output$array_results$html)
    expect_match(results, "LC_24h vs HC")
    expect_match(results, "Study GSE82120 has no samples of")
    expect_match(as.character(output$`array_c1-num_studies`$html), "limma")
    script <- readLines(output$array_script)
    expect_match(script, "tma_microarray_demo_data()", all = FALSE, fixed = TRUE)
    expect_match(script, "run_array_analysis(", all = FALSE, fixed = TRUE)
    expect_no_error(parse(text = script))
  })
})

test_that("uploaded microarray files are checked and analysed", {
  mock_organism_lists()
  demo <- tma_microarray_demo_data()
  dir <- withr::local_tempdir()
  files <- vapply(names(demo$matrices), function(s) {
    f <- file.path(dir, paste0(s, ".txt"))
    utils::write.table(demo$matrices[[s]], f, sep = "\t", quote = FALSE, row.names = FALSE)
    f
  }, character(1))
  design_file <- file.path(dir, "design.txt")
  utils::write.table(demo$design, design_file, sep = "\t", quote = FALSE, row.names = FALSE)
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(array_source = "upload",
                      array_files = data.frame(name = basename(files), datapath = unname(files)),
                      array_design = list(datapath = design_file, name = "design.txt"),
                      `array_settings-batch_method` = "model", array_analyze = 1)
    # a batch option without a batch column is refused
    expect_equal(array_task$status(), "initial")

    session$setInputs(`array_settings-batch_method` = "none", `array_settings-reference` = "HC",
                      `array_settings-meta_method` = "invnorm", array_analyze = 2)
    wait_for_task(array_task, session)
    expect_equal(array_task$status(), "success")
    expect_equal(array_analysis()$results[[1]]$meta_method, "invnorm")
    script <- readLines(output$array_script)
    expect_match(script, "PLACEHOLDERS", all = FALSE)
    expect_match(script, 'read_array_matrix("PATH/TO/GSE82120.txt")', all = FALSE, fixed = TRUE)
  })
})

test_that("GEO series are imported with the chosen columns and labels", {
  mock_organism_lists()
  local_mocked_bindings(fetch_geo = function(gse, ...) {
    list(GPL123 = fake_eset(prefix = paste0(gse, "_"), n_genes = 40))
  })
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(array_source = "geo", geo_ids = "GSE1, gse2", geo_fetch = 1)
    expect_named(geo_sets(), c("GSE1", "GSE2"))
    expect_match(as.character(output$geo_config$html), "GSE2")
    for (i in 1:2) {
      do.call(session$setInputs, stats::setNames(list("characteristics_ch1", "ORF"),
                                                 paste0(c("geo_cond_", "geo_gene_"), i)))
      do.call(session$setInputs, stats::setNames(list("Control", "Stress", ""),
                                                 paste0("geo_label_", i, "_", 1:3)))
    }
    expect_match(as.character(output$geo_labels_1$html), "treatment: stress \\(2 samples\\)")
    design <- array_design_input()
    expect_equal(nrow(design), 8)
    expect_setequal(design$condition, c("Control", "Stress"))

    session$setInputs(`array_settings-reference` = "Control", array_analyze = 1)
    wait_for_task(array_task, session)
    expect_equal(array_task$status(), "success")
    expect_equal(array_analysis()$contrasts, "Stress vs Control")
    script <- paste(readLines(output$array_script), collapse = "\n")
    expect_match(script, 'fetch_geo("GSE1")[["GPL123"]]', fixed = TRUE)
    expect_match(script, 'condition_column = "characteristics_ch1"', fixed = TRUE)
    expect_match(script, '"treatment: stress" = "Stress"', fixed = TRUE)
  })
})
