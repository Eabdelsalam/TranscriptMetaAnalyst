test_that("mod_settings_server returns the global methods by default", {
  shiny::testServer(mod_settings_server, {
    session$setInputs(de_method = "edgeR", meta_method = "invnorm")
    expect_equal(session$returned(), list(de_method = "edgeR", meta_method = "invnorm",
                                          batch_method = "none", batch_plots = FALSE,
                                          reference = NULL))
  })
})

test_that("mod_settings_server applies per-study overrides", {
  design <- tma_demo_data()$design
  shiny::testServer(mod_settings_server, args = list(design = shiny::reactive(design)), {
    session$setInputs(de_method = "DESeq2", meta_method = "fisher")
    # no override chosen: one method for all studies
    expect_equal(session$returned()$de_method, "DESeq2")

    session$setInputs(study_method_B = "limma-voom", study_method_C = "")
    expect_equal(session$returned()$de_method,
                 c(A = "DESeq2", B = "limma-voom", C = "DESeq2"))
    expect_match(as.character(output$per_study$html), "Same as above")
  })
})

test_that("the demo is recomputed for new settings and cached per combination", {
  res <- demo_result()
  calls <- character(0)
  local_mocked_bindings(run_tma_pipeline = function(counts, design, de_method, meta_method, ...) {
    calls <<- c(calls, paste(de_method, meta_method))
    out <- res
    out$meta_method <- meta_method
    out
  })
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(`demo_settings-de_method` = "DESeq2", `demo_settings-meta_method` = "fisher",
                      main_tabs = "demo")
    expect_equal(demo_result()$meta_method, "fisher")
    session$setInputs(`demo_settings-meta_method` = "invnorm")
    expect_equal(demo_result()$meta_method, "invnorm")
    session$setInputs(`demo_settings-meta_method` = "fisher")
    expect_equal(demo_result()$meta_method, "fisher")
  })
  expect_equal(calls, c("DESeq2 fisher", "DESeq2 invnorm"))
})

test_that("the RNA-seq run uses the chosen settings", {
  res <- demo_result()
  local_mocked_bindings(run_tma_analysis = function(counts, design, reference, de_method,
                                                     meta_method, batch_method, batch_plots, ...) {
    out <- res
    out$de_methods <- de_method
    out$meta_method <- meta_method
    out$batch <- list(method = batch_method, plots = batch_plots)
    out$reference <- reference
    fake_analysis(out)
  })
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(rnaseq_counts = list(datapath = tma_example_file("counts")),
                      rnaseq_design = list(datapath = tma_example_file("design")),
                      `rnaseq_settings-de_method` = "edgeR",
                      `rnaseq_settings-meta_method` = "invnorm",
                      `rnaseq_settings-study_method_C` = "limma-voom",
                      `rnaseq_settings-reference` = "Treated",
                      rnaseq_analyze = 1)
    for (i in 1:200) {
      later::run_now(0.05)
      session$flushReact()
      if (rnaseq_task$status() != "running") break
    }
    res_run <- rnaseq_analysis()$results[[1]]
    expect_equal(res_run$de_methods, c(A = "edgeR", B = "edgeR", C = "limma-voom"))
    expect_equal(res_run$meta_method, "invnorm")
    expect_equal(res_run$batch, list(method = "none", plots = FALSE))
    expect_equal(res_run$reference, "Treated")
  })
})

test_that("the batch options are always shown on the RNA-seq tab, defaulting to Ignore", {
  rnaseq_ui <- as.character(mod_settings_ui("x", per_study = TRUE))
  expect_match(rnaseq_ui, "Batch in the DE model")
  expect_match(rnaseq_ui, "ComBat-seq on the raw counts")
  expect_match(rnaseq_ui, 'value="none" checked')
  expect_match(rnaseq_ui, "Use batch-corrected values")
  # the demo data have no batch column, so the Demo tab has no batch options
  expect_false(grepl("ComBat-seq on the raw counts", as.character(mod_settings_ui("y"))))

  shiny::testServer(mod_settings_server, {
    session$setInputs(batch_method = "combat_seq", batch_plots = TRUE)
    expect_equal(session$returned()$batch_method, "combat_seq")
    expect_true(session$returned()$batch_plots)
  })
})

test_that("validate_batch_choice needs a batch column for any batch option", {
  design <- tma_demo_data()$design
  expect_identical(validate_batch_choice(design, "none", FALSE), character(0))
  expect_match(validate_batch_choice(design, "model"), "'Batch in the DE model.*needs a 'batch' column")
  expect_match(validate_batch_choice(design, "combat_seq"), "'ComBat-seq on the raw counts' needs")
  expect_match(validate_batch_choice(design, "none", TRUE), "'Use batch-corrected values.*needs")
  design$batch <- rep(c("b1", "b2"), 6)
  expect_identical(validate_batch_choice(design, "model", TRUE), character(0))
})

test_that("the RNA-seq run is blocked when a batch option is chosen without a batch column", {
  called <- FALSE
  local_mocked_bindings(run_tma_pipeline = function(...) {
    called <<- TRUE
    demo_result()
  })
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(rnaseq_counts = list(datapath = tma_example_file("counts")),
                      rnaseq_design = list(datapath = tma_example_file("design")),
                      `rnaseq_settings-batch_method` = "model",
                      rnaseq_analyze = 1)
    expect_error(rnaseq_inputs(), class = "shiny.silent.error")
    expect_equal(rnaseq_task$status(), "initial")
  })
  expect_false(called)
})

test_that("the reference drop-down lists the design's conditions and guesses the control", {
  design <- multi_data_set()$design
  shiny::testServer(mod_settings_server, args = list(design = shiny::reactive(design)), {
    html <- as.character(output$reference_ui$html)
    expect_match(html, '<option value="Cont" selected>Cont</option>')
    expect_match(html, '<option value="Heat">Heat</option>')
    session$setInputs(reference = "Heat")
    expect_equal(session$returned()$reference, "Heat")
  })
})
