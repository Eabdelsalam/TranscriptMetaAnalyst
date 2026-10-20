test_that("mod_settings_server returns the global methods by default", {
  shiny::testServer(mod_settings_server, {
    session$setInputs(de_method = "edgeR", meta_method = "invnorm")
    expect_equal(session$returned(), list(de_method = "edgeR", meta_method = "invnorm",
                                          batch_method = "model", batch_plots = FALSE))
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
  local_mocked_bindings(run_tma_pipeline = function(counts, design, de_method, meta_method,
                                                     batch_method, batch_plots, ...) {
    out <- res
    out$de_methods <- de_method
    out$meta_method <- meta_method
    out$batch <- list(method = batch_method, plots = batch_plots)
    out
  })
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(rnaseq_counts = list(datapath = tma_example_file("counts")),
                      rnaseq_design = list(datapath = tma_example_file("design")),
                      `rnaseq_settings-de_method` = "edgeR",
                      `rnaseq_settings-meta_method` = "invnorm",
                      `rnaseq_settings-study_method_C` = "limma-voom",
                      rnaseq_analyze = 1)
    for (i in 1:200) {
      later::run_now(0.05)
      session$flushReact()
      if (rnaseq_task$status() != "running") break
    }
    expect_equal(rnaseq_result()$de_methods, c(A = "edgeR", B = "edgeR", C = "limma-voom"))
    expect_equal(rnaseq_result()$meta_method, "invnorm")
    expect_equal(rnaseq_result()$batch, list(method = "model", plots = FALSE))
  })
})

test_that("mod_settings_server offers batch options only for a design with batch", {
  design <- tma_demo_data()$design
  shiny::testServer(mod_settings_server, args = list(design = shiny::reactive(design)), {
    expect_error(output$batch, class = "shiny.silent.error")
  })
  design$batch <- rep(c("b1", "b2"), 6)
  shiny::testServer(mod_settings_server, args = list(design = shiny::reactive(design)), {
    expect_match(as.character(output$batch$html), "ComBat-seq")
    session$setInputs(batch_method = "combat_seq", batch_plots = TRUE)
    expect_equal(session$returned()$batch_method, "combat_seq")
    expect_true(session$returned()$batch_plots)
  })
})
