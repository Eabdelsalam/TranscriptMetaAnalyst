test_that("mod_results_server renders a tma_result", {
  res <- demo_result()
  shiny::testServer(mod_results_server,
                    args = list(result = shiny::reactive(res), organism = list(go = 1111708, kegg = "syn")), {
    expect_match(as.character(output$num_studies$html), "3 studies")
    expect_match(as.character(output$num_degs$html),
                 paste(length(res$meta_degs), "out of", length(res$genes)))

    session$setInputs(heatmap_k = 4)
    expect_equal(heat_k(), 4L)
    expect_equal(rownames(heat()), res$meta_degs)
    heat_csv <- utils::read.csv(output$download_heatmap_data, row.names = 1)
    expect_equal(dim(heat_csv), c(length(res$meta_degs), nrow(res$design)))

    expect_equal(volcano()$fc, unname(mean_fold_change(res$fold_changes)))
    expect_named(sets(), c(res$studies, "Meta"))

    session$setInputs(meta_table_rows_selected = 4)
    expect_equal(selected_gene(), rownames(res$meta)[4])
    expect_s3_class(counts_plot(), "ggplot")
    expect_equal(output$counts_message, "")
  })
})

test_that("GO and KEGG enrichment get the meta-analysis DEGs and the organism", {
  res <- demo_result()
  seen <- list()
  local_mocked_bindings(
    run_go_enrichment = function(genes, organism) {
      seen$go <<- list(genes = genes, organism = organism)
      data.frame()
    },
    run_kegg_enrichment = function(genes, organism) {
      seen$kegg <<- list(genes = genes, organism = organism)
      data.frame()
    }
  )
  shiny::testServer(mod_results_server,
                    args = list(result = shiny::reactive(res), organism = list(go = 1111708, kegg = "syn")), {
    session$setInputs(run_go = 1)
    go_table()
    session$setInputs(run_kegg = 1)
    kegg_table()
  })
  expect_equal(seen$go, list(genes = res$meta_degs, organism = 1111708))
  expect_equal(seen$kegg, list(genes = res$meta_degs, organism = "syn"))
})

test_that("enrichment is not run when there are no meta-analysis DEGs", {
  res <- demo_result()
  res$meta_degs <- character(0)
  called <- FALSE
  local_mocked_bindings(run_go_enrichment = function(...) {
    called <<- TRUE
    data.frame()
  })
  shiny::testServer(mod_results_server,
                    args = list(result = shiny::reactive(res), organism = list(go = 1, kegg = "syn")), {
    session$setInputs(run_go = 1)
    expect_error(go_table(), class = "shiny.silent.error")
  })
  expect_false(called)
})

test_that("the app computes the demo once, when the Demo tab is first opened", {
  res <- demo_result()
  runs <- 0
  local_mocked_bindings(run_tma_pipeline = function(counts, design, ...) {
    runs <<- runs + 1
    res
  })
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(main_tabs = "home")
    expect_equal(runs, 0)
    session$setInputs(main_tabs = "demo")
    expect_s3_class(demo_result(), "tma_result")
    session$setInputs(main_tabs = "rnaseq")
    session$setInputs(main_tabs = "demo")
    demo_result()
    expect_equal(runs, 1)
  })
})

test_that("the RNA-seq tab validates the uploads and runs the pipeline on click", {
  res <- demo_result()
  # the mock runs inside the background task, so it reports its input through
  # the returned result
  local_mocked_bindings(run_tma_pipeline = function(counts, design, ...) {
    out <- res
    out$design <- design
    out
  })
  bad_design <- withr::local_tempfile(fileext = ".txt")
  utils::write.table(tma_demo_data()$design[, c("sample", "condition")], bad_design,
                     sep = "\t", quote = FALSE, row.names = FALSE)
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(rnaseq_counts = list(datapath = tma_example_file("counts")),
                      rnaseq_design = list(datapath = bad_design),
                      rnaseq_analyze = 1)
    expect_error(rnaseq_result(), class = "shiny.silent.error")
    expect_equal(rnaseq_task$status(), "initial")

    session$setInputs(rnaseq_design = list(datapath = tma_example_file("design")),
                      rnaseq_analyze = 2)
    # the analysis runs as an ExtendedTask; wait for its promise to resolve
    for (i in 1:200) {
      later::run_now(0.05)
      session$flushReact()
      if (rnaseq_task$status() != "running") break
    }
    expect_equal(rnaseq_task$status(), "success")
    expect_s3_class(rnaseq_result(), "tma_result")
    expect_equal(rnaseq_result()$design, tma_demo_data()$design)
  })
})
