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
    expect_equal(fname("pca.pdf"), paste0(session$ns(NULL), "_pca.pdf"))
    expect_equal(dim(heat_csv), c(length(res$meta_degs), nrow(res$design)))

    expect_equal(volcano()$fc, unname(mean_fold_change(res$fold_changes)))
    expect_named(sets(), c(res$studies, "Meta"))

    session$setInputs(meta_table_rows_selected = 4)
    expect_equal(selected_gene(), rownames(res$meta)[4])
    expect_s3_class(counts_plot(), "ggplot")
    expect_equal(output$counts_message, "")
  })
})

fake_enrichment <- function(service, genes, ...) {
  structure(list(service = service,
                 table = data.frame(term.id = "GO:1", term.label = "t", fold_enrichment = 2, fdr = 0.01,
                                    number_in_list = 3, plus_minus = "+", category = "BP"),
                 mapping = mapping_report(genes, genes[-1], service),
                 version = "test version", background = "all tested genes",
                 organism = "test organism", id_type = "kegg"),
            class = "tma_enrichment")
}

test_that("GO, KEGG and STRING get the DEGs, organism, ID type and tested-genes background", {
  res <- demo_result()
  seen <- list()
  record <- function(service) function(genes, organism = NULL, species = NULL, background, id_type, ...) {
    seen[[service]] <<- list(genes = genes, organism = c(organism, species),
                             background = background, id_type = id_type)
    fake_enrichment(service, genes)
  }
  local_mocked_bindings(run_go_enrichment = record("PANTHER"), run_kegg_enrichment = record("KEGG"),
                        run_string_enrichment = record("STRING"))
  fixed <- list(go = 1111708, kegg = "syn", string = 1148, id_type = "kegg")
  shiny::testServer(mod_results_server, args = list(result = shiny::reactive(res), organism = fixed), {
    session$setInputs(run_go = 1, run_kegg = 1, run_string = 1)
    expect_s3_class(go_result(), "tma_enrichment")
    kegg_result()
    string_result()
    report <- as.character(output$go_report$html)
    expect_match(report, "Mapping rate")
    expect_match(report, "Download the IDs that were not mapped")
    expect_match(report, "test version")
    table_file <- output$download_go_table
    expect_match(readLines(table_file, n = 1), "^# PANTHER \\(test version\\)")
    expect_equal(readLines(output$download_go_unmapped), res$meta_degs[1])
  })
  for (service in c("PANTHER", "KEGG", "STRING")) {
    expect_equal(seen[[service]]$genes, res$meta_degs)
    expect_equal(seen[[service]]$background, res$genes)
    expect_equal(seen[[service]]$id_type, "kegg")
  }
  expect_equal(seen$PANTHER$organism, 1111708)
  expect_equal(seen$KEGG$organism, "syn")
  expect_equal(seen$STRING$organism, 1148)
})

test_that("a failed enrichment is shown as a message, not an error", {
  res <- demo_result()
  alerts <- character(0)
  local_mocked_bindings(run_string_enrichment = function(...) stop("Timeout was reached"))
  local_mocked_bindings(shinyalert = function(title, text = "", ...) alerts <<- c(alerts, text),
                        .package = "shinyalert")
  shiny::testServer(mod_results_server,
                    args = list(result = shiny::reactive(res),
                                organism = list(go = 1, kegg = "syn", string = 1148, id_type = "kegg")), {
    session$setInputs(run_string = 1)
    expect_error(string_result(), class = "shiny.silent.error")
  })
  expect_match(alerts, "The STRING enrichment could not be completed \\(Timeout was reached\\)", all = FALSE)
})

test_that("enrichment is not run without DEGs, ID type or organism", {
  res <- demo_result()
  called <- FALSE
  local_mocked_bindings(run_go_enrichment = function(...) {
    called <<- TRUE
    fake_enrichment("PANTHER", "g")
  })
  mock_organism_lists()
  no_degs <- res
  no_degs$meta_degs <- character(0)
  shiny::testServer(mod_results_server,
                    args = list(result = shiny::reactive(no_degs),
                                organism = list(go = 1, kegg = "syn", string = 1148, id_type = "kegg")), {
    session$setInputs(run_go = 1)
    expect_error(go_result(), class = "shiny.silent.error")
  })
  # RNA-seq tab: the user must choose the ID type and the organism
  shiny::testServer(mod_results_server, args = list(result = shiny::reactive(res)), {
    session$setInputs(organism_go = "1111708", run_go = 1)
    expect_error(go_result(), class = "shiny.silent.error")
    session$setInputs(id_type = "kegg", organism_go = "", run_go = 2)
    expect_error(go_result(), class = "shiny.silent.error")
  })
  expect_false(called)
})

test_that("the KEGG organism code is checked while typing", {
  mock_organism_lists()
  shiny::testServer(mod_results_server, args = list(result = shiny::reactive(demo_result())), {
    session$setInputs(organism_kegg = "syn")
    expect_equal(output$kegg_organism_check, "Organism: Synechocystis sp. PCC 6803")
    session$setInputs(organism_kegg = "xyz9")
    expect_equal(output$kegg_organism_check, "Unknown KEGG organism code.")
  })
})

test_that("the app computes the demo once, when the Demo tab is first opened", {
  res <- demo_result()
  runs <- 0
  local_mocked_bindings(run_tma_pipeline = function(counts, design, ...) {
    runs <<- runs + 1
    res
  })
  mock_organism_lists()
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
  local_mocked_bindings(run_tma_analysis = function(counts, design, ...) {
    out <- res
    out$design <- design
    fake_analysis(out)
  })
  bad_design <- withr::local_tempfile(fileext = ".txt")
  utils::write.table(tma_demo_data()$design[, c("sample", "condition")], bad_design,
                     sep = "\t", quote = FALSE, row.names = FALSE)
  mock_organism_lists()
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(rnaseq_counts = list(datapath = tma_example_file("counts")),
                      rnaseq_design = list(datapath = bad_design),
                      rnaseq_analyze = 1)
    expect_error(rnaseq_analysis(), class = "shiny.silent.error")
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
    expect_s3_class(rnaseq_analysis()$results[[1]], "tma_result")
    expect_equal(rnaseq_analysis()$results[[1]]$design, tma_demo_data()$design)
  })
})

test_that("the meta-analysis table shows 4 significant digits but keeps the data", {
  res <- demo_result()
  widget <- meta_table_widget(res$meta)
  expect_s3_class(widget, "datatables")
  defs <- widget$x$options$columnDefs
  rendered <- Filter(function(d) !is.null(d$render), defs)
  # columns are counted after the row-name column: TestStatistic, rawpval, adjpval
  expect_setequal(unlist(lapply(rendered, `[[`, "targets")), 1:3)
  for (d in rendered) {
    # 4 significant digits on display; sorting and searching use the full value
    expect_match(as.character(d$render), "formatSignif\\(data, 4,")
    expect_match(as.character(d$render), "type !== 'display' \\? data")
  }
  expect_equal(widget$x$data$adjpval, res$meta$adjpval)
})

test_that("the Batch effects tab is hidden for a result without a batch column", {
  # (the browser-only part, running the update outside the reactive domain,
  # cannot be reproduced in testServer)
  calls <- list()
  local_mocked_bindings(
    hideTab = function(inputId, target, session = shiny::getDefaultReactiveDomain()) {
      calls[[length(calls) + 1]] <<- list(target = target)
    },
    .package = "shiny"
  )
  res <- demo_result() # no batch column
  shiny::testServer(mod_results_server, args = list(result = shiny::reactive(res),
                                                    organism = list(go = 1, kegg = "syn")), {
    session$flushReact()
  })
  hidden <- vapply(calls, `[[`, "", "target")
  expect_true("batch" %in% hidden)
})

test_that("the organism lists are filled again for every new result", {
  # the results UI is drawn again for every analysis, so its organism inputs
  # start empty each time
  calls <- list()
  local_mocked_bindings(
    updateSelectizeInput = function(session, inputId, label = NULL, choices = NULL, selected = NULL, ...) {
      calls[[length(calls) + 1]] <<- list(id = inputId, choices = choices, selected = selected)
    },
    .package = "shiny"
  )
  mock_organism_lists()
  res <- demo_result()
  current <- shiny::reactiveVal(res)
  shiny::testServer(mod_results_server, args = list(result = current), {
    session$flushReact()
    second <- res
    second$contrast <- "High vs Cont"
    current(second)
    session$flushReact()
    # the updates are sent after the flush in which the result arrived
    session$flushReact()
  })
  go_updates <- Filter(function(x) x$id == "organism_go" && length(x$choices) > 0, calls)
  expect_length(go_updates, 2)
})
