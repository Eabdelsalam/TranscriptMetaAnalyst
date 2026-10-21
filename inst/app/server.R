function(input, output, session) {
  #increase maximum file upload size to 30 MB
  options(shiny.maxRequestSize = 30 * 1024^2)

  #download the example files
  output$example_counts <- downloadHandler(
    filename = "example_rnaseq_counts.txt",
    content = function(file) file.copy(tma_example_file("counts"), file),
    contentType = "text/plain"
  )
  output$example_design <- downloadHandler(
    filename = "example_rnaseq_design.txt",
    content = function(file) file.copy(tma_example_file("design"), file),
    contentType = "text/plain"
  )

  ####### RNA-Seq Analysis #######
  #the data that "Start Analysis!" analyses: the uploaded files, or the demo
  #data after "Load demo data" (until a file is uploaded again)
  rnaseq_source <- reactiveVal("upload")
  observeEvent(list(input$rnaseq_counts, input$rnaseq_design), rnaseq_source("upload"),
               ignoreInit = TRUE)

  #the current design lists the studies and conditions for the settings
  rnaseq_design <- reactive({
    if (rnaseq_source() == "demo") return(tma_demo_data()$design)
    req(input$rnaseq_design)
    design <- read_tma_table(input$rnaseq_design$datapath)
    req(length(validate_design(design)) == 0)
    design
  })
  rnaseq_settings <- mod_settings_server("rnaseq_settings", design = rnaseq_design)

  #read and check the uploaded files; shows an alert and stops on a problem
  read_uploads <- function() {
    if (is.null(input$rnaseq_counts)) {
      shinyalert("Missed Counts File!!!", "Please upload the raw counts file first!")
    } else if (is.null(input$rnaseq_design)) {
      shinyalert("Missed Experimental Desgin File!!!", "Please upload the experimental design file first!")
    }
    req(input$rnaseq_counts, input$rnaseq_design)
    counts <- read_tma_table(input$rnaseq_counts$datapath)
    design <- read_tma_table(input$rnaseq_design$datapath)
    problems <- c(validate_counts(counts), validate_design(design))
    if (length(problems) == 0) problems <- validate_match(counts, design)
    if (length(problems) > 0) {
      shinyalert(text = problems[1])
      req(FALSE)
    }
    list(counts = counts, design = design)
  }
  rnaseq_data <- function() if (rnaseq_source() == "demo") tma_demo_data() else read_uploads()

  #check the current data for batch effects before the analysis
  mod_batch_check_server("rnaseq_batch_check", data = rnaseq_data)

  #run the meta-analysis once per click, in a background process so the app
  #stays responsive while it runs
  rnaseq_task <- ExtendedTask$new(function(counts, design, de_method, meta_method,
                                            batch_method, batch_plots, reference) {
    promises::future_promise(
      run_tma_analysis(counts, design, reference = reference, de_method = de_method,
                       meta_method = meta_method, batch_method = batch_method,
                       batch_plots = batch_plots),
      seed = TRUE
    )
  })
  #disable the buttons while an analysis is running
  observe({
    shinyjs::toggleState("rnaseq_analyze", condition = rnaseq_task$status() != "running")
    shinyjs::toggleState("rnaseq_demo", condition = rnaseq_task$status() != "running")
  })
  #what the last run used, for the R script and the enrichment choices
  rnaseq_run <- reactiveVal(NULL)
  start_rnaseq_analysis <- function(data, settings) {
    problems <- validate_batch_choice(data$design, settings$batch_method, settings$batch_plots)
    if (length(problems) > 0) {
      shinyalert(text = problems[1])
      req(FALSE)
    }
    shinyalert("Running...", "RNA-Seq meta-analysis is in progress... \nResults will appear soon!",
               type = "info", closeOnEsc = TRUE, closeOnClickOutside = TRUE, timer = 7000)
    rnaseq_run(list(data = rnaseq_source(), settings = settings,
                    counts_name = input$rnaseq_counts$name, design_name = input$rnaseq_design$name))
    rnaseq_task$invoke(data$counts, data$design, settings$de_method, settings$meta_method,
                       settings$batch_method, settings$batch_plots, settings$reference)
  }
  #the demo data have no batch column, and the reference is guessed unless one
  #of their conditions is chosen
  observeEvent(input$rnaseq_demo, {
    rnaseq_source("demo")
    demo <- tma_demo_data()
    settings <- rnaseq_settings()
    settings$batch_method <- "none"
    settings$batch_plots <- FALSE
    if (!isTRUE(settings$reference %in% demo$design$condition)) settings$reference <- NULL
    start_rnaseq_analysis(demo, settings)
  })
  observeEvent(input$rnaseq_analyze, {
    start_rnaseq_analysis(rnaseq_data(), rnaseq_settings())
  })
  rnaseq_analysis <- reactive({
    #an error in the analysis is shown to the user instead of stopping the app
    tryCatch(rnaseq_task$result(), error = function(e) {
      if (inherits(e, "shiny.silent.error")) stop(e)
      shinyalert("Analysis failed", conditionMessage(e), type = "error")
      req(FALSE)
    })
  })

  #results of an analysis: a summary of the comparisons and one results tab per
  #contrast (condition vs reference); used by the RNA-seq and microarray tabs
  analysis_results <- function(prefix, analysis, script_id, selected = NULL) {
    output[[paste0(prefix, "_results")]] <- renderUI({
      res <- analysis()
      results <- res$results
      messages <- analysis_messages(res)
      tabs <- lapply(seq_along(results), function(i) {
        tabPanel(names(results)[i],
                 mod_results_ui(paste0(prefix, "_c", i), select_organism = TRUE,
                                venn = overlap_plot_type(deg_sets(results[[i]])) == "venn",
                                batch = !is.null(results[[i]]$batch)))
      })
      tagList(
        fluidRow(box(
          title = "Comparisons", status = "primary", solidHeader = TRUE, width = 12,
          p(sprintf("Every condition is compared with the reference condition '%s'. Each comparison that could be meta-analysed has its own results tab below.",
                    res$reference)),
          tableOutput(paste0(prefix, "_summary")),
          if (length(messages) > 0) tags$ul(lapply(messages, tags$li)),
          downloadButton(script_id, "Download R script"),
          helpText("An R script that repeats this analysis with the TranscriptMetaAnalyst package. For uploaded files, the input file paths in the script are placeholders: replace them with the paths of your files before running it.")
        )),
        do.call(tabsetPanel, c(list(id = paste0(prefix, "_contrasts")), unname(tabs)))
      )
    })
    output[[paste0(prefix, "_summary")]] <- renderTable({
      summary <- analysis()$summary
      data.frame(Comparison = summary$contrast, `Studies used` = summary$studies,
                 `Meta-analysed` = ifelse(summary$analysed, "yes", "no"), check.names = FALSE)
    })
    #a results module per contrast tab, created once per tab position
    n_modules <- 0
    observeEvent(analysis(), {
      res <- analysis()
      while (n_modules < length(res$results)) {
        n_modules <<- n_modules + 1
        local({
          i <- n_modules
          result_i <- reactive({
            results <- analysis()$results
            req(length(results) >= i)
            results[[i]]
          })
          #the tabs that do not apply are left out of the UI drawn for the result
          mod_results_server(paste0(prefix, "_c", i), result_i,
                             file_prefix = reactive(paste0(prefix, "_", gsub("[^A-Za-z0-9]+", "_", result_i()$contrast))),
                             toggle_tabs = FALSE, selected = selected)
        })
      }
      #tell the user which studies or comparisons were left out (replicate rule)
      messages <- analysis_messages(res)
      if (length(messages) > 0) {
        shinyalert("Some studies were left out", paste(messages, collapse = "\n\n"), type = "warning")
      }
    })
  }
  #the demo data use Synechocystis locus tags: preselect them for the enrichment
  demo_enrichment <- list(id_type = "kegg", go = 1111708, kegg = "syn", string = 1148)
  analysis_results("rnaseq", rnaseq_analysis, "rnaseq_script",
                   selected = reactive(if (identical(rnaseq_run()$data, "demo")) demo_enrichment))

  #reproducible R script of the last analysis (placeholder file paths for
  #uploaded files)
  output$rnaseq_script <- downloadHandler(
    filename = "tma_rnaseq_analysis.R",
    content = function(file) {
      run <- rnaseq_run()
      settings <- run$settings[c("de_method", "meta_method", "batch_method", "batch_plots")]
      settings$reference <- rnaseq_analysis()$reference
      #enrichment choices of the first results tab, if made
      choice <- function(name) input[[paste0("rnaseq_c1-", name)]]
      write_tma_script(file, settings = settings, data = run$data,
                       counts_name = run$counts_name, design_name = run$design_name,
                       enrichment = list(id_type = choice("id_type"), go = choice("organism_go"),
                                         kegg = choice("organism_kegg"), string = choice("organism_string")))
    },
    contentType = "text/plain"
  )

  ####### Microarray Analysis #######
  #GEO series fetched in this session: a list (by GSE) of ExpressionSets by platform
  geo_sets <- reactiveVal(list())
  observeEvent(input$geo_fetch, {
    ids <- unique(toupper(trimws(unlist(strsplit(input$geo_ids, "[,; ]+")))))
    ids <- ids[nzchar(ids)]
    if (length(ids) == 0) {
      shinyalert("No accession", "Please enter at least one GEO series accession (e.g. GSE82120).")
      return()
    }
    fetched <- list()
    withProgress(message = "Downloading from GEO", value = 0, {
      for (gse in ids) {
        incProgress(1 / length(ids), detail = gse)
        sets <- tryCatch(fetch_geo(gse), error = function(e) {
          shinyalert("GEO import", enrichment_message(e, "GEO"), type = "error")
          NULL
        })
        if (!is.null(sets)) fetched[[gse]] <- sets
      }
    })
    geo_sets(fetched)
  })
  #the chosen platform of each fetched series
  geo_eset <- function(i) {
    sets <- geo_sets()[[i]]
    platform <- input[[paste0("geo_platform_", i)]]
    if (is.null(platform) || !platform %in% names(sets)) platform <- names(sets)[1]
    sets[[platform]]
  }
  #one panel per series: platform, condition column, condition labels, gene column
  output$geo_config <- renderUI({
    sets <- geo_sets()
    req(length(sets) > 0)
    tagList(
      helpText("For each series, choose the column that holds the condition and give each of its values a condition label. Samples whose value has no label are not used. Use the same labels in all series (e.g. 'Control', 'Treated')."),
      lapply(seq_along(sets), function(i) {
        wellPanel(
          h4(names(sets)[i]),
          if (length(sets[[i]]) > 1) {
            selectInput(paste0("geo_platform_", i), "Platform (this series uses several)", choices = names(sets[[i]]))
          },
          uiOutput(paste0("geo_detail_", i))
        )
      })
    )
  })
  observeEvent(geo_sets(), {
    for (i in seq_along(geo_sets())) {
      local({
        idx <- i
        output[[paste0("geo_detail_", idx)]] <- renderUI({
          es <- geo_eset(idx)
          tagList(
            if (geo_is_two_colour(es)) {
              div(class = "alert alert-warning",
                  "This looks like a two-colour experiment: its channels are not independent samples and should not be meta-analysed as single-channel data.")
            },
            fluidRow(
              column(6, selectInput(paste0("geo_cond_", idx), "Condition column", choices = geo_condition_columns(es))),
              column(6, selectInput(paste0("geo_gene_", idx), "Gene column",
                                    choices = c("None (use the probe IDs)" = "", geo_gene_columns(es))))
            ),
            uiOutput(paste0("geo_labels_", idx))
          )
        })
        output[[paste0("geo_labels_", idx)]] <- renderUI({
          es <- geo_eset(idx)
          col <- input[[paste0("geo_cond_", idx)]]
          req(col, col %in% names(Biobase::pData(es)))
          values <- unique(as.character(Biobase::pData(es)[[col]]))
          counts <- table(as.character(Biobase::pData(es)[[col]]))[values]
          tagList(lapply(seq_along(values), function(j) {
            textInput(paste0("geo_label_", idx, "_", j),
                      sprintf("%s (%d samples)", values[j], counts[[j]]), value = "",
                      placeholder = "condition label, or empty to leave out")
          }))
        })
      })
    }
  })
  #the import choices of each series, and the studies built from them
  geo_choices <- reactive({
    sets <- geo_sets()
    req(length(sets) > 0)
    lapply(seq_along(sets), function(i) {
      es <- geo_eset(i)
      col <- input[[paste0("geo_cond_", i)]]
      req(col)
      values <- unique(as.character(Biobase::pData(es)[[col]]))
      labels <- vapply(seq_along(values), function(j) {
        x <- input[[paste0("geo_label_", i, "_", j)]]
        if (is.null(x)) "" else trimws(x)
      }, character(1))
      list(gse = names(sets)[i], platform = Biobase::annotation(es), condition_column = col,
           labels = stats::setNames(labels, values), gene_column = input[[paste0("geo_gene_", i)]])
    })
  })
  geo_studies <- reactive({
    choices <- geo_choices()
    lapply(choices, function(g) {
      geo_study(geo_sets()[[g$gse]][[g$platform]], g$gse, g$condition_column, g$labels,
                if (is.null(g$gene_column) || g$gene_column == "") NULL else g$gene_column)
    })
  })

  #the design of the current input, for the reference-condition list
  array_design_input <- reactive({
    if (identical(input$array_source, "geo")) {
      design <- tryCatch(do.call(rbind, lapply(geo_studies(), `[[`, "design")), error = function(e) NULL)
      req(design)
      return(design)
    }
    req(input$array_design)
    design <- read_tma_table(input$array_design$datapath)
    req(length(validate_design(design)) == 0)
    design
  })
  array_settings <- mod_settings_server("array_settings", design = array_design_input)

  #what the last run used, for the R script
  array_run <- reactiveVal(NULL)
  array_task <- ExtendedTask$new(function(matrices, design, reference, meta_method, batch_method, batch_plots) {
    promises::future_promise(
      run_array_analysis(matrices, design, reference = reference, meta_method = meta_method,
                         batch_method = batch_method, batch_plots = batch_plots),
      seed = TRUE
    )
  })
  observe({
    shinyjs::toggleState("array_analyze", condition = array_task$status() != "running")
    shinyjs::toggleState("array_demo", condition = array_task$status() != "running")
  })
  start_array_analysis <- function(matrices, design, settings, run) {
    shinyalert("Running...", "Microarray meta-analysis is in progress... \nResults will appear soon!",
               type = "info", closeOnEsc = TRUE, closeOnClickOutside = TRUE, timer = 7000)
    array_run(c(run, list(settings = settings)))
    array_task$invoke(matrices, design, settings$reference, settings$meta_method,
                      settings$batch_method, settings$batch_plots)
  }
  observeEvent(input$array_demo, {
    demo <- tma_microarray_demo_data()
    settings <- array_settings()
    settings$batch_method <- "none"
    settings$batch_plots <- FALSE
    settings$reference <- "HC"
    start_array_analysis(demo$matrices, demo$design, settings, list(data = "demo"))
  })
  observeEvent(input$array_analyze, {
    settings <- array_settings()
    if (identical(input$array_source, "geo")) {
      studies <- tryCatch(geo_studies(), error = function(e) {
        shinyalert("GEO import", if (inherits(e, "shiny.silent.error")) "Please fetch the series and choose their condition labels first." else conditionMessage(e))
        NULL
      })
      req(studies)
      matrices <- stats::setNames(lapply(studies, `[[`, "matrix"), vapply(geo_choices(), `[[`, "", "gse"))
      design <- do.call(rbind, unname(lapply(studies, `[[`, "design")))
      run <- list(data = "geo", geo = geo_choices())
    } else {
      if (is.null(input$array_files) || is.null(input$array_design)) {
        shinyalert("Missing files", "Please upload one expression file per study and the experimental design file first!")
        req(FALSE)
      }
      matrices <- lapply(input$array_files$datapath, read_array_matrix)
      names(matrices) <- input$array_files$name
      design <- read_tma_table(input$array_design$datapath)
      run <- list(data = "upload", files = input$array_files$name, design_name = input$array_design$name)
    }
    problems <- c(validate_design(design),
                  validate_batch_choice(design, settings$batch_method, settings$batch_plots))
    if (length(problems) > 0) {
      shinyalert(text = problems[1])
      req(FALSE)
    }
    start_array_analysis(matrices, design, settings, run)
  })
  array_analysis <- reactive({
    tryCatch(array_task$result(), error = function(e) {
      if (inherits(e, "shiny.silent.error")) stop(e)
      shinyalert("Analysis failed", conditionMessage(e), type = "error")
      req(FALSE)
    })
  })
  #what was done to the data before the analysis (log2, probes, common genes)
  output$array_preparation <- renderUI({
    prep <- array_analysis()$preparation
    fluidRow(box(
      title = "Data preparation", status = "info", solidHeader = TRUE, width = 12, collapsible = TRUE,
      tableOutput("array_preparation_table"),
      tags$ul(lapply(prep$notes, tags$li))
    ))
  })
  output$array_preparation_table <- renderTable({
    report <- array_analysis()$preparation$report
    data.frame(Study = report$study, Samples = report$samples, Probes = report$probes,
               `Probes without a gene` = report$no_gene,
               `Probes with several genes` = report$multiple_genes,
               Genes = report$genes, `log2 applied` = ifelse(report$log2_applied, "yes", "no"),
               check.names = FALSE)
  })
  analysis_results("array", array_analysis, "array_script")
  output$array_script <- downloadHandler(
    filename = "tma_microarray_analysis.R",
    content = function(file) {
      run <- array_run()
      settings <- run$settings
      settings$reference <- array_analysis()$reference
      choice <- function(name) input[[paste0("array_c1-", name)]]
      write_tma_script(file, settings = settings, data = run$data, platform = "array",
                       array_files = run$files, design_name = run$design_name, geo = run$geo,
                       enrichment = list(id_type = choice("id_type"), go = choice("organism_go"),
                                         kegg = choice("organism_kegg"), string = choice("organism_string")))
    },
    contentType = "text/plain"
  )

  ####### Demo #######
  #computed the first time the Demo tab is opened with given settings, then
  #cached for all sessions
  demo_settings <- mod_settings_server("demo_settings")
  mod_batch_check_server("demo_batch_check", data = tma_demo_data)
  demo_result <- reactive({
    shinyalert("Demo Running!", "Demo RNA-Seq meta-analysis is in progress... \nResults will appear soon!",
               type = "info", closeOnEsc = TRUE, closeOnClickOutside = TRUE, timer = 7000)
    demo <- tma_demo_data()
    run_tma_pipeline(demo$counts, demo$design,
                     de_method = demo_settings()$de_method,
                     meta_method = demo_settings()$meta_method)
  }) |>
    bindCache("tma-demo", utils::packageVersion("TranscriptMetaAnalyst"), demo_settings(),
              cache = "app") |>
    bindEvent(req(input$main_tabs == "demo"), demo_settings())
  mod_results_server("demo", demo_result,
                     organism = list(go = 1111708, kegg = "syn", string = 1148, id_type = "kegg"))
  output$demo_script <- downloadHandler(
    filename = "tma_demo_analysis.R",
    content = function(file) {
      settings <- c(demo_settings()[c("de_method", "meta_method")],
                    list(batch_method = "none", batch_plots = FALSE, reference = "Cont"))
      write_tma_script(file, settings = settings, data = "demo",
                       enrichment = list(id_type = "kegg", go = 1111708, kegg = "syn", string = 1148))
    },
    contentType = "text/plain"
  )

  #stop the app after the browser closes, but only for local runs started with
  #runTMA(); on a shared server this would shut the app down for every user
  if (isTRUE(getOption("tma.stop_on_session_end", FALSE))) {
    session$onSessionEnded(function() {
      stopApp()
    })
  }
}
