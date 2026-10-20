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
  #the uploaded design lists the studies for the per-study DE method setting
  uploaded_design <- reactive({
    req(input$rnaseq_design)
    design <- read_tma_table(input$rnaseq_design$datapath)
    req(length(validate_design(design)) == 0)
    design
  })
  rnaseq_settings <- mod_settings_server("rnaseq_settings", design = uploaded_design)

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

  #check the uploaded data for batch effects before the analysis
  mod_batch_check_server("rnaseq_batch_check", data = read_uploads)

  #read and check the uploaded files when "Start Analysis!" is hit
  rnaseq_inputs <- eventReactive(input$rnaseq_analyze, {
    uploads <- read_uploads()
    settings <- rnaseq_settings()
    problems <- validate_batch_choice(uploads$design, settings$batch_method, settings$batch_plots)
    if (length(problems) > 0) {
      shinyalert(text = problems[1])
      req(FALSE)
    }
    c(uploads, settings)
  })

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
  #disable "Start Analysis!" while an analysis is running
  observe({
    shinyjs::toggleState("rnaseq_analyze", condition = rnaseq_task$status() != "running")
  })
  observeEvent(rnaseq_inputs(), {
    inputs <- rnaseq_inputs()
    shinyalert("Running...", "RNA-Seq meta-analysis is in progress... \nResults will appear soon!",
               type = "info", closeOnEsc = TRUE, closeOnClickOutside = TRUE, timer = 7000)
    shinyjs::hide("rnaseq_steps")
    rnaseq_task$invoke(inputs$counts, inputs$design, inputs$de_method, inputs$meta_method,
                       inputs$batch_method, inputs$batch_plots, inputs$reference)
  })
  rnaseq_analysis <- reactive({
    #an error in the analysis is shown to the user instead of stopping the app
    tryCatch(rnaseq_task$result(), error = function(e) {
      if (inherits(e, "shiny.silent.error")) stop(e)
      shinyalert("Analysis failed", conditionMessage(e), type = "error")
      req(FALSE)
    })
  })

  #one results tab per contrast (condition vs reference), and a summary
  output$rnaseq_results <- renderUI({
    analysis <- rnaseq_analysis()
    results <- analysis$results
    messages <- analysis_messages(analysis)
    tabs <- lapply(seq_along(results), function(i) {
      tabPanel(names(results)[i],
               mod_results_ui(paste0("rnaseq_c", i), select_organism = TRUE,
                              venn = overlap_plot_type(deg_sets(results[[i]])) == "venn",
                              batch = !is.null(results[[i]]$batch)))
    })
    tagList(
      fluidRow(box(
        title = "Comparisons", status = "primary", solidHeader = TRUE, width = 12,
        p(sprintf("Every condition is compared with the reference condition '%s'. Each comparison that could be meta-analysed has its own results tab below.",
                  analysis$reference)),
        tableOutput("rnaseq_summary"),
        if (length(messages) > 0) tags$ul(lapply(messages, tags$li))
      )),
      do.call(tabsetPanel, c(list(id = "rnaseq_contrasts"), unname(tabs)))
    )
  })
  output$rnaseq_summary <- renderTable({
    summary <- rnaseq_analysis()$summary
    data.frame(Comparison = summary$contrast, `Studies used` = summary$studies,
               `Meta-analysed` = ifelse(summary$analysed, "yes", "no"), check.names = FALSE)
  })
  #a results module per contrast tab, created once per tab position
  n_result_modules <- 0
  observeEvent(rnaseq_analysis(), {
    analysis <- rnaseq_analysis()
    while (n_result_modules < length(analysis$results)) {
      n_result_modules <<- n_result_modules + 1
      local({
        i <- n_result_modules
        result_i <- reactive({
          results <- rnaseq_analysis()$results
          req(length(results) >= i)
          results[[i]]
        })
        mod_results_server(paste0("rnaseq_c", i), result_i,
                           file_prefix = reactive(paste0("rnaseq_", gsub("[^A-Za-z0-9]+", "_", result_i()$contrast))))
      })
    }
    #tell the user which studies or comparisons were left out (replicate rule)
    messages <- analysis_messages(analysis)
    if (length(messages) > 0) {
      shinyalert("Some studies were left out", paste(messages, collapse = "\n\n"), type = "warning")
    }
  })

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
  mod_results_server("demo", demo_result, organism = list(go = 1111708, kegg = "syn"))

  #stop the app after the browser closes, but only for local runs started with
  #runTMA(); on a shared server this would shut the app down for every user
  if (isTRUE(getOption("tma.stop_on_session_end", FALSE))) {
    session$onSessionEnded(function() {
      stopApp()
    })
  }
}
