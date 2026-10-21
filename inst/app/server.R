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
  #read and check the uploaded files when "Start Analysis!" is hit
  rnaseq_inputs <- eventReactive(input$rnaseq_analyze, {
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
  })

  #run the meta-analysis once per click, in a background process so the app
  #stays responsive while it runs
  rnaseq_task <- ExtendedTask$new(function(counts, design) {
    promises::future_promise(run_tma_pipeline(counts, design), seed = TRUE)
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
    rnaseq_task$invoke(inputs$counts, inputs$design)
  })
  rnaseq_result <- reactive({
    #an error in the analysis is shown to the user instead of stopping the app
    tryCatch(rnaseq_task$result(), error = function(e) {
      if (inherits(e, "shiny.silent.error")) stop(e)
      shinyalert("Analysis failed", conditionMessage(e), type = "error")
      req(FALSE)
    })
  })
  mod_results_server("rnaseq", rnaseq_result)

  ####### Demo #######
  #computed the first time the Demo tab is opened, then cached for all sessions
  demo_result <- reactive({
    shinyalert("Demo Running!", "Demo RNA-Seq meta-analysis is in progress... \nResults will appear soon!",
               type = "info", closeOnEsc = TRUE, closeOnClickOutside = TRUE, timer = 7000)
    demo <- tma_demo_data()
    run_tma_pipeline(demo$counts, demo$design)
  }) |>
    bindCache("tma-demo", utils::packageVersion("TranscriptMetaAnalyst"), cache = "app") |>
    bindEvent(req(input$main_tabs == "demo"))
  mod_results_server("demo", demo_result, organism = list(go = 1111708, kegg = "syn"))

  #stop the app after the browser closes, but only for local runs started with
  #runTMA(); on a shared server this would shut the app down for every user
  if (isTRUE(getOption("tma.stop_on_session_end", FALSE))) {
    session$onSessionEnded(function() {
      stopApp()
    })
  }
}
