#' Analysis settings module
#'
#' A box to choose the differential expression method (for all studies, with
#' an optional override per study) and the meta-analysis method, with short
#' guidance on when to use each. With `per_study = TRUE` (the RNA-seq tab) it
#' also offers the batch-effect options; they need a `batch` column in the
#' design, which is checked when the analysis starts (see
#' [validate_batch_choice()]).
#'
#' @param id Module id.
#' @param per_study If `TRUE`, show the reference condition, a DE-method
#'   override for each study of the design passed to `mod_settings_server()`
#'   (RNA-seq only), and the batch-effect options.
#' @param platform `"rnaseq"` or `"array"`. For microarray data the DE method
#'   is always limma, and the batch options are the model term or none.
#' @param design A reactive returning the experimental design, used to list
#'   the studies for the per-study override (`NULL` if `per_study = FALSE`).
#' @return `mod_settings_ui()` returns UI. `mod_settings_server()` returns a
#'   reactive list with `de_method` (one method, or a named vector with one
#'   method per study), `meta_method`, `batch_method` (see [batch_methods()]),
#'   `batch_plots` and `reference` (the reference condition, `NULL` until a
#'   design is available).
#' @name mod_settings
#' @examples
#' # the UI part; the server part runs inside a Shiny app
#' ui <- mod_settings_ui("settings", per_study = TRUE)
#' if (interactive()) {
#'   shiny::shinyApp(shinydashboard::dashboardPage(
#'     shinydashboard::dashboardHeader(), shinydashboard::dashboardSidebar(),
#'     shinydashboard::dashboardBody(ui)), function(input, output, session) NULL)
#' }
#' @export
mod_settings_ui <- function(id, per_study = FALSE, platform = c("rnaseq", "array")) {
  ns <- shiny::NS(id)
  platform <- match.arg(platform)
  de_input <- if (platform == "rnaseq") {
    shiny::selectInput(ns("de_method"), "Differential expression method",
                       choices = de_methods(), selected = "DESeq2")
  } else {
    shiny::div(shiny::strong("Differential expression method"),
               shiny::p("limma (lmFit, contrasts.fit, eBayes with an intensity trend)"))
  }
  shinydashboard::box(
    title = "Analysis settings", status = "info", solidHeader = TRUE, width = 12,
    collapsible = TRUE,
    shiny::fluidRow(
      shiny::column(4, de_input, if (per_study) shiny::uiOutput(ns("reference_ui"))),
      shiny::column(4, shiny::radioButtons(ns("meta_method"), "Meta-analysis method",
                                           choices = meta_methods(), selected = "fisher")),
      if (per_study && platform == "rnaseq") shiny::column(4, shiny::uiOutput(ns("per_study")))
    ),
    if (per_study) batch_inputs(ns, platform),
    method_guidance(platform)
  )
}

#' @rdname mod_settings
#' @export
mod_settings_server <- function(id, design = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    # the design may not be available (yet); the settings are still usable
    current_design <- shiny::reactive({
      if (is.null(design)) return(NULL)
      tryCatch(design(), error = function(e) NULL)
    })
    studies <- shiny::reactive({
      d <- current_design()
      if (is.null(d)) character(0) else unique(d$study)
    })
    study_input <- function(study) paste0("study_method_", make.names(study))
    output$reference_ui <- shiny::renderUI({
      shiny::req(current_design())
      conditions <- sort(unique(as.character(current_design()$condition)))
      shiny::req(length(conditions) >= 2)
      shiny::selectInput(session$ns("reference"),
                         "Reference condition (every other condition is compared with it)",
                         choices = conditions, selected = default_reference(conditions))
    })

    output$per_study <- shiny::renderUI({
      shiny::req(length(studies()) > 0)
      shiny::tagList(
        shiny::strong("DE method per study (optional)"),
        lapply(studies(), function(s) {
          shiny::selectInput(session$ns(study_input(s)), label = s,
                             choices = c("Same as above" = "", de_methods()))
        })
      )
    })

    shiny::reactive({
      # fall back to the defaults until the inputs have been set
      global <- if (is.null(input$de_method)) "DESeq2" else input$de_method
      meta_method <- if (is.null(input$meta_method)) "fisher" else input$meta_method
      de_method <- global
      if (length(studies()) > 0) {
        per_study <- vapply(studies(), function(s) {
          choice <- input[[study_input(s)]]
          if (is.null(choice) || choice == "") global else choice
        }, character(1))
        if (any(per_study != global)) de_method <- stats::setNames(per_study, studies())
      }
      list(de_method = de_method, meta_method = meta_method,
           batch_method = if (is.null(input$batch_method)) "none" else input$batch_method,
           batch_plots = isTRUE(input$batch_plots),
           reference = input$reference)
    })
  })
}

# Batch-effect inputs; "Ignore batch" is selected by default because most
# designs have no batch column
batch_inputs <- function(ns, platform = "rnaseq") {
  choices <- batch_methods()
  if (platform == "array") choices <- choices[choices != "combat_seq"]
  shiny::fluidRow(
    shiny::column(6, shiny::radioButtons(ns("batch_method"),
                                         "Batch effects (needs a 'batch' column in the design file)",
                                         choices = choices, selected = "none")),
    shiny::column(6, shiny::checkboxInput(ns("batch_plots"),
                                          "Use batch-corrected values in the PCA and heatmap",
                                          value = FALSE),
                  shiny::helpText("Batch is handled within each study. The Batch effects tab of the results shows PCA and RLE plots before and after the correction."))
  )
}

# Short in-app guidance on choosing the methods
method_guidance <- function(platform = "rnaseq") {
  array <- platform == "array"
  shiny::tags$details(
    # styled as a link so it is clear that it can be clicked to open
    shiny::tags$summary(
      style = "color: #337ab7; text-decoration: underline; cursor: pointer; font-weight: bold;",
      "Which method should I use?"
    ),
    shiny::tags$p(shiny::strong("Differential expression"), " (run in each study):"),
    if (array) shiny::tags$ul(
      shiny::tags$li(shiny::strong("limma"), ": linear models on the log2 expression values, every condition compared with the reference through contrasts, and empirical Bayes moderation with an intensity trend (eBayes(trend = TRUE)), because the variance of microarray data depends on the expression level (Ritchie et al., 2015, Nucleic Acids Res 43:e47).")
    ),
    if (!array) shiny::tags$ul(
      shiny::tags$li(shiny::strong("DESeq2"), " (default): negative-binomial model with shrunken dispersion estimates. Works well with few replicates per group."),
      shiny::tags$li(shiny::strong("edgeR"), " (quasi-likelihood F-test): also a negative-binomial model; the quasi-likelihood test keeps the false discovery rate under control with few replicates."),
      shiny::tags$li(shiny::strong("limma-voom"), ": linear models on log-counts with precision weights. Fast, and a good choice for larger studies (several replicates per group) or when library sizes differ a lot.")
    ),
    if (!array) shiny::tags$p("In a benchmark with few replicates, DESeq2 and edgeR were among the best-performing tools (Schurch et al., 2016, RNA 22:839-851). Use the per-study setting if one study differs from the others, e.g. it has many more samples."),
    shiny::tags$p(shiny::strong("Meta-analysis"), " (combines the p-values of the studies):"),
    shiny::tags$ul(
      shiny::tags$li(shiny::strong("Fisher"), ": sums -log(p) over the studies. A very strong signal in a single study can be enough to make a gene significant."),
      shiny::tags$li(shiny::strong("Inverse normal"), ": combines the studies' z-scores, weighting each study by the square root of its number of samples. It favours genes with consistent evidence across studies and gives larger studies more weight (Rau et al., 2014, BMC Bioinformatics 15:91).")
    ),
    shiny::tags$p(shiny::strong("Batch effects"), " (needs a 'batch' column in the design file):"),
    if (array) shiny::tags$ul(
      shiny::tags$li(shiny::strong("Batch in the DE model"), ": adds batch to the limma model of each study, so the p-values account for it."),
      shiny::tags$li(shiny::strong("Ignore batch"), " (default): the batch column, if any, is not used."),
      shiny::tags$li("Batch can only be used in a study with at least two batches that are not confounded with the condition.")
    ),
    if (!array) shiny::tags$ul(
      shiny::tags$li(shiny::strong("Batch in the DE model"), " (recommended when there is a batch column): fits ~ batch + condition in each study. The counts are not changed; the batch effect is estimated together with the condition effect, so the p-values account for it."),
      shiny::tags$li(shiny::strong("ComBat-seq"), ": corrects the raw counts for batch first (keeping the condition effect), then fits ~ condition. The DE test then treats the corrected counts as raw data, which can make p-values too small in small studies (Nygaard et al., 2016, Biostatistics 17:29-39); it also needs at least two samples per batch."),
      shiny::tags$li(shiny::strong("Ignore batch"), " (default): the batch column, if any, is not used."),
      shiny::tags$li("Batch can only be used in a study with at least two batches that are not confounded with the condition; otherwise that study is analysed without batch, and the Batch effects tab says why.")
    ),
    shiny::tags$p(if (array) "Click Start Analysis! again after changing a setting." else
      "On the RNA-Seq tab, click Start Analysis! again after changing a setting; the Demo tab updates by itself."),
    shiny::tags$p("More details and references are in the Methods Guide tab.")
  )
}
