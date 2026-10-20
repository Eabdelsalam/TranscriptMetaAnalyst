#' Analysis settings module
#'
#' A box to choose the differential expression method (for all studies, with
#' an optional override per study) and the meta-analysis method, with short
#' guidance on when to use each.
#'
#' @param id Module id.
#' @param per_study If `TRUE`, show a DE-method override for each study of
#'   the design passed to `mod_settings_server()`.
#' @param design A reactive returning the experimental design, used to list
#'   the studies for the per-study override (`NULL` if `per_study = FALSE`).
#' @return `mod_settings_ui()` returns UI. `mod_settings_server()` returns a
#'   reactive list with `de_method` (one method, or a named vector with one
#'   method per study) and `meta_method`.
#' @name mod_settings
#' @export
mod_settings_ui <- function(id, per_study = FALSE) {
  ns <- shiny::NS(id)
  shinydashboard::box(
    title = "Analysis settings", status = "info", solidHeader = TRUE, width = 12,
    collapsible = TRUE,
    shiny::fluidRow(
      shiny::column(4, shiny::selectInput(ns("de_method"), "Differential expression method",
                                          choices = de_methods(), selected = "DESeq2")),
      shiny::column(4, shiny::radioButtons(ns("meta_method"), "Meta-analysis method",
                                           choices = meta_methods(), selected = "fisher")),
      if (per_study) shiny::column(4, shiny::uiOutput(ns("per_study")))
    ),
    method_guidance()
  )
}

#' @rdname mod_settings
#' @export
mod_settings_server <- function(id, design = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    studies <- shiny::reactive({
      if (is.null(design)) return(character(0))
      unique(design()$study)
    })
    study_input <- function(study) paste0("study_method_", make.names(study))

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
      list(de_method = de_method, meta_method = meta_method)
    })
  })
}

# Short in-app guidance on choosing the methods
method_guidance <- function() {
  shiny::tags$details(
    shiny::tags$summary(shiny::strong("Which method should I use?")),
    shiny::tags$p(shiny::strong("Differential expression"), " (run in each study):"),
    shiny::tags$ul(
      shiny::tags$li(shiny::strong("DESeq2"), " (default): negative-binomial model with shrunken dispersion estimates. Works well with few replicates per group."),
      shiny::tags$li(shiny::strong("edgeR"), " (quasi-likelihood F-test): also a negative-binomial model; the quasi-likelihood test keeps the false discovery rate under control with few replicates."),
      shiny::tags$li(shiny::strong("limma-voom"), ": linear models on log-counts with precision weights. Fast, and a good choice for larger studies (several replicates per group) or when library sizes differ a lot.")
    ),
    shiny::tags$p("In a benchmark with few replicates, DESeq2 and edgeR were among the best-performing tools (Schurch et al., 2016, RNA 22:839-851). Use the per-study setting if one study differs from the others, e.g. it has many more samples."),
    shiny::tags$p(shiny::strong("Meta-analysis"), " (combines the p-values of the studies):"),
    shiny::tags$ul(
      shiny::tags$li(shiny::strong("Fisher"), ": sums -log(p) over the studies. A very strong signal in a single study can be enough to make a gene significant."),
      shiny::tags$li(shiny::strong("Inverse normal"), ": combines the studies' z-scores, weighting each study by the square root of its number of samples. It favours genes with consistent evidence across studies and gives larger studies more weight (Rau et al., 2014, BMC Bioinformatics 15:91).")
    ),
    shiny::tags$p("On the RNA-Seq tab, click Start Analysis! again after changing a setting; the Demo tab updates by itself.")
  )
}
