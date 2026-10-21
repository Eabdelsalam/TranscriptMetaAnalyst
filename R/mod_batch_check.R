#' Batch-effect check module
#'
#' A box that runs [batch_check()] on the input data before the analysis and
#' shows, for each study, the variance explained by condition and batch with
#' advice, and the PCA and RLE plots.
#'
#' @param id Module id.
#' @param data A reactive returning a list with `counts` and `design` (it may
#'   stop with [shiny::req()] while the inputs are not ready).
#' @return `mod_batch_check_ui()` returns UI; `mod_batch_check_server()`
#'   returns a reactive with the `tma_batch_check`.
#' @name mod_batch_check
#' @examples
#' # the UI part; the server part runs inside a Shiny app
#' ui <- mod_batch_check_ui("rnaseq_batch_check")
#' if (interactive()) {
#'   shiny::shinyApp(shinydashboard::dashboardPage(
#'     shinydashboard::dashboardHeader(), shinydashboard::dashboardSidebar(),
#'     shinydashboard::dashboardBody(ui)), function(input, output, session) NULL)
#' }
#' @export
mod_batch_check_ui <- function(id) {
  ns <- shiny::NS(id)
  shinydashboard::box(
    title = "Check for batch effects (before the analysis)", status = "warning",
    solidHeader = TRUE, width = 12, collapsible = TRUE,
    shiny::p("Looks for batch effects in each study before the differential expression analysis: how much of the variation between samples follows the condition and how much follows the batch (if the design file has a 'batch' column), with PCA and RLE plots. Use it to decide on the batch setting above."),
    shiny::actionButton(ns("run"), "Check for batch effects", icon = shiny::icon("magnifying-glass")),
    shiny::tableOutput(ns("table")),
    # the plots take space only once the check has run
    shiny::uiOutput(ns("plots"))
  )
}

#' @rdname mod_batch_check
#' @export
mod_batch_check_server <- function(id, data) {
  shiny::moduleServer(id, function(input, output, session) {
    check <- shiny::eventReactive(input$run, {
      d <- data()
      batch_check(d$counts, d$design)
    })

    output$table <- shiny::renderTable({
      tab <- check()$table
      data.frame(Study = tab$study,
                 Samples = tab$samples,
                 `Variance explained by condition (%)` = sprintf("%.0f", tab$condition_var),
                 `Variance explained by batch (%)` = ifelse(is.na(tab$batch_var), "no batch column",
                                                           sprintf("%.0f", tab$batch_var)),
                 Advice = tab$message,
                 check.names = FALSE)
    })
    pca <- shiny::reactive(plot_batch_pca(check()))
    rle <- shiny::reactive(plot_rle(check()))
    output$pca <- shiny::renderPlot(pca())
    output$rle <- shiny::renderPlot(rle())
    output$plots <- shiny::renderUI({
      check()
      ns <- session$ns
      shiny::fluidRow(
        shiny::column(6, shiny::plotOutput(ns("pca"), height = "380px"),
                      shiny::downloadButton(ns("download_pca"), label = "Download PCA")),
        shiny::column(6, shiny::plotOutput(ns("rle"), height = "380px"),
                      shiny::downloadButton(ns("download_rle"), label = "Download RLE plot"))
      )
    })
    output$download_pca <- shiny::downloadHandler(paste0(id, "_pca.pdf"), content = function(file) {
      ggplot2::ggsave(file, pca(), device = "pdf", width = 25, height = 12, units = "cm")
    })
    output$download_rle <- shiny::downloadHandler(paste0(id, "_rle.pdf"), content = function(file) {
      ggplot2::ggsave(file, rle(), device = "pdf", width = 25, height = 12, units = "cm")
    })
    check
  })
}
