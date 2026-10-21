#' Results module
#'
#' The results area shared by the RNA-seq and Demo tabs: information boxes,
#' main plots, the meta-analysis table, functional annotation and gene count
#' plots. The server part renders a `tma_result`; every output is declared once
#' and updates when the result changes.
#'
#' @param id Module id.
#' @param select_organism If `TRUE`, the functional annotation tabs show inputs
#'   to choose the GO organism and the KEGG organism code. If `FALSE`, the
#'   organism given to `mod_results_server()` is used.
#' @param hidden If `TRUE`, the results stay hidden until the first result is
#'   available.
#' @param result A reactive returning a `tma_result` (see [run_tma_pipeline()]).
#' @param organism `NULL` to let the user choose, or a list with `go` (PANTHER
#'   taxon ID) and `kegg` (KEGG organism code).
#' @return `mod_results_ui()` returns UI; `mod_results_server()` is called for
#'   its side effects.
#' @name mod_results
#' @export
mod_results_ui <- function(id, select_organism = FALSE, hidden = FALSE) {
  ns <- shiny::NS(id)
  plot_tab <- function(title, output, ...) {
    shiny::tabPanel(title = title, shiny::plotOutput(ns(output), height = "480px"), ...)
  }
  go_inputs <- if (select_organism) {
    shiny::selectizeInput(ns("organism_go"), "Select organism", choices = "", multiple = FALSE)
  }
  kegg_inputs <- if (select_organism) {
    shiny::textInput(ns("organism_kegg"),
                     label = shiny::p("Enter organism code found at", shiny::tags$a(
                       "KEGG database", target = "_blank",
                       href = "https://www.genome.jp/kegg-bin/get_htext?br08601")))
  }

  ui <- shiny::tagList(
    shiny::fluidRow(
      shinydashboard::infoBoxOutput(ns("num_studies"), width = 3),
      shinydashboard::box(title = "Download DEGs of individual studies", solidHeader = TRUE,
                          width = 4, height = 90, status = "primary",
                          shiny::column(7, shiny::selectInput(ns("individual_study"), label = NULL,
                                                              choices = "", multiple = FALSE)),
                          shiny::column(4, shiny::downloadButton(ns("download_individual"),
                                                                 label = "Download DEGs"))),
      shinydashboard::infoBoxOutput(ns("num_degs"), width = 5)
    ),
    shiny::fluidRow(
      shinydashboard::tabBox(
        title = "Main Plots", id = ns("main_plots"),
        plot_tab("p-Values Dist", "p_dist",
                 shiny::downloadButton(ns("download_pdist"), label = "Download plot")),
        plot_tab("PCA", "pca",
                 shiny::downloadButton(ns("download_pca"), label = "Download plot")),
        shiny::tabPanel(title = "Heatmap",
                        shiny::numericInput(ns("heatmap_k"), label = "Number of gene clusters",
                                            value = 2, min = 1, step = 1, width = "200px"),
                        plotly::plotlyOutput(ns("heatmap"), height = "480px"),
                        shiny::downloadButton(ns("download_heatmap"), label = "Download plot"),
                        shiny::downloadButton(ns("download_heatmap_data"), label = "Download data")),
        plot_tab("Venn diagram", "venn",
                 shiny::downloadButton(ns("download_venn"), label = "Download plot")),
        plot_tab("UpSet Plot", "upset",
                 shiny::downloadButton(ns("download_upset"), label = "Download plot")),
        plot_tab("Volcano Plot", "volcano",
                 shiny::downloadButton(ns("download_volcano"), label = "Download plot"),
                 shiny::downloadButton(ns("download_volcano_data"), label = "Download data")),
        shiny::tabPanel(title = "Batch effects", value = "batch",
                        shiny::tableOutput(ns("batch_status")),
                        shiny::plotOutput(ns("batch_pca"), height = "420px"),
                        shiny::downloadButton(ns("download_batch_pca"), label = "Download PCA"),
                        shiny::plotOutput(ns("batch_rle"), height = "420px"),
                        shiny::downloadButton(ns("download_batch_rle"), label = "Download RLE plot"))
      ),
      shinydashboard::box(title = "Meta-analysis DEGs", status = "success", solidHeader = TRUE,
                          DT::DTOutput(ns("meta_table")),
                          shiny::downloadButton(ns("download_meta_table"), label = "Download table"))
    ),
    shiny::fluidRow(
      shinydashboard::tabBox(
        title = "Functional Annotation",
        shiny::tabPanel(title = "Gene Ontology",
                        go_inputs,
                        shiny::actionButton(ns("run_go"), "Run GO enrichment", icon = shiny::icon("person-running")),
                        shinycssloaders::withSpinner(shiny::plotOutput(ns("go_plot"), height = "450px"), type = 6),
                        shiny::downloadButton(ns("download_go_table"), label = "Download GO table"),
                        shiny::downloadButton(ns("download_go_plot"), label = "Download GO plot")),
        shiny::tabPanel(title = "KEGG pathway",
                        kegg_inputs,
                        shiny::actionButton(ns("run_kegg"), "Run KEGG enrichment", icon = shiny::icon("person-running")),
                        shinycssloaders::withSpinner(shiny::plotOutput(ns("kegg_plot"), height = "450px"), type = 6),
                        shiny::downloadButton(ns("download_kegg_table"), label = "Download KEGG table"),
                        shiny::downloadButton(ns("download_kegg_plot"), label = "Download KEGG plot"))
      ),
      shinydashboard::box(title = "Count plots", status = "primary", solidHeader = TRUE,
                          shiny::textOutput(ns("counts_message")),
                          shiny::plotOutput(ns("counts_plot"), height = "550px"),
                          shiny::downloadButton(ns("download_counts_plot"), label = "Download counts plot"))
    )
  )
  if (hidden) shinyjs::hidden(shiny::div(id = ns("results"), ui)) else ui
}

#' @rdname mod_results
#' @export
mod_results_server <- function(id, result, organism = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    # download file names start with the tab, e.g. "rnaseq_pca.pdf"
    fname <- function(name) paste0(id, "_", name)
    shinycssloaders::hideSpinner("go_plot")
    shinycssloaders::hideSpinner("kegg_plot")

    # update the parts of the page that depend on the studies ----------------
    shiny::observeEvent(result(), {
      res <- result()
      shinyjs::show("results")
      shiny::updateSelectInput(session, "individual_study", choices = res$studies)
      # Venn diagrams are only drawn for up to 3 studies (+ the meta-analysis)
      if (length(res$studies) > 3) {
        shiny::hideTab("main_plots", target = "Venn diagram")
      } else {
        shiny::showTab("main_plots", target = "Venn diagram")
      }
      # the batch diagnostics only exist when the design has a batch column
      if (is.null(res$batch)) {
        shiny::hideTab("main_plots", target = "batch")
      } else {
        shiny::showTab("main_plots", target = "batch")
      }
    })

    # information boxes ------------------------------------------------------
    output$num_studies <- shinydashboard::renderInfoBox({
      res <- result()
      meta <- names(meta_methods())[meta_methods() == res$meta_method]
      shinydashboard::infoBox(title = "Number of studies",
                              paste(length(res$studies), "studies"),
                              paste0(paste(unique(res$de_methods), collapse = ", "), "; ", meta, " meta-analysis"),
                              width = 3, icon = shiny::icon("list-ol"))
    })
    output$num_degs <- shinydashboard::renderInfoBox({
      res <- result()
      diff <- res$idd_irr
      shinydashboard::infoBox(
        title = "Siginficant DEGs (meta-analysis)",
        paste(length(res$meta_degs), "out of", length(res$genes), "genes were identified as DEGs"),
        paste0("The meta-analysis identified ", diff[["IDD"]], "(", diff[["IDR"]], "%) new gene(s), while ",
               diff[["Loss"]], "(", diff[["IRR"]], "%) are sidetracked."),
        width = 4, icon = shiny::icon("chart-simple"))
    })
    output$download_individual <- shiny::downloadHandler(
      filename = function() paste0(input$individual_study, "_DEGs.csv"),
      content = function(file) {
        utils::write.csv(study_significant(result(), input$individual_study), file,
                         quote = FALSE, row.names = FALSE)
      },
      contentType = "text/csv"
    )

    # main plots -------------------------------------------------------------
    output$p_dist <- shiny::renderPlot(plot_pvalue_hist(result()))
    output$download_pdist <- shiny::downloadHandler(fname("pdist.pdf"), content = function(file) {
      grDevices::pdf(file, width = 6, height = 4.5)
      on.exit(grDevices::dev.off())
      plot_pvalue_hist(result())
    })

    pca_plot <- shiny::reactive(plot_pca(result()))
    output$pca <- shiny::renderPlot(pca_plot())
    output$download_pca <- shiny::downloadHandler(fname("pca.pdf"), content = function(file) {
      ggplot2::ggsave(file, pca_plot(), device = "pdf", width = 20, height = 16, units = "cm")
    })

    heat <- shiny::reactive(heatmap_matrix(result()))
    heat_k <- shiny::reactive(heatmap_row_clusters(nrow(heat()), input$heatmap_k))
    need_heat <- function() {
      shiny::validate(shiny::need(nrow(heat()) >= 2,
                                  "At least two meta-analysis DEGs are needed to draw a heatmap."))
    }
    output$heatmap <- plotly::renderPlotly({
      need_heat()
      plot_heatmap(heat(), heat_k(), length(result()$studies))
    })
    output$download_heatmap <- shiny::downloadHandler(fname("heatmap.pdf"), content = function(file) {
      need_heat()
      grDevices::pdf(file, height = 10)
      on.exit(grDevices::dev.off())
      plot_heatmap_static(heat(), result()$design, heat_k(), length(result()$studies))
    })
    output$download_heatmap_data <- shiny::downloadHandler(fname("heatmap_data.csv"), content = function(file) {
      utils::write.csv(heat(), file, quote = FALSE)
    }, contentType = "text/csv")

    sets <- shiny::reactive(deg_sets(result()))
    venn_plot <- shiny::reactive(plot_venn(sets()))
    output$venn <- shiny::renderPlot(venn_plot())
    output$download_venn <- shiny::downloadHandler(fname("venndiagram.pdf"), content = function(file) {
      ggplot2::ggsave(file, venn_plot(), device = "pdf", width = 15, height = 15, units = "cm")
    })

    upset_plot <- shiny::reactive(plot_upset(sets()))
    output$upset <- shiny::renderPlot(ComplexHeatmap::draw(upset_plot()))
    output$download_upset <- shiny::downloadHandler(fname("upset.pdf"), content = function(file) {
      n <- length(result()$studies)
      grDevices::pdf(file, width = max(6, n * 2), height = max(4, n))
      on.exit(grDevices::dev.off())
      ComplexHeatmap::draw(upset_plot())
    })

    volcano <- shiny::reactive(volcano_data(result()))
    volcano_plot <- shiny::reactive(plot_volcano(volcano()))
    output$volcano <- shiny::renderPlot(volcano_plot())
    output$download_volcano <- shiny::downloadHandler(fname("volcano_plot.pdf"), content = function(file) {
      ggplot2::ggsave(file, volcano_plot(), device = "pdf", width = 30, height = 25, units = "cm")
    })
    output$download_volcano_data <- shiny::downloadHandler(fname("volcano_data.csv"), content = function(file) {
      utils::write.csv(volcano(), file, quote = FALSE, row.names = FALSE)
    }, contentType = "text/csv")

    # batch effects ----------------------------------------------------------
    output$batch_status <- shiny::renderTable({
      batch <- result()$batch
      shiny::req(batch)
      method <- names(batch_methods())[batch_methods() == batch$method]
      data.frame(Study = names(batch$status),
                 `Batch handling` = ifelse(batch$status == "applied", method,
                                           paste0("not used (", batch$status, ")")),
                 check.names = FALSE)
    })
    batch_pca <- shiny::reactive({
      shiny::req(result()$batch)
      plot_batch_pca(result())
    })
    batch_rle <- shiny::reactive({
      shiny::req(result()$batch)
      plot_rle(result())
    })
    output$batch_pca <- shiny::renderPlot(batch_pca())
    output$batch_rle <- shiny::renderPlot(batch_rle())
    output$download_batch_pca <- shiny::downloadHandler(fname("batch_pca.pdf"), content = function(file) {
      ggplot2::ggsave(file, batch_pca(), device = "pdf", width = 25, height = 15, units = "cm")
    })
    output$download_batch_rle <- shiny::downloadHandler(fname("batch_rle.pdf"), content = function(file) {
      ggplot2::ggsave(file, batch_rle(), device = "pdf", width = 25, height = 15, units = "cm")
    })

    # meta-analysis table and count plots ------------------------------------
    output$meta_table <- DT::renderDT(meta_table_widget(result()$meta))
    output$download_meta_table <- shiny::downloadHandler(fname("meta_analysis_results.csv"), content = function(file) {
      utils::write.csv(result()$meta, file, quote = FALSE)
    }, contentType = "text/csv")

    selected_gene <- shiny::reactive({
      row <- input$meta_table_rows_selected
      shiny::req(row)
      rownames(result()$meta)[row]
    })
    count_rows <- shiny::reactive(max(1, length(result()$studies) %/% 2))
    count_cols <- shiny::reactive(ceiling(length(result()$studies) / count_rows()))
    counts_plot <- shiny::reactive(plot_gene_counts(result(), selected_gene(), count_rows()))
    output$counts_message <- shiny::renderText({
      if (is.null(input$meta_table_rows_selected)) {
        "Click on any gene in the above DEGs table to obtain a plot of the gene's normalized counts across the included studies"
      } else {
        ""
      }
    })
    output$counts_plot <- shiny::renderPlot(counts_plot())
    output$download_counts_plot <- shiny::downloadHandler(
      filename = function() paste(selected_gene(), fname("counts_plot.pdf"), sep = "_"),
      content = function(file) {
        ggplot2::ggsave(file, counts_plot(), device = "pdf",
                        height = count_rows() * 4, width = count_cols() * 4)
      }
    )

    # functional annotation --------------------------------------------------
    go_organisms <- shiny::reactive(organisms_go())
    if (is.null(organism)) {
      shiny::observeEvent(result(), {
        shiny::updateSelectizeInput(session, "organism_go", choices = go_organisms()$long_name,
                                    server = TRUE)
      }, once = TRUE)
    }
    go_taxon <- shiny::reactive({
      if (!is.null(organism)) return(organism$go)
      shiny::req(input$organism_go)
      go_organisms()$taxon_id[go_organisms()$long_name == input$organism_go]
    })
    kegg_code <- shiny::reactive({
      if (!is.null(organism)) return(organism$kegg)
      shiny::req(input$organism_kegg)
      input$organism_kegg
    })
    have_degs <- function() {
      if (length(result()$meta_degs) == 0) {
        shinyalert::shinyalert("No DEGs", "The meta-analysis found no significant, direction-consistent DEGs to test for enrichment.")
        return(FALSE)
      }
      TRUE
    }

    go_table <- shiny::eventReactive(input$run_go, {
      shiny::req(have_degs())
      shinycssloaders::showSpinner("go_plot")
      on.exit(shinycssloaders::hideSpinner("go_plot"))
      run_go_enrichment(result()$meta_degs, organism = go_taxon())
    })
    go_plot <- shiny::reactive(plot_go_lollipop(go_table()))
    output$go_plot <- shiny::renderPlot(go_plot())
    output$download_go_table <- shiny::downloadHandler(fname("GO_results.txt"), content = function(file) {
      utils::write.table(go_table(), file, quote = FALSE, sep = "\t", row.names = FALSE)
    }, contentType = "text/plain")
    output$download_go_plot <- shiny::downloadHandler(fname("go_plot.pdf"), content = function(file) {
      ggplot2::ggsave(file, go_plot(), device = "pdf", width = 30, height = 20, units = "cm")
    })

    kegg_table <- shiny::eventReactive(input$run_kegg, {
      shiny::req(have_degs())
      shinycssloaders::showSpinner("kegg_plot")
      on.exit(shinycssloaders::hideSpinner("kegg_plot"))
      run_kegg_enrichment(result()$meta_degs, organism = kegg_code())
    })
    kegg_plot <- shiny::reactive(plot_kegg_lollipop(kegg_table()))
    output$kegg_plot <- shiny::renderPlot(kegg_plot())
    output$download_kegg_table <- shiny::downloadHandler(fname("KEGG_results.csv"), content = function(file) {
      utils::write.csv(kegg_table(), file, quote = FALSE)
    }, contentType = "text/csv")
    output$download_kegg_plot <- shiny::downloadHandler(fname("kegg_plot.pdf"), content = function(file) {
      ggplot2::ggsave(file, kegg_plot(), device = "pdf", width = 30, height = 20, units = "cm")
    })
  })
}

# The meta-analysis table as shown in the app: numbers rounded to 4
# significant digits for reading (so very small p-values stay visible); the
# downloaded table keeps full precision
meta_table_widget <- function(meta) {
  numeric_cols <- names(meta)[vapply(meta, is.numeric, logical(1))]
  DT::formatSignif(DT::datatable(meta, selection = "single"), columns = numeric_cols, digits = 4)
}
