#' Results module
#'
#' The results area shared by the RNA-seq and Demo tabs: information boxes,
#' main plots, the meta-analysis table, functional annotation and gene count
#' plots. The server part renders a `tma_result`; every output is declared once
#' and updates when the result changes.
#'
#' @param id Module id.
#' @param select_organism If `TRUE`, the functional annotation tabs show inputs
#'   to choose the gene ID type and the PANTHER, KEGG and STRING organisms. If
#'   `FALSE`, the ones given to `mod_results_server()` are used.
#' @param hidden If `TRUE`, the results stay hidden until the first result is
#'   available.
#' @param venn,batch Include the Venn diagram and Batch effects tabs. When the
#'   UI is drawn for a known result, leave out the tabs that do not apply (Venn
#'   for more than 3 studies, Batch effects without a batch column).
#' @param result A reactive returning a `tma_result` (see [run_tma_pipeline()]).
#' @param organism `NULL` to let the user choose, or a list with `go` (PANTHER
#'   taxon ID), `kegg` (KEGG organism code), `string` (STRING taxon ID) and
#'   `id_type` (see [id_types()]).
#' @param file_prefix A reactive giving the start of the download file names;
#'   `NULL` uses the module id (e.g. "demo_pca.pdf").
#' @param toggle_tabs If `TRUE`, the server hides the Venn diagram and Batch
#'   effects tabs when they do not apply to the result. Use `FALSE` when the UI
#'   was drawn for the result with `venn` and `batch` already set: the tabs
#'   then do not exist, and hiding them would fail in the browser.
#' @return `mod_results_ui()` returns UI; `mod_results_server()` is called for
#'   its side effects.
#' @name mod_results
#' @export
mod_results_ui <- function(id, select_organism = FALSE, hidden = FALSE, venn = TRUE, batch = TRUE) {
  ns <- shiny::NS(id)
  plot_tab <- function(title, output, ...) {
    shiny::tabPanel(title = title, shiny::plotOutput(ns(output), height = "480px"), ...)
  }
  go_inputs <- if (select_organism) {
    shiny::selectizeInput(ns("organism_go"), "Organism (PANTHER)", choices = NULL,
                          options = list(placeholder = "Type to search the PANTHER organisms"))
  }
  kegg_inputs <- if (select_organism) {
    shiny::tagList(
      shiny::textInput(ns("organism_kegg"),
                       label = shiny::p("KEGG organism code, from the", shiny::tags$a(
                         "KEGG organism list", target = "_blank",
                         href = "https://www.genome.jp/kegg-bin/get_htext?br08601"))),
      shiny::textOutput(ns("kegg_organism_check"))
    )
  }
  string_inputs <- if (select_organism) {
    shiny::selectizeInput(ns("organism_string"), "Organism (STRING)", choices = NULL,
                          options = list(placeholder = "Type to search the STRING organisms"))
  }
  id_type_input <- if (select_organism) {
    shiny::selectizeInput(ns("id_type"), "Gene ID type of the counts file ('Geneid' column)",
                          choices = c("Choose the ID type" = "", id_types()))
  }
  enrichment_tab <- function(title, key, inputs, run_label, table_label) {
    shiny::tabPanel(title = title,
                    inputs,
                    shiny::actionButton(ns(paste0("run_", key)), run_label, icon = shiny::icon("person-running")),
                    shiny::uiOutput(ns(paste0(key, "_report"))),
                    shinycssloaders::withSpinner(shiny::plotOutput(ns(paste0(key, "_plot")), height = "450px"), type = 6),
                    shiny::downloadButton(ns(paste0("download_", key, "_table")), label = table_label),
                    shiny::downloadButton(ns(paste0("download_", key, "_plot")), label = "Download plot"))
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
      shinydashboard::box(
        title = "Plot download options", status = "primary", solidHeader = FALSE, width = 12,
        collapsible = TRUE, collapsed = TRUE,
        shiny::fluidRow(
          shiny::column(4, shiny::selectInput(ns("download_format"), "File format",
                                              choices = plot_formats(), selected = "pdf")),
          shiny::column(4, shiny::numericInput(ns("download_width"), "Width (cm, empty = default)",
                                               value = NA, min = 2, step = 1)),
          shiny::column(4, shiny::numericInput(ns("download_height"), "Height (cm, empty = default)",
                                               value = NA, min = 2, step = 1))
        ),
        shiny::helpText("Used by every 'Download plot' button in this results tab.")
      )
    ),
    shiny::fluidRow(
      shinydashboard::tabBox(
        title = "Main Plots", id = ns("main_plots"),
        plot_tab("p-Values Dist", "p_dist",
                 shiny::downloadButton(ns("download_pdist"), label = "Download plot")),
        plot_tab("PCA", "pca",
                 shiny::downloadButton(ns("download_pca"), label = "Download plot")),
        shiny::tabPanel(title = "Heatmap",
                        heatmap_options(ns),
                        plotly::plotlyOutput(ns("heatmap"), height = "480px"),
                        shiny::downloadButton(ns("download_heatmap"), label = "Download plot"),
                        shiny::downloadButton(ns("download_heatmap_data"), label = "Download data")),
        if (venn) plot_tab("Venn diagram", "venn",
                           shiny::downloadButton(ns("download_venn"), label = "Download plot")),
        plot_tab("UpSet Plot", "upset",
                 shiny::downloadButton(ns("download_upset"), label = "Download plot")),
        shiny::tabPanel(title = "Volcano Plot",
                        volcano_options(ns),
                        shiny::plotOutput(ns("volcano"), height = "480px"),
                        shiny::downloadButton(ns("download_volcano"), label = "Download plot"),
                        shiny::downloadButton(ns("download_volcano_data"), label = "Download data")),
        if (batch) shiny::tabPanel(title = "Batch effects", value = "batch",
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
      shiny::column(
        6,
        shinydashboard::tabBox(
          title = "Functional Annotation", width = NULL,
          enrichment_tab("Gene Ontology", "go", go_inputs, "Run GO enrichment (PANTHER)", "Download GO table"),
          enrichment_tab("KEGG pathway", "kegg", kegg_inputs, "Run KEGG enrichment", "Download KEGG table"),
          enrichment_tab("STRING", "string", string_inputs, "Run STRING enrichment", "Download STRING table")
        ),
        shinydashboard::box(
          title = "Enrichment input", width = NULL, status = "primary", solidHeader = FALSE,
          id_type_input,
          shiny::textOutput(ns("enrichment_input")),
          shiny::helpText("The meta-analysis DEGs are tested against all genes tested in the analysis (the background). PANTHER and STRING recognise most ID types; KEGG accepts locus tags / KEGG gene IDs, NCBI Entrez Gene IDs and UniProt accessions.")
        )
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
mod_results_server <- function(id, result, organism = NULL, file_prefix = NULL, toggle_tabs = TRUE) {
  shiny::moduleServer(id, function(input, output, session) {
    # download file names start with the tab, e.g. "demo_pca.pdf"
    fname <- function(name) {
      prefix <- if (is.null(file_prefix)) id else file_prefix()
      paste0(prefix, "_", name)
    }
    shinycssloaders::hideSpinner("go_plot")
    shinycssloaders::hideSpinner("kegg_plot")

    # update the parts of the page that depend on the studies ----------------
    shiny::observeEvent(result(), {
      res <- result()
      # sent after the page is updated, so it also works when the results UI is
      # drawn at the same time (one tab per contrast)
      # the callback runs outside the module's reactive domain, so set it again
      session$onFlushed(function() shiny::withReactiveDomain(session, {
        shinyjs::show("results")
        shiny::updateSelectInput(session, "individual_study", choices = res$studies)
        if (toggle_tabs) {
          # Venn diagrams only for up to 4 sets (studies + meta-analysis)
          if (overlap_plot_type(deg_sets(res)) == "upset") {
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
        }
      }), once = TRUE)
    })

    # information boxes ------------------------------------------------------
    output$num_studies <- shinydashboard::renderInfoBox({
      res <- result()
      contrast <- if (is.null(res$contrast)) "" else paste0(res$contrast, "; ")
      meta <- names(meta_methods())[meta_methods() == res$meta_method]
      shinydashboard::infoBox(title = "Number of studies",
                              paste(length(res$studies), "studies"),
                              paste0(contrast, paste(unique(res$de_methods), collapse = ", "), "; ", meta, " meta-analysis"),
                              width = 3, icon = shiny::icon("list-ol"))
    })
    output$num_degs <- shinydashboard::renderInfoBox({
      res <- result()
      diff <- res$idd_irr
      shinydashboard::infoBox(
        title = "Significant DEGs (meta-analysis)",
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
    # every plot download uses the format and size of "Plot download options";
    # the default size (cm, a number or a function) is the plot's usual size
    download_format <- function() {
      if (is.null(input$download_format)) "pdf" else input$download_format
    }
    plot_download <- function(name, plot, width, height) {
      shiny::downloadHandler(
        filename = function() fname(paste0(name, ".", download_format())),
        content = function(file) {
          size <- function(x) if (is.function(x)) x() else x
          save_plot(file, plot(), format = download_format(),
                    width = input$download_width, height = input$download_height,
                    default_width = size(width), default_height = size(height))
        }
      )
    }

    output$p_dist <- shiny::renderPlot(plot_pvalue_hist(result()))
    output$download_pdist <- plot_download("pdist", function() function() plot_pvalue_hist(result()),
                                           15.24, 11.43)

    pca_plot <- shiny::reactive(plot_pca(result()))
    output$pca <- shiny::renderPlot(pca_plot())
    output$download_pca <- plot_download("pca", pca_plot, 20, 16)

    # the heatmap options only change what is shown, not the DEGs
    heat_genes <- shiny::reactive({
      heatmap_genes(result(), padj_cutoff = number_or(input$heatmap_padj, result()$alpha),
                    lfc_cutoff = number_or(input$heatmap_lfc, 0),
                    top_n = number_or(input$heatmap_top, NA))
    })
    heat <- shiny::reactive(heatmap_matrix(result(), genes = heat_genes()))
    heat_k <- shiny::reactive(heatmap_row_clusters(nrow(heat()), input$heatmap_k))
    heat_options <- shiny::reactive(list(
      palette = if (is.null(input$heatmap_palette)) "yellow_black_purple" else input$heatmap_palette,
      cluster_rows = !isFALSE(input$heatmap_cluster_rows),
      cluster_cols = !isFALSE(input$heatmap_cluster_cols),
      method = if (is.null(input$heatmap_method)) "complete" else input$heatmap_method,
      show_rownames = isTRUE(input$heatmap_names)
    ))
    need_heat <- function() {
      shiny::validate(shiny::need(nrow(heat()) >= 2,
                                  "At least two DEGs are needed to draw a heatmap (check the gene filters above)."))
    }
    output$heatmap <- plotly::renderPlotly({
      need_heat()
      o <- heat_options()
      plot_heatmap(heat(), heat_k(), length(result()$studies), palette = o$palette,
                   cluster_rows = o$cluster_rows, cluster_cols = o$cluster_cols,
                   method = o$method, show_rownames = o$show_rownames)
    })
    output$download_heatmap <- plot_download("heatmap", function() {
      need_heat()
      o <- heat_options()
      mat <- heat()
      design <- result()$design
      k_row <- heat_k()
      k_col <- length(result()$studies)
      function() {
        plot_heatmap_static(mat, design, k_row, k_col, palette = o$palette,
                            cluster_rows = o$cluster_rows, cluster_cols = o$cluster_cols,
                            method = o$method, show_rownames = o$show_rownames)
      }
    }, 17.78, 25.4)
    output$download_heatmap_data <- shiny::downloadHandler(function() fname("heatmap_data.csv"), content = function(file) {
      utils::write.csv(heat(), file, quote = FALSE)
    }, contentType = "text/csv")

    sets <- shiny::reactive(deg_sets(result()))
    venn_plot <- shiny::reactive(plot_venn(sets()))
    output$venn <- shiny::renderPlot(venn_plot())
    output$download_venn <- plot_download("venndiagram", venn_plot, 15, 15)

    upset_plot <- shiny::reactive(plot_upset(sets()))
    output$upset <- shiny::renderPlot(ComplexHeatmap::draw(upset_plot()))
    output$download_upset <- plot_download("upset", upset_plot,
                                           function() max(6, length(result()$studies) * 2) * 2.54,
                                           function() max(4, length(result()$studies)) * 2.54)

    volcano <- shiny::reactive(volcano_data(result()))
    volcano_plot <- shiny::reactive({
      plot_volcano(volcano(),
                   p_cutoff = number_or(input$volcano_p, 0.05),
                   fc_cutoff = number_or(input$volcano_fc, 1),
                   labels = if (is.null(input$volcano_labels)) "significant" else input$volcano_labels,
                   top_n = number_or(input$volcano_top, 10),
                   y_max = number_or(input$volcano_ymax, NA))
    })
    output$volcano <- shiny::renderPlot(volcano_plot())
    output$download_volcano <- plot_download("volcano_plot", volcano_plot, 30, 25)
    output$download_volcano_data <- shiny::downloadHandler(function() fname("volcano_data.csv"), content = function(file) {
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
    output$download_batch_pca <- plot_download("batch_pca", batch_pca, 25, 15)
    output$download_batch_rle <- plot_download("batch_rle", batch_rle, 25, 15)

    # meta-analysis table and count plots ------------------------------------
    output$meta_table <- DT::renderDT(meta_table_widget(result()$meta))
    output$download_meta_table <- shiny::downloadHandler(function() fname("meta_analysis_results.csv"), content = function(file) {
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
      filename = function() paste(selected_gene(), fname(paste0("counts_plot.", download_format())), sep = "_"),
      content = function(file) {
        save_plot(file, counts_plot(), format = download_format(),
                  width = input$download_width, height = input$download_height,
                  default_width = count_cols() * 4 * 2.54, default_height = count_rows() * 4 * 2.54)
      }
    )

    # functional annotation --------------------------------------------------
    # organism lists, fetched once and shared by every tab and session
    go_organisms <- shiny::reactive(organisms_go()) |>
      shiny::bindCache("panther-organisms", cache = "app")
    string_organisms <- shiny::reactive(organisms_string()) |>
      shiny::bindCache("string-organisms", cache = "app")
    if (is.null(organism)) {
      shiny::observeEvent(result(), {
        lists <- list(
          organism_go = function() {
            o <- go_organisms()
            stats::setNames(o$taxon_id, paste0(o$long_name, " (", o$taxon_id, ")"))
          },
          organism_string = function() {
            o <- string_organisms()
            stats::setNames(o$taxon_id, paste0(o$name, " (", o$taxon_id, ")"))
          }
        )
        for (input_id in names(lists)) {
          choices <- tryCatch(lists[[input_id]](), error = function(e) NULL)
          local({
            id <- input_id
            ch <- choices
            session$onFlushed(function() shiny::withReactiveDomain(session, {
              shiny::updateSelectizeInput(session, id, choices = ch, selected = character(0),
                                          server = TRUE)
            }), once = TRUE)
          })
        }
      }, once = TRUE)
    }
    chosen <- function(fixed, value, what) {
      if (!is.null(organism)) return(organism[[fixed]])
      if (is.null(value) || !nzchar(value)) {
        shinyalert::shinyalert("Missing choice", paste("Please choose the", what, "first."))
        shiny::req(FALSE)
      }
      value
    }
    id_type <- function() chosen("id_type", input$id_type, "gene ID type of the counts file")
    output$kegg_organism_check <- shiny::renderText({
      code <- trimws(if (is.null(input$organism_kegg)) "" else input$organism_kegg)
      if (!nzchar(code)) return("")
      name <- kegg_organism_name(code)
      if (is.na(name)) "Unknown KEGG organism code." else paste("Organism:", name)
    })
    output$enrichment_input <- shiny::renderText({
      res <- result()
      sprintf("Input: %d meta-analysis DEGs; background: all %d tested genes.",
              length(res$meta_degs), length(res$genes))
    })
    have_degs <- function() {
      if (length(result()$meta_degs) == 0) {
        shinyalert::shinyalert("No DEGs", "The meta-analysis found no significant, direction-consistent DEGs to test for enrichment.")
        return(FALSE)
      }
      TRUE
    }

    # one enrichment service: run on click, show the mapping report, plot and
    # downloads; failures are shown as a message instead of an error
    enrichment <- function(key, service, run, plot_fun) {
      res <- shiny::eventReactive(input[[paste0("run_", key)]], {
        shiny::req(have_degs())
        args <- run()
        shinycssloaders::showSpinner(paste0(key, "_plot"))
        on.exit(shinycssloaders::hideSpinner(paste0(key, "_plot")))
        out <- tryCatch(do.call(args$fun, args$args), error = function(e) {
          shinyalert::shinyalert(paste(service, "enrichment"), enrichment_message(e, service), type = "error")
          NULL
        })
        shiny::req(out)
        out
      })
      plot <- shiny::reactive({
        shiny::validate(shiny::need(nrow(res()$table) > 0,
                                    sprintf("No enriched terms were found with %s.", service)))
        plot_fun(res()$table)
      })
      output[[paste0(key, "_plot")]] <- shiny::renderPlot(plot())
      output[[paste0(key, "_report")]] <- shiny::renderUI(mapping_report_ui(res(), session$ns(paste0("download_", key, "_unmapped"))))
      output[[paste0("download_", key, "_table")]] <- shiny::downloadHandler(
        function() fname(paste0(key, "_enrichment.txt")),
        content = function(file) write_enrichment(res(), file),
        contentType = "text/plain"
      )
      output[[paste0("download_", key, "_unmapped")]] <- shiny::downloadHandler(
        function() fname(paste0(key, "_unmapped_ids.txt")),
        content = function(file) writeLines(res()$mapping$unmapped_ids, file),
        contentType = "text/plain"
      )
      output[[paste0("download_", key, "_plot")]] <- plot_download(paste0(key, "_plot"), plot, 30, 20)
      res
    }

    go_result <- enrichment("go", "PANTHER", function() {
      list(fun = run_go_enrichment,
           args = list(genes = result()$meta_degs,
                       organism = chosen("go", input$organism_go, "PANTHER organism"),
                       background = result()$genes, id_type = id_type()))
    }, plot_go_lollipop)
    kegg_result <- enrichment("kegg", "KEGG", function() {
      list(fun = run_kegg_enrichment,
           args = list(genes = result()$meta_degs,
                       organism = trimws(chosen("kegg", input$organism_kegg, "KEGG organism code")),
                       background = result()$genes, id_type = id_type()))
    }, plot_kegg_lollipop)
    string_result <- enrichment("string", "STRING", function() {
      list(fun = run_string_enrichment,
           args = list(genes = result()$meta_degs,
                       species = chosen("string", input$organism_string, "STRING organism"),
                       background = result()$genes, id_type = id_type()))
    }, plot_string_lollipop)
  })
}

# The meta-analysis table as shown in the app: numbers rounded to 4
# significant digits for reading (so very small p-values stay visible); the
# downloaded table keeps full precision
meta_table_widget <- function(meta) {
  numeric_cols <- names(meta)[vapply(meta, is.numeric, logical(1))]
  DT::formatSignif(DT::datatable(meta, selection = "single"), columns = numeric_cols, digits = 4)
}

# A numeric input, or `default` while it is empty or invalid
number_or <- function(x, default) {
  if (length(x) == 1 && is.numeric(x) && !is.na(x)) x else default
}

# A clickable toggle (styled as a link) that shows or hides plot options
options_toggle <- function(label, ...) {
  shiny::tags$details(
    shiny::tags$summary(
      style = "color: #337ab7; text-decoration: underline; cursor: pointer; margin-bottom: 8px;",
      label
    ),
    ...
  )
}

# Display options of the heatmap (they do not change the DEGs)
heatmap_options <- function(ns) {
  half <- function(...) shiny::column(6, ...)
  shiny::tagList(
    shiny::fluidRow(
      half(shiny::numericInput(ns("heatmap_k"), "Number of gene clusters", value = 2, min = 1, step = 1)),
      half(shiny::numericInput(ns("heatmap_top"), "Top N genes (empty = all)", value = NA, min = 2, step = 10))
    ),
    options_toggle(
      "More heatmap options",
      shiny::fluidRow(
        half(shiny::numericInput(ns("heatmap_padj"), "Show DEGs with adj. p <=",
                                 value = 0.05, min = 0, max = 1, step = 0.01)),
        half(shiny::numericInput(ns("heatmap_lfc"), "and |mean log2FC| >=", value = 0, min = 0, step = 0.5))
      ),
      shiny::fluidRow(
        half(shiny::selectInput(ns("heatmap_palette"), "Colours",
                                choices = heatmap_palettes(), selected = "yellow_black_purple")),
        half(shiny::selectInput(ns("heatmap_method"), "Clustering method",
                                choices = c("complete", "average", "ward.D2", "single"),
                                selected = "complete"))
      ),
      shiny::fluidRow(
        half(shiny::checkboxInput(ns("heatmap_cluster_rows"), "Cluster genes", value = TRUE),
             shiny::checkboxInput(ns("heatmap_cluster_cols"), "Cluster samples", value = TRUE)),
        half(shiny::checkboxInput(ns("heatmap_names"), "Show gene names", value = FALSE))
      ),
      shiny::helpText("These options only change the heatmap; the DEGs are those of the analysis.")
    )
  )
}

# Display options of the volcano plot
volcano_options <- function(ns) {
  half <- function(...) shiny::column(6, ...)
  options_toggle(
    "Volcano plot options",
    shiny::fluidRow(
      half(shiny::numericInput(ns("volcano_p"), "Adj. p cutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
      half(shiny::numericInput(ns("volcano_fc"), "|log2FC| cutoff", value = 1, min = 0, step = 0.5))
    ),
    shiny::fluidRow(
      half(shiny::selectInput(ns("volcano_labels"), "Gene labels",
                              choices = c("All past both cutoffs" = "significant",
                                          "Top N most significant" = "top",
                                          "None" = "none"),
                              selected = "significant")),
      half(shiny::numericInput(ns("volcano_top"), "N (for top N)", value = 10, min = 1, step = 5))
    ),
    shiny::fluidRow(
      half(shiny::numericInput(ns("volcano_ymax"), "Y-axis maximum (-log10 adj. p, empty = all)",
                               value = NA, min = 1, step = 5)),
      half(shiny::helpText("With very small p-values the axis is long and the cutoff line sits near the bottom; a lower maximum spreads that region out."))
    )
  )
}

# The ID-mapping report and the label shown above an enrichment result
mapping_report_ui <- function(enrichment, download_id) {
  m <- enrichment$mapping
  shiny::div(
    style = "margin: 10px 0;",
    shiny::tags$table(
      class = "table table-condensed", style = "width: auto; margin-bottom: 5px;",
      shiny::tags$tr(shiny::tags$th("IDs submitted"), shiny::tags$th("Mapped"),
                     shiny::tags$th("Not mapped"), shiny::tags$th("Mapping rate")),
      shiny::tags$tr(shiny::tags$td(m$submitted), shiny::tags$td(m$mapped),
                     shiny::tags$td(m$unmapped), shiny::tags$td(sprintf("%.0f%%", 100 * m$rate)))
    ),
    if (!is.null(m$warning)) shiny::div(class = "alert alert-warning", style = "padding: 6px 10px;", m$warning),
    if (m$unmapped > 0) shiny::downloadLink(download_id, "Download the IDs that were not mapped"),
    shiny::p(shiny::em(enrichment_label(enrichment)), style = "margin-top: 6px;")
  )
}

# The enrichment table, preceded by the service, version and background
write_enrichment <- function(enrichment, file) {
  header <- c(paste("#", enrichment_label(enrichment)),
              sprintf("# IDs submitted: %d; mapped: %d (%.0f%%)", enrichment$mapping$submitted,
                      enrichment$mapping$mapped, 100 * enrichment$mapping$rate))
  writeLines(header, file)
  suppressWarnings(utils::write.table(enrichment$table, file, sep = "\t", quote = FALSE,
                                      row.names = FALSE, append = TRUE))
}
