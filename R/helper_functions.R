#' Make a sidebar menu item with sub-items behave as a tab
#'
#' Lets a [shinydashboard::menuItem()] that contains inputs also switch the
#' body to its tab when clicked (from Stack Overflow).
#'
#' @param mi A `menuItem`.
#' @param tabName The tab name the item should open.
#' @return The modified `menuItem`.
#' @export
convertMenuItem <- function(mi,tabName) {
  mi$children[[1]]$attribs['data-toggle']="tab"
  mi$children[[1]]$attribs['data-value'] = tabName
  mi
}

#' Plot the normalised counts of a gene from DESeq2 objects
#'
#' Superseded by [plot_gene_counts()]; kept until the app uses the results
#' module.
#'
#' @param dds_list Named list of `DESeqDataSet` objects, one per study.
#' @param gene_name Gene to plot.
#' @param nrows,ncols Facet layout.
#' @return A ggplot.
#' @keywords internal
#' @export
draw_counts <- function(dds_list, gene_name, nrows, ncols) {
  #extract counts of the selected gene
  counts_list <- lapply(dds_list, function(x) {
    plotCounts(x, gene=gene_name, intgroup="condition", returnData=TRUE)
  })
  #create a dataframe of counts
  counts_df <- dplyr::bind_rows(counts_list, .id = "study")
  #create the plot
  plt <- ggplot(counts_df, aes(x=condition, y=count))
  plt <- plt + geom_point(position=position_jitter(w=0.1,h=0))
  plt <- plt + labs(title = paste("Counts of", gene_name, "across studies"),
                    x = "Condition", y = "Normalized counts")
  plt <- plt + scale_y_log10()
  plt <- plt + facet_wrap(~study, nrow = nrows)
  plt
}