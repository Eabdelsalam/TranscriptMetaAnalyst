#' Make a sidebar menu item with sub-items behave as a tab
#'
#' Lets a [shinydashboard::menuItem()] that contains inputs also switch the
#' body to its tab when clicked (from Stack Overflow).
#'
#' @param mi A `menuItem`.
#' @param tabName The tab name the item should open.
#' @return The modified `menuItem`.
#' @export
#' @examples
#' item <- shinydashboard::menuItem("RNA-Seq", tabName = "rnaseq",
#'                                  shiny::fileInput("counts", "Counts file"))
#' convertMenuItem(item, "rnaseq")
convertMenuItem <- function(mi, tabName) {
  mi$children[[1]]$attribs["data-toggle"] <- "tab"
  mi$children[[1]]$attribs["data-value"] <- tabName
  mi
}
