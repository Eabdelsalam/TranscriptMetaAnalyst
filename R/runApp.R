#' Run the TranscriptMetaAnalyst app
#'
#' Starts the TMA Shiny app locally. The app stops when its browser window is
#' closed.
#'
#' @return Called for its side effect of running the app.
#' @export
#' @examples
#' \dontrun{
#' runTMA()
#' }
runTMA <- function() {
  appDir <- tma_app_dir()
  if (appDir == "") {
    stop("Could not find the TMA app. Try re-installing `TranscriptMetaAnalyst`.", call. = FALSE)
  }
  
  #local run: let the app stop when the browser window is closed
  old_opts <- options(tma.stop_on_session_end = TRUE)
  on.exit(options(old_opts), add = TRUE)
  shiny::runApp(appDir, display.mode = "normal")
}

# Folder of the Shiny app. Looked up from inside the package, so that under
# devtools::load_all() (tests, development) it is the app of the source tree
# and not that of an installed copy
tma_app_dir <- function() {
  system.file("app", package = "TranscriptMetaAnalyst")
}
