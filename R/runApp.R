runTMA <- function() {
  appDir <- system.file("app", package = "TranscriptMetaAnalyst")
  if (appDir == "") {
    stop("Could not find the TMA app. Try re-installing `TranscriptMetaAnalyst`.", call. = FALSE)
  }
  
  #local run: let the app stop when the browser window is closed
  old_opts <- options(tma.stop_on_session_end = TRUE)
  on.exit(options(old_opts), add = TRUE)
  shiny::runApp(appDir, display.mode = "normal")
}
