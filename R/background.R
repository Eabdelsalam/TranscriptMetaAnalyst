#' Set up background workers for the app
#'
#' Long analyses run in a background R process (via `future`) so the app stays
#' responsive. This sets a `future` plan if none has been chosen: forked
#' workers (`multicore`) where supported, otherwise separate R sessions
#' (`multisession`). An existing non-sequential plan is left unchanged.
#'
#' @param workers Maximum number of background analyses at the same time.
#' @return The previous `future` plan, invisibly.
#' @export
#' @examples
#' \dontrun{
#' tma_use_background_workers()
#' }
tma_use_background_workers <- function(workers = 2) {
  if (!inherits(future::plan(), "sequential")) {
    return(invisible(future::plan()))
  }
  if (future::supportsMulticore()) {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::multisession, workers = workers)
  }
}
