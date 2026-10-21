# End-to-end smoke tests: the app runs in a headless browser (shinytest2) and
# every tab is opened. Skipped when shinytest2 or a Chrome/Chromium browser is
# not available, and on CRAN.

# a snap-packaged Chromium can take longer than chromote's default 10 s to start
withr::local_options(chromote.timeout = 60)
# Chrome writes a profile folder to TMPDIR; keep it inside R's session temp
# folder, which is removed when R exits
withr::local_envvar(TMPDIR = tempdir())

skip_if_no_browser <- function() {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")
  skip_if(is.null(tryCatch(suppressMessages(chromote::find_chrome()), error = function(e) NULL)),
          "no Chrome or Chromium browser found")
}

# The app, with the online organism lists replaced by offline stubs. Runs in
# the background R process of shinytest2, which loads the package source.
smoke_app <- function() {
  library(TranscriptMetaAnalyst)
  utils::assignInNamespace("organisms_go", function() {
    data.frame(long_name = "Synechocystis", taxon_id = 1111708)
  }, "TranscriptMetaAnalyst")
  utils::assignInNamespace("organisms_string", function() {
    data.frame(taxon_id = 1148, name = "Synechocystis sp. PCC 6803", domain = "Bacteria")
  }, "TranscriptMetaAnalyst")
  future::plan(future::sequential)
  shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
}

start_app <- function(name) {
  shinytest2::AppDriver$new(smoke_app, name = name, height = 1000, width = 1400,
                            load_timeout = 60 * 1000, timeout = 120 * 1000)
}

# open a tab by clicking its sidebar link, as a user does
open_tab <- function(app, tab) {
  app$run_js(sprintf("document.querySelector('.sidebar-menu a[data-value=\\'%s\\']').click();", tab))
  app$wait_for_idle()
}

page_text <- function(app) app$get_js("document.body.innerText")

# wait until the page shows `text`
wait_for_text <- function(app, text, timeout = 300) {
  app$wait_for_js(sprintf("document.body.innerText.includes(%s)", encodeString(text, quote = "\"")),
                  timeout = timeout * 1000)
}

expect_no_app_errors <- function(app) {
  logs <- as.data.frame(app$get_logs())
  errors <- logs$message[logs$level %in% c("error", "severe") |
                           grepl("^(Error|Warning: Error)", logs$message)]
  expect_identical(errors, character(0))
}

test_that("Home, Methods Guide and User Manual tabs open", {
  skip_if_no_browser()
  app <- start_app("static-tabs")
  withr::defer(app$stop())
  expect_match(page_text(app), "Quick overview")
  open_tab(app, "methods")
  wait_for_text(app, "Study requirements", timeout = 30)
  expect_match(page_text(app), "Conesa")
  open_tab(app, "manual")
  wait_for_text(app, "Contact us", timeout = 30)
  expect_no_app_errors(app)
})

test_that("the RNA-Seq tab runs the demo data", {
  skip_if_no_browser()
  app <- start_app("rnaseq-demo")
  withr::defer(app$stop())
  open_tab(app, "rnaseq")
  app$click("rnaseq_demo")
  wait_for_text(app, "genes were identified as DEGs")
  expect_match(page_text(app), "Meta-analysis DEGs")
  expect_match(page_text(app), "Treated vs Cont")
  expect_no_app_errors(app)
})

test_that("the RNA-Seq tab analyses uploaded files", {
  skip_if_no_browser()
  app <- start_app("rnaseq-tab")
  withr::defer(app$stop())
  open_tab(app, "rnaseq")
  app$upload_file(rnaseq_counts = tma_example_file("counts"))
  app$upload_file(rnaseq_design = tma_example_file("design"))
  app$click("rnaseq_analyze")
  wait_for_text(app, "Download R script")
  expect_match(page_text(app), "Comparisons")
  wait_for_text(app, "genes were identified as DEGs", timeout = 60)
  expect_no_app_errors(app)
})

test_that("the Microarray tab runs the microarray demo", {
  skip_if_no_browser()
  app <- start_app("microarray-tab")
  withr::defer(app$stop())
  open_tab(app, "microarray")
  app$click("array_demo")
  wait_for_text(app, "Data preparation")
  wait_for_text(app, "genes were identified as DEGs")
  expect_match(page_text(app), "measured in all 3 studies")
  expect_match(page_text(app), "LC_24h vs HC")
  expect_no_app_errors(app)
})
