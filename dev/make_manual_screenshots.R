# Takes the app screenshots of manual/TMA_Manual.md with shinytest2 (needs
# Chrome or Chromium, and an internet connection for the organism lists of the
# enrichment tabs). Run from the package root:
#
#   NOT_CRAN=true Rscript -e 'devtools::load_all(); source("dev/make_manual_screenshots.R")'
#
# The result figures drawn with the plot functions are made by
# dev/make_manual_figures.R.

options(chromote.timeout = 60)
fig_dir <- normalizePath(file.path("manual", "figures"))
width <- 1600
height <- 1000

app <- shinytest2::AppDriver$new(
  function() {
    library(TranscriptMetaAnalyst)
    shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  },
  width = width, height = height, load_timeout = 120 * 1000, timeout = 300 * 1000
)

# path of a figure; an existing file is removed first, as screenshots do not
# overwrite files
fig <- function(name) {
  file <- file.path(fig_dir, paste0(name, ".png"))
  unlink(file)
  file
}

# wait until the page shows `text`
wait_for_text <- function(text, timeout = 600) {
  app$wait_for_js(sprintf("document.body.innerText.includes(%s)", encodeString(text, quote = "\"")),
                  timeout = timeout * 1000)
}

# give the box (in the visible tab) whose title starts with `title` an id, so
# it can be captured
mark_box <- function(title, id) {
  app$run_js(sprintf(
    "var boxes = document.querySelectorAll('.tab-pane.active .box, .tab-pane.active .nav-tabs-custom');
     for (var b of boxes) {
       var t = b.querySelector('.box-title, .nav-tabs .header');
       if (t && t.innerText.trim().startsWith(%s) && b.offsetParent !== null) { b.id = %s; break; }
     }", encodeString(title, quote = "'"), encodeString(id, quote = "'")))
}

# viewport with the element `selector` at its top
shot_from <- function(file, selector) {
  app$run_js(sprintf("document.querySelector(%s).scrollIntoView();", encodeString(selector, quote = "'")))
  Sys.sleep(1)
  app$get_screenshot(file, selector = "viewport")
}

# the element `selector`: scrolled into view and cut out of a viewport
# screenshot (screenshots of elements below the first screen were misplaced)
shot_of <- function(file, selector) {
  sel <- encodeString(selector, quote = "'")
  app$run_js(sprintf("document.querySelector(%s).scrollIntoView();", sel))
  Sys.sleep(1)
  rect <- unlist(app$get_js(sprintf(
    "(() => { const r = document.querySelector(%s).getBoundingClientRect();
              return [r.left, r.top, r.right, r.bottom]; })()", sel)))
  tmp <- tempfile(fileext = ".png")
  app$get_screenshot(tmp, selector = "viewport")
  img <- png::readPNG(tmp)
  rows <- seq(max(1, floor(rect[2]) + 1), min(nrow(img), ceiling(rect[4])))
  cols <- seq(max(1, floor(rect[1]) + 1), min(ncol(img), ceiling(rect[3])))
  png::writePNG(img[rows, cols, , drop = FALSE], file)
}

# close the "Some studies were left out" alert, which follows the "Running..."
# notice when the analysis has finished
close_alert <- function() {
  app$wait_for_js("(() => { const a = document.querySelector('.sweet-alert.visible h2');
                           return a !== null && a.innerText.includes('Some studies were left out'); })()",
                  timeout = 60 * 1000)
  Sys.sleep(1)
  app$run_js("swalService.close();")
  # the open alert also blocks scrolling of the page
  app$wait_for_js("document.querySelector('.sweet-alert.visible') === null &&
                   !document.body.classList.contains('stop-scrolling')", timeout = 30 * 1000)
  Sys.sleep(1)
}

# open a tab by clicking its sidebar link (setting the input does not open the
# RNA-Seq item)
open_tab <- function(tab) {
  app$run_js(sprintf("document.querySelector('.sidebar-menu a[data-value=\\'%s\\']').click();", tab))
  app$wait_for_idle()
  app$run_js("window.scrollTo(0, 0);")
  Sys.sleep(1)
}

# Home (Figure 1)
app$get_screenshot(fig("fig01_home_tab"), selector = "viewport")

# RNA-Seq tab before an upload: settings, batch check and quick start (Figure 2)
open_tab("rnaseq")
app$get_screenshot(fig("fig02_rnaseq_tab"), selector = "viewport")

# Microarray tab with the demo results (Figures 5 and 6)
open_tab("microarray")
app$click("array_demo")
wait_for_text("Data preparation")
wait_for_text("genes were identified as DEGs")
close_alert()
mark_box("Microarray data", "fig-array-input")
shot_of(fig("fig05_microarray_input"), "#fig-array-input")
mark_box("Data preparation", "fig-array-preparation")
shot_from(fig("fig06_microarray_preparation"), "#fig-array-preparation")

# Demo tab with the results (Figures 7 and 11)
open_tab("demo")
wait_for_text("genes were identified as DEGs")
Sys.sleep(3)
app$get_screenshot(fig("fig07_demo_tab"), selector = "viewport")
mark_box("Meta-analysis DEGs", "fig-degs")
shot_from(fig("fig11_outputs_overview"), "#shiny-tab-demo .info-box")

# DEGs table (Figure 17a)
shot_of(fig("fig17a_degs_table"), "#fig-degs")

# Functional annotation inputs of an RNA-Seq analysis of the example files:
# GO and KEGG tabs with the organism and ID-type inputs (Figure 10)
open_tab("rnaseq")
app$upload_file(rnaseq_counts = tma_example_file("counts"))
app$upload_file(rnaseq_design = tma_example_file("design"))
app$click("rnaseq_analyze")
wait_for_text("Download R script")
wait_for_text("genes were identified as DEGs")
Sys.sleep(5)
mark_box("Functional Annotation", "fig-enrichment")
# the tab box and the Enrichment input box below it
app$run_js("document.querySelector('#fig-enrichment').closest('.col-sm-6').id = 'fig-enrichment-column';")
shot_of(fig("fig10a_go_input"), "#fig-enrichment-column")
app$run_js("for (const a of document.querySelectorAll('#fig-enrichment .nav-tabs a')) if (a.innerText.trim() === 'KEGG pathway') a.click();")
Sys.sleep(1)
shot_of(fig("fig10b_kegg_input"), "#fig-enrichment-column")

# Methods Guide (Figure 8) and User Manual (Figure 9) tabs
open_tab("methods")
app$get_screenshot(fig("fig08_methods_guide"), selector = "viewport")
open_tab("manual")
app$get_screenshot(fig("fig09_manual_tab"), selector = "viewport")

app$stop()
message("Screenshots written to ", fig_dir)
