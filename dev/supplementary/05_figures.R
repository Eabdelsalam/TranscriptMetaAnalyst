# Manuscript figures 1 and 2 (300 dpi PNG and PDF, 17.4 cm wide).
#
#   fig1_interface: app screenshots, (a) RNA-Seq tab with the Analysis settings
#     and data boxes, (b) Check for batch effects on the simulated batch data
#     (sdlog 0.5, see 01_batch_case_study.R), (c) GEO import of a microarray
#     series, (d) Methods Guide. Needs Chrome or Chromium and an internet
#     connection (GEO).
#   fig2_results: RNA-Seq demo results (DESeq2 + Fisher), (a) meta-analysis
#     p-value histogram, (b) PCA, (c) Venn diagram, (d) UpSet plot, (e) heatmap,
#     (f) volcano plot with the cutoff lines.
# Figure 3 is made by 01_batch_case_study.R.
#
#   NOT_CRAN=true Rscript dev/supplementary/05_figures.R

source("dev/supplementary/00_setup.R")
Sys.setenv(NOT_CRAN = "true")
options(chromote.timeout = 60)

# Figure 2: demo results --------------------------------------------------------
demo <- tma_demo_data()
res <- run_tma_analysis(demo$counts, demo$design, reference = "Cont", de_method = "DESeq2",
                        meta_method = "fisher", batch_method = "none")$results[["Treated vs Cont"]]
# each panel is drawn with save_plot() at its final size and 300 dpi; `scale`
# draws it on a larger canvas so that the text stays readable when printed
panels <- file.path(tempdir(), "fig2")
dir.create(panels, showWarnings = FALSE)
column <- width_double / 2
render_panel <- function(name, plot, height, scale) {
  file <- file.path(panels, paste0(name, ".png"))
  save_plot(file, plot, format = "png", width = column * scale, height = height * scale, dpi = 300 / scale)
  file
}
heat <- heatmap_display(res, tma_plot_options())
files2 <- c(
  render_panel("a", function() plot_pvalue_hist(res), 6.5, 1.6),
  render_panel("b", plot_pca(res) + ggplot2::theme_bw(), 6.5, 1.6),
  render_panel("c", plot_venn(deg_sets(res)), 7, 1.7),
  render_panel("d", plot_upset(deg_sets(res)), 7, 1.7),
  render_panel("e", function() plot_heatmap_static(heat$mat, res$design, heat$k_row, heat$k_col), 11, 1.6),
  render_panel("f", volcano_display(res, tma_plot_options(volcano_labels = "top", volcano_top = 8)), 11, 2.2))
raster_panel <- function(file) {
  patchwork::wrap_elements(full = grid::rasterGrob(png::readPNG(file), interpolate = TRUE))
}
fig2 <- patchwork::wrap_plots(lapply(files2, raster_panel), ncol = 2, heights = c(6.5, 7, 11)) +
  patchwork::plot_annotation(tag_levels = "a") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 9, face = "bold"))
save_figure(fig2, "fig2_results", width = width_double, height = 6.5 + 7 + 11 + 1.5)

# Figure 1: app screenshots ------------------------------------------------------
shots <- file.path(tempdir(), "fig1")
dir.create(shots, showWarnings = FALSE)
batch_files <- local({
  sim <- simulate_studies(seed = seed)
  design <- with_batches(sim$design)
  counts <- add_batch_effect(sim$counts, design, sdlog = 0.5, seed = seed + 1)
  files <- c(counts = file.path(shots, "simulated_counts.txt"), design = file.path(shots, "simulated_design.txt"))
  utils::write.table(counts, files[["counts"]], sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(design, files[["design"]], sep = "\t", quote = FALSE, row.names = FALSE)
  files
})

pkg <- normalizePath(".")
app_fun <- function() {
  pkgload::load_all(PKG, quiet = TRUE)
  options(shiny.maxRequestSize = 100 * 1024^2)
  shiny::shinyAppDir(TranscriptMetaAnalyst:::tma_app_dir())
}
body(app_fun) <- do.call(substitute, list(body(app_fun), list(PKG = pkg)))
app <- shinytest2::AppDriver$new(app_fun, width = 1400, height = 2200,
                                 load_timeout = 120 * 1000, timeout = 600 * 1000)
open_tab <- function(tab) {
  app$run_js(sprintf("document.querySelector('.sidebar-menu a[data-value=\\'%s\\']').click();", tab))
  app$wait_for_idle()
  app$run_js("window.scrollTo(0, 0);")
  Sys.sleep(1)
}
wait_for_text <- function(text, timeout = 600) {
  app$wait_for_js(sprintf("document.body.innerText.includes(%s)", encodeString(text, quote = "\"")),
                  timeout = timeout * 1000)
}
# screenshot of the region covering the boxes (of the visible tab) whose titles
# start with the given texts
shot_boxes <- function(file, titles) {
  rect <- unlist(app$get_js(sprintf(
    "(() => { const titles = %s; let r = null;
       for (const b of document.querySelectorAll('.tab-pane.active .box')) {
         const t = b.querySelector('.box-title');
         if (!t || b.offsetParent === null || !titles.some(x => t.innerText.trim().startsWith(x))) continue;
         const q = b.getBoundingClientRect();
         r = r ? [Math.min(r[0], q.left), Math.min(r[1], q.top), Math.max(r[2], q.right), Math.max(r[3], q.bottom)]
               : [q.left, q.top, q.right, q.bottom];
       }
       window.scrollTo(0, window.scrollY + r[1] - 10);
       return r; })()", jsonlite::toJSON(titles))))
  Sys.sleep(1)
  rect <- unlist(app$get_js(sprintf(
    "(() => { const titles = %s; let r = null;
       for (const b of document.querySelectorAll('.tab-pane.active .box')) {
         const t = b.querySelector('.box-title');
         if (!t || b.offsetParent === null || !titles.some(x => t.innerText.trim().startsWith(x))) continue;
         const q = b.getBoundingClientRect();
         r = r ? [Math.min(r[0], q.left), Math.min(r[1], q.top), Math.max(r[2], q.right), Math.max(r[3], q.bottom)]
               : [q.left, q.top, q.right, q.bottom];
       }
       return r; })()", jsonlite::toJSON(titles))))
  tmp <- tempfile(fileext = ".png")
  app$get_screenshot(tmp, selector = "viewport")
  img <- png::readPNG(tmp)
  rows <- seq(max(1, floor(rect[2]) + 1), min(nrow(img), ceiling(rect[4])))
  cols <- seq(max(1, floor(rect[1]) + 1), min(ncol(img), ceiling(rect[3])))
  png::writePNG(img[rows, cols, , drop = FALSE], file)
}

open_tab("rnaseq")
shot_boxes(file.path(shots, "a.png"), c("Analysis settings", "RNA-Seq data"))

app$upload_file(rnaseq_counts = batch_files[["counts"]])
app$upload_file(rnaseq_design = batch_files[["design"]])
app$click("rnaseq_batch_check-run")
wait_for_text("Variance explained by batch")
Sys.sleep(5)
shot_boxes(file.path(shots, "b.png"), "Check for batch effects")

open_tab("microarray")
app$set_inputs(array_source = "geo")
app$set_inputs(geo_ids = "GSE63352")
app$click("geo_fetch")
wait_for_text("Condition column")
Sys.sleep(3)
shot_boxes(file.path(shots, "c.png"), "Microarray data")

open_tab("methods")
shot_boxes(file.path(shots, "d.png"), c("Scope", "Study requirements"))
app$stop()

# a screenshot panel, with room above it for the panel letter
panel <- function(file) {
  img <- png::readPNG(file)
  pad <- array(1, dim = c(round(0.06 * ncol(img)), ncol(img), dim(img)[3]))
  patchwork::wrap_elements(full = grid::rasterGrob(abind_rows(pad, img), interpolate = TRUE))
}
abind_rows <- function(a, b) {
  out <- array(0, dim = c(dim(a)[1] + dim(b)[1], dim(a)[2], dim(a)[3]))
  out[seq_len(dim(a)[1]), , ] <- a
  out[dim(a)[1] + seq_len(dim(b)[1]), , ] <- b
  out
}
ratio <- function(file) { d <- dim(png::readPNG(file)); (d[1] + round(0.06 * d[2])) / d[2] }
files <- file.path(shots, paste0(c("a", "b", "c", "d"), ".png"))
fig1 <- patchwork::wrap_plots(lapply(files, panel), ncol = 2,
                              heights = c(max(ratio(files[1]), ratio(files[2])),
                                          max(ratio(files[3]), ratio(files[4])))) +
  patchwork::plot_annotation(tag_levels = "a") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 9, face = "bold"))
height <- width_double / 2 * (max(ratio(files[1]), ratio(files[2])) + max(ratio(files[3]), ratio(files[4]))) + 1.5
save_figure(fig1, "fig1_interface", width = width_double, height = height)
file.copy(files, file.path(fig_dir, paste0("fig1_panel_", c("a", "b", "c", "d"), ".png")), overwrite = TRUE)
message("Figures done")
