test_that("heatmap_genes filters the DEGs for display only", {
  res <- demo_result()
  expect_equal(heatmap_genes(res), res$meta_degs)
  strict <- heatmap_genes(res, padj_cutoff = 1e-5)
  expect_true(all(res$meta[strict, "adjpval"] <= 1e-5))
  big <- heatmap_genes(res, lfc_cutoff = 2)
  expect_true(all(abs(mean_fold_change(res$fold_changes[big, , drop = FALSE])) >= 2))
  top <- heatmap_genes(res, top_n = 20)
  expect_length(top, 20)
  expect_equal(sort(res$meta[top, "adjpval"]), sort(res$meta[res$meta_degs, "adjpval"])[1:20])
  expect_equal(top, res$meta_degs[res$meta_degs %in% top]) # keeps the DEG order
  expect_equal(heatmap_genes(res, top_n = NA), res$meta_degs)
  # the DEGs themselves do not change
  expect_length(res$meta_degs, 373)
})

test_that("the heatmap plots take palette, clustering and gene-name options", {
  res <- demo_result()
  heat <- heatmap_matrix(res, genes = heatmap_genes(res, top_n = 30))
  expect_equal(nrow(heat), 30)
  expect_length(heatmap_palette("blue_white_red", 7), 7)
  expect_error(heatmap_palette("rainbow"))
  for (pal in heatmap_palettes()) {
    expect_no_error(plot_heatmap(heat, 2, 3, palette = pal))
  }
  expect_no_error(plot_heatmap(heat, 2, 3, cluster_rows = FALSE, cluster_cols = FALSE,
                               show_rownames = TRUE))
  expect_no_error(plot_heatmap(heat, 3, 3, cluster_cols = FALSE, method = "ward.D2"))

  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  expect_no_error(plot_heatmap_static(heat, res$design, 2, 3, palette = "viridis",
                                      cluster_rows = FALSE, show_rownames = TRUE))
  expect_no_error(plot_heatmap_static(heat, res$design, 2, 3, method = "average"))
  grDevices::dev.off()
})

test_that("the volcano plot uses the chosen cutoffs and labels", {
  vol <- volcano_data(demo_result())
  p <- plot_volcano(vol, p_cutoff = 0.01, fc_cutoff = 2)
  lines <- ggplot2::ggplot_build(p)$data
  hlines <- unlist(lapply(lines, function(d) d$yintercept))
  vlines <- unlist(lapply(lines, function(d) d$xintercept))
  expect_true(any(abs(hlines - 2) < 1e-8))  # -log10(0.01)
  expect_true(all(c(-2, 2) %in% round(vlines, 8)))

  labelled <- function(plot) {
    d <- ggplot2::ggplot_build(plot)$data
    text <- Filter(function(x) "label" %in% names(x), d)
    if (length(text) == 0) return(character(0))
    labs <- unlist(lapply(text, `[[`, "label"))
    # gene labels only, not the label of the cutoff line
    labs[!is.na(labs) & labs != "" & !grepl("^adj\\. p =", labs)]
  }
  expect_length(labelled(plot_volcano(vol, labels = "none")), 0)
  top <- labelled(plot_volcano(vol, labels = "top", top_n = 5))
  expect_length(top, 5)
  passing <- vol[vol$pval <= 0.05 & abs(vol$fc) >= 1, ]
  expect_setequal(top, passing$genes[order(passing$pval)][1:5])
})

test_that("Venn diagrams are offered for up to 4 sets, otherwise only UpSet", {
  sets <- function(n) stats::setNames(lapply(seq_len(n), function(i) paste0("g", i)), paste0("S", seq_len(n)))
  expect_equal(overlap_plot_type(sets(2)), "venn")
  expect_equal(overlap_plot_type(sets(4)), "venn")
  expect_equal(overlap_plot_type(sets(5)), "upset")
  expect_equal(overlap_plot_type(deg_sets(demo_result())), "venn") # 3 studies + Meta
})

test_that("save_plot writes PDF, PNG and TIFF files of the chosen size", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  signature <- function(file, n) readBin(file, "raw", n)
  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  save_plot(pdf_file, p, "pdf")
  expect_equal(rawToChar(signature(pdf_file, 4)), "%PDF")
  png_file <- withr::local_tempfile(fileext = ".png")
  save_plot(png_file, p, "png", width = 10, height = 5)
  expect_equal(signature(png_file, 4), as.raw(c(0x89, 0x50, 0x4e, 0x47)))
  # 10 x 5 cm at 300 dpi
  dims <- readBin(png_file, "integer", n = 6, size = 4, endian = "big")[5:6]
  expect_equal(dims, round(c(10, 5) / 2.54 * 300), tolerance = 2)
  tiff_file <- withr::local_tempfile(fileext = ".tiff")
  save_plot(tiff_file, function() plot(1:10), "tiff", width = 8, height = 8)
  expect_true(rawToChar(signature(tiff_file, 2)) %in% c("II", "MM"))
  upset_file <- withr::local_tempfile(fileext = ".png")
  expect_no_error(save_plot(upset_file, plot_upset(deg_sets(demo_result())), "png"))
  expect_error(save_plot(pdf_file, p, "svg"))
})

test_that("the results module applies the plot options and download format", {
  res <- demo_result()
  shiny::testServer(mod_results_server,
                    args = list(result = shiny::reactive(res), organism = list(go = 1, kegg = "syn")), {
    expect_equal(heat_genes(), res$meta_degs)
    session$setInputs(heatmap_top = 25, heatmap_lfc = 1)
    expect_length(heat_genes(), 25)
    expect_equal(nrow(heat()), 25)
    session$setInputs(heatmap_top = NA, heatmap_lfc = 0, heatmap_padj = 1e-10)
    expect_true(all(res$meta[heat_genes(), "adjpval"] <= 1e-10))

    session$setInputs(heatmap_palette = "viridis", heatmap_cluster_rows = FALSE, heatmap_names = TRUE)
    expect_equal(heat_options()$palette, "viridis")
    expect_false(heat_options()$cluster_rows)
    expect_true(heat_options()$show_rownames)

    session$setInputs(volcano_p = 0.01, volcano_fc = 2, volcano_labels = "none")
    expect_s3_class(volcano_plot(), "ggplot")

    session$setInputs(download_format = "png", download_width = 12, download_height = 9)
    expect_equal(fname(paste0("pca.", download_format())), paste0(session$ns(NULL), "_pca.png"))
    png_file <- output$download_pca
    expect_equal(readBin(png_file, "raw", 4), as.raw(c(0x89, 0x50, 0x4e, 0x47)))
    heat_png <- output$download_heatmap
    expect_equal(readBin(heat_png, "raw", 4), as.raw(c(0x89, 0x50, 0x4e, 0x47)))
  })
})

test_that("the volcano cutoff line is labelled and the y-axis can be capped", {
  vol <- volcano_data(demo_result())
  built <- ggplot2::ggplot_build(plot_volcano(vol, p_cutoff = 1e-5))
  texts <- unlist(lapply(built$data, function(d) if ("label" %in% names(d)) d$label))
  expect_true("adj. p = 1e-05" %in% texts)

  capped <- plot_volcano(vol, p_cutoff = 1e-5, y_max = 20)
  y <- ggplot2::ggplot_build(capped)$data[[1]]$y
  expect_lte(max(y, na.rm = TRUE), 20 + 1e-8)
  n_above <- sum(-log10(vol$pval) > 20, na.rm = TRUE)
  expect_gt(n_above, 0)
  expect_match(capped$labels$caption, sprintf("^%d gene\\(s\\) above the y-axis maximum", n_above))
  expect_null(plot_volcano(vol)$labels$caption)
})
