test_that("the result plots build without errors on the demo data", {
  res <- demo_result()

  expect_s3_class(plot_pca(res), "ggplot")
  expect_s3_class(plot_venn(deg_sets(res)), "ggplot")
  expect_s4_class(plot_upset(deg_sets(res)), "Heatmap")
  vol <- plot_volcano(volcano_data(res))
  expect_s3_class(vol, "ggplot")
  expect_no_error(ggplot2::ggplot_build(vol))

  counts_plot <- plot_gene_counts(res, res$meta_degs[1])
  built <- ggplot2::ggplot_build(counts_plot)
  expect_equal(nrow(built$data[[1]]), nrow(res$design))

  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  expect_no_error(plot_pvalue_hist(res))
  grDevices::dev.off()
})

test_that("plot_gene_counts shows normalised counts + 0.5 like DESeq2::plotCounts", {
  res <- demo_result()
  gene <- res$meta_degs[1]
  p <- plot_gene_counts(res, gene)
  expect_equal(p$data$count, unname(res$norm_counts[gene, res$design$sample]) + 0.5)
  expect_equal(p$data$study, res$design$study)
})

test_that("deg_sets holds each study's DEGs and the meta DEGs", {
  res <- demo_result()
  sets <- deg_sets(res)
  expect_named(sets, c(res$studies, "Meta"))
  expect_equal(sets$Meta, res$meta_degs)
})

test_that("volcano_data has no zero p-values for the plot", {
  res <- demo_result()
  vol <- volcano_data(res)
  expect_false(any(vol$pval == 0, na.rm = TRUE))
  expect_false(any(vol$capped))

  # a p-value that is exactly 0 is capped at the smallest non-zero one
  res$meta$adjpval[1:2] <- 0
  vol <- volcano_data(res)
  smallest <- min(res$meta$adjpval[res$meta$adjpval > 0], na.rm = TRUE)
  expect_equal(vol$pval[1:2], c(smallest, smallest))
  expect_equal(which(vol$capped), 1:2)
})

test_that("the PCA plot has no fixed aspect ratio and labels the variance", {
  res <- demo_result()
  p <- plot_pca(res)
  expect_false(inherits(p$coordinates, "CoordFixed"))
  expect_match(p$labels$x, "^PC1: [0-9]+% variance$")
  expect_equal(nrow(p$data), nrow(res$design))
  expect_setequal(p$data$study, res$studies)
})
