# Regression tests on the bundled demo data (3 studies x 4 samples, Synechocystis)

test_that("demo heatmap matrix has the meta DEGs as rows and every sample as a column", {
  demo <- demo_pipeline()
  heat <- scale_deg_counts(demo$norm_counts, demo$meta_degs, demo$design$sample)

  expect_equal(rownames(heat), demo$meta_degs)
  expect_equal(colnames(heat), demo$design$sample)
  expect_false(anyNA(heat))
})

test_that("demo heatmap values belong to the right gene and sample", {
  demo <- demo_pipeline()
  heat <- scale_deg_counts(demo$norm_counts, demo$meta_degs, demo$design$sample)

  for (gene in utils::head(demo$meta_degs, 25)) {
    for (study in demo$studies) {
      x <- demo$norm_counts[[study]][gene, ]
      expected <- if (stats::sd(x) > 0) (x - mean(x)) / stats::sd(x) else x * 0
      expect_equal(heat[gene, names(x)], expected, tolerance = 1e-10)
    }
  }
})

test_that("demo heatmap renders with heatmaply and pheatmap", {
  skip_if_not_installed("heatmaply")
  skip_if_not_installed("pheatmap")
  demo <- demo_pipeline()
  heat <- scale_deg_counts(demo$norm_counts, demo$meta_degs, demo$design$sample)
  k_row <- heatmap_row_clusters(nrow(heat))
  k_col <- length(demo$studies)

  expect_no_error(suppressWarnings(heatmaply::heatmaply(heat, k_row = k_row, k_col = k_col)))
  # fewer DEGs than samples: the case that used to crash
  expect_no_error(suppressWarnings(heatmaply::heatmaply(
    heat[1:3, ], k_row = heatmap_row_clusters(3), k_col = k_col
  )))

  annot <- demo$design[, c("condition", "study")]
  rownames(annot) <- demo$design$sample
  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  expect_no_error(pheatmap::pheatmap(heat, cutree_rows = k_row, cutree_cols = k_col,
                                     show_rownames = FALSE, annotation_col = annot))
  grDevices::dev.off()
})

test_that("demo PCA matrix keeps the normalised counts of every sample", {
  demo <- demo_pipeline()
  full <- combine_study_matrices(demo$norm_counts, demo$design$sample)

  expect_equal(dim(full), c(length(demo$genenames), nrow(demo$design)))
  expect_equal(rownames(full), demo$genenames)
  for (study in demo$studies) {
    m <- demo$norm_counts[[study]]
    expect_equal(full[, colnames(m)], m)
  }
})

test_that("demo volcano mean fold change matches the per-gene average", {
  demo <- demo_pipeline()
  mfc <- mean_fold_change(demo$fold_changes)
  fc <- do.call(cbind, demo$fold_changes)

  expect_length(mfc, length(demo$genenames))
  idx <- c(1, 4, 100, 1000, length(demo$genenames))
  for (i in idx) {
    expect_equal(mfc[i], mean(fc[i, ], na.rm = TRUE))
  }
})

test_that("demo enrichment genes are the meta-significant, direction-consistent DEGs", {
  demo <- demo_pipeline()
  degs <- demo$meta_degs
  consistent <- rownames(demo$fc_select)[which(abs(demo$fc_select$signFC) == 1)]
  meta_adj <- stats::setNames(demo$fishcomb$adjpval, demo$genenames)

  expect_gt(length(degs), 0)
  expect_true(all(meta_adj[degs] <= 0.05))
  expect_true(all(degs %in% consistent))
  # the old code sent every direction-consistent gene, significant or not
  expect_lt(length(degs), length(consistent))
})

test_that("demo meta-analysis finds the recorded number of DEGs", {
  demo <- demo_pipeline()
  # recorded with DESeq2 1.52.0 / metaRNASeq 1.0.8 on the bundled demo data; a change
  # here means the statistical pipeline changed and should be reviewed
  expect_length(demo$meta_degs, 373)
  expect_equal(sum(abs(demo$fc_select$signFC) == 1, na.rm = TRUE), 908)
})
