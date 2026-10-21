# Regression tests on the bundled demo data (3 studies x 4 samples, Synechocystis)

test_that("demo heatmap matrix has the meta DEGs as rows and every sample as a column", {
  res <- demo_result()
  heat <- heatmap_matrix(res)

  expect_equal(rownames(heat), res$meta_degs)
  expect_equal(colnames(heat), res$design$sample)
  expect_false(anyNA(heat))
})

test_that("demo heatmap values are within-study z-scores of the right gene and sample", {
  res <- demo_result()
  heat <- heatmap_matrix(res)

  for (gene in utils::head(res$meta_degs, 25)) {
    for (study in res$studies) {
      samples <- res$design$sample[res$design$study == study]
      x <- res$norm_counts[gene, samples]
      expected <- if (stats::sd(x) > 0) (x - mean(x)) / stats::sd(x) else x * 0
      expect_equal(heat[gene, samples], expected, tolerance = 1e-10)
    }
  }
})

test_that("demo heatmap renders, including fewer DEGs than samples and chosen clusters", {
  skip_if_not_installed("heatmaply")
  skip_if_not_installed("pheatmap")
  res <- demo_result()
  heat <- heatmap_matrix(res)
  k_col <- length(res$studies)

  expect_no_error(plot_heatmap(heat, heatmap_row_clusters(nrow(heat)), k_col))
  expect_no_error(plot_heatmap(heat[1:3, ], heatmap_row_clusters(3), k_col))
  expect_no_error(plot_heatmap(heat, heatmap_row_clusters(nrow(heat), 5), k_col))

  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  expect_no_error(plot_heatmap_static(heat, res$design, heatmap_row_clusters(nrow(heat)), k_col))
  expect_no_error(plot_heatmap_static(heat, res$design, heatmap_row_clusters(nrow(heat), 4), k_col))
  grDevices::dev.off()
})

test_that("demo normalised counts and DE table match a direct DESeq2 run", {
  res <- demo_result()
  demo <- tma_demo_data()
  design_a <- demo$design[demo$design$study == "A", ]
  rownames(design_a) <- design_a$sample
  counts_a <- as.matrix(demo$counts[, design_a$sample])
  rownames(counts_a) <- demo$counts$Geneid
  dds <- suppressMessages(suppressWarnings(
    DESeq2::DESeqDataSetFromMatrix(counts_a, design_a, ~condition)
  ))
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  direct <- DESeq2::results(dds, alpha = 0.05)

  expect_equal(res$norm_counts[, design_a$sample], DESeq2::counts(dds, normalized = TRUE))
  expect_equal(res$de$A$log2_fc, direct$log2FoldChange)
  expect_equal(res$de$A$padj, direct$padj)
})

test_that("demo volcano mean fold change matches the per-gene average", {
  res <- demo_result()
  vol <- volcano_data(res)

  expect_equal(nrow(vol), length(res$genes))
  for (i in c(1, 4, 100, 1000, length(res$genes))) {
    expect_equal(vol$fc[i], mean(res$fold_changes[i, ], na.rm = TRUE))
  }
  expect_equal(vol$pval, res$meta$adjpval)
})

test_that("demo enrichment genes are the meta-significant, direction-consistent DEGs", {
  res <- demo_result()
  degs <- res$meta_degs
  consistent <- rownames(res$fc_select)[which(abs(res$fc_select$signFC) == 1)]

  expect_gt(length(degs), 0)
  expect_true(all(res$meta[degs, "adjpval"] <= 0.05))
  expect_true(all(degs %in% consistent))
  expect_lt(length(degs), length(consistent))
  expect_setequal(rownames(res$meta)[which(res$meta$significant == "yes")], degs)
})

test_that("demo meta-analysis reproduces the recorded DEGs", {
  res <- demo_result()
  ref <- read_demo_reference()
  same_versions <- as.character(utils::packageVersion("DESeq2")) == ref$deseq2 &&
    as.character(utils::packageVersion("metaRNASeq")) == ref$metarnaseq
  regenerate <- "If the change is expected, regenerate the reference with tests/testthat/fixtures/make-demo-reference.R."

  if (same_versions) {
    # same package versions: the result must be identical
    expect_setequal(res$meta_degs, ref$degs)
    expect_equal(sum(abs(res$fc_select$signFC) == 1, na.rm = TRUE), ref$consistent_genes)
  } else {
    # newer DESeq2/metaRNASeq: small numerical changes are allowed, a changed
    # pipeline is not (>= 95% overlap and <= 5% change in the number of DEGs)
    overlap <- length(intersect(res$meta_degs, ref$degs)) /
      length(union(res$meta_degs, ref$degs))
    size_change <- abs(length(res$meta_degs) - length(ref$degs)) / length(ref$degs)
    expect(overlap >= 0.95, sprintf(
      "DEG overlap with the reference (DESeq2 %s) is %.1f%%, below 95%%. %s",
      ref$deseq2, 100 * overlap, regenerate
    ))
    expect(size_change <= 0.05, sprintf(
      "%d DEGs vs %d in the reference (DESeq2 %s), a %.1f%% change. %s",
      length(res$meta_degs), length(ref$degs), ref$deseq2, 100 * size_change, regenerate
    ))
  }
})
