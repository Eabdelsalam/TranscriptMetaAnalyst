test_that("split_studies orders columns by the design and ignores annotation columns", {
  counts <- data.frame(Geneid = c("g1", "g2"), Length = c(100, 200),
                       s2 = c(5, 6), s1 = c(1, 2), s3 = c(9, 9))
  design <- data.frame(sample = c("s1", "s2", "s3"), study = c("A", "A", "B"),
                       condition = c("c", "t", "c"))
  studies <- split_studies(counts, design)

  expect_named(studies, c("A", "B"))
  expect_equal(colnames(studies$A$counts), c("s1", "s2"))
  expect_equal(studies$A$counts["g2", "s2"], 6)
  expect_equal(rownames(studies$B$design), "s3")
})

test_that("run_fisher_meta matches fishercomb and drops filtered p-values", {
  p <- matrix(c(0.001, 0.2, 0.5, 0.04, 0.004, 0.3, 0.9, 0.01), ncol = 2)
  padj <- p
  padj[2, 1] <- NA
  meta <- run_fisher_meta(p, padj)

  manual <- p
  manual[2, 1] <- NA
  fish <- metaRNASeq::fishercomb(list(manual[, 1], manual[, 2]), BHth = 0.05)
  expect_equal(meta$rawpval, fish$rawpval)
  expect_equal(meta$adjpval, fish$adjpval)
  # the filtered p-value is left out, so the result differs from using it
  unfiltered <- run_fisher_meta(p, p)
  expect_false(isTRUE(all.equal(meta$rawpval[2], unfiltered$rawpval[2])))
})

test_that("consistent_direction flags conflicting and missing fold changes", {
  fc <- matrix(c(1, -2, 0.5, NA,
                 3, -1, -0.5, 1), ncol = 2)
  expect_equal(consistent_direction(fc), c(1, -1, 0, NA))
})

test_that("run_tma_pipeline returns a complete tma_result on the demo data", {
  res <- demo_result()

  expect_s3_class(res, "tma_result")
  expect_equal(res$studies, c("A", "B", "C"))
  expect_equal(dim(res$norm_counts), c(length(res$genes), 12L))
  expect_equal(dim(res$fold_changes), c(length(res$genes), 3L))
  expect_equal(nrow(res$meta), length(res$genes))
  expect_named(res$study_degs, res$studies)
  expect_named(res$idd_irr, c("DE", "IDD", "Loss", "IDR", "IRR"))
  expect_equal(unname(res$idd_irr["DE"]), length(res$meta_degs))
  expect_output(print(res), "meta-analysis DEGs")
})

test_that("study_significant applies the padj and fold-change thresholds", {
  res <- demo_result()
  sig <- study_significant(res, "B")
  expect_true(all(sig$padj < 0.05 & abs(sig$log2_fc) > 1))
  all_b <- res$de$B
  expect_equal(nrow(sig), sum(all_b$padj < 0.05 & abs(all_b$log2_fc) > 1, na.rm = TRUE))
})
