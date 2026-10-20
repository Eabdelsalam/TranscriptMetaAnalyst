serial <- BiocParallel::SerialParam()

test_that("batch_status explains when batch cannot be used", {
  base <- data.frame(condition = rep(c("c", "t"), each = 3), batch = rep(c("b1", "b2"), 3))
  expect_equal(batch_status(base, "model"), "applied")
  expect_equal(batch_status(base, "combat_seq"), "applied")
  expect_equal(batch_status(base, "none"), "ignored")
  expect_equal(batch_status(base[, "condition", drop = FALSE], "model"), "no batch column")
  expect_equal(batch_status(transform(base, batch = "b1"), "model"), "single batch")
  expect_equal(batch_status(transform(base, batch = condition), "model"), "confounded with condition")
  four <- data.frame(condition = c("c", "c", "t", "t"), batch = c("b1", "b2", "b1", "b2"))
  expect_equal(batch_status(four[1:3, ], "model"), "no residual degrees of freedom")
  one <- data.frame(condition = c("c", "c", "t", "t"), batch = c("b1", "b2", "b2", "b2"))
  expect_equal(batch_status(one, "combat_seq"), "a batch has only one sample")
})

test_that("DE with batch fits ~ batch + condition and reports the condition effect", {
  sim <- batch_data_set()
  study <- split_studies(sim$counts, sim$design)$S1
  de <- run_deseq2(study$counts, study$design, use_batch = TRUE)
  d <- study$design
  d$batch <- factor(d$batch)
  d$condition <- factor(d$condition)
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(study$counts, d, ~ batch + condition))
  direct <- DESeq2::results(DESeq2::DESeq(dds, quiet = TRUE), name = "condition_Treated_vs_Cont")
  expect_equal(de$table$log2_fc, direct$log2FoldChange)
  expect_equal(de$table$pvalue, direct$pvalue)

  model <- stats::model.matrix(~ factor(d$batch) + d$condition)
  y <- edgeR::estimateDisp(edgeR::normLibSizes(edgeR::DGEList(study$counts)), model)
  fit <- edgeR::glmQLFit(y, model)
  direct_edger <- edgeR::glmQLFTest(fit, coef = ncol(model))$table
  expect_equal(run_edger(study$counts, study$design, use_batch = TRUE)$table$pvalue,
               direct_edger$PValue)
  v <- limma::voom(edgeR::normLibSizes(edgeR::DGEList(study$counts)), model)
  direct_limma <- limma::eBayes(limma::lmFit(v, model))
  expect_equal(run_limma_voom(study$counts, study$design, use_batch = TRUE)$table$log2_fc,
               unname(direct_limma$coefficients[, ncol(model)]))
})

test_that("putting batch in the model finds more of the true DEGs", {
  sim <- batch_data_set()
  with_batch <- run_tma_pipeline(sim$counts, sim$design, BPPARAM = serial)
  without <- run_tma_pipeline(sim$counts, sim$design, batch_method = "none", BPPARAM = serial)
  expect_equal(unname(with_batch$batch$status), c("applied", "applied"))
  expect_equal(unname(without$batch$status), c("ignored", "ignored"))
  tp <- function(res) sum(res$meta_degs %in% sim$de_genes)
  expect_gt(tp(with_batch), tp(without))
})

test_that("ComBat-seq corrects the counts before a ~ condition analysis", {
  sim <- batch_data_set()
  study <- filter_low_counts(split_studies(sim$counts, sim$design))$S1
  corrected <- combat_seq_counts(study$counts, study$design)
  expect_equal(dimnames(corrected), dimnames(study$counts))
  expect_false(isTRUE(all.equal(corrected, study$counts)))

  res <- run_tma_pipeline(sim$counts, sim$design, batch_method = "combat_seq", BPPARAM = serial)
  expect_equal(unname(res$batch$status), c("applied", "applied"))
  expect_equal(res$de$S1, run_deseq2(corrected, study$design)$table)
  # the counts kept for the plots are the uncorrected ones
  expect_equal(res$norm_counts[, colnames(study$counts)],
               run_deseq2(study$counts, study$design)$norm_counts)
})

test_that("a study whose batch is confounded with condition falls back to ~ condition", {
  sim <- batch_data_set()
  design <- sim$design
  design$batch[design$study == "S2"] <- design$condition[design$study == "S2"]
  res <- run_tma_pipeline(sim$counts, design, BPPARAM = serial)
  expect_equal(res$batch$status, c(S1 = "applied", S2 = "confounded with condition"))
  study2 <- filter_low_counts(split_studies(sim$counts, design))$S2
  expect_equal(res$de$S2, run_deseq2(study2$counts, study2$design)$table)
})

test_that("the corrected values remove the batch difference but keep the condition effect", {
  sim <- batch_data_set()
  res <- run_tma_pipeline(sim$counts, sim$design, BPPARAM = serial)
  s1 <- sim$design[sim$design$study == "S1", ]
  b2 <- s1$sample[s1$batch == "b2"]
  b1 <- s1$sample[s1$batch == "b1"]
  genes <- setdiff(sim$batch_genes, sim$de_genes)
  gap <- function(m) mean(rowMeans(m[genes, b2]) - rowMeans(m[genes, b1]))
  expect_gt(gap(log2(res$norm_counts + 1)), 1.5)
  expect_lt(abs(gap(res$batch$corrected_log)), 0.1)

  treated <- s1$sample[s1$condition == "Treated"]
  cont <- s1$sample[s1$condition == "Cont"]
  de_gap <- mean(rowMeans(res$batch$corrected_log[sim$de_genes, treated]) -
                   rowMeans(res$batch$corrected_log[sim$de_genes, cont]))
  expect_gt(de_gap, 0.6) # log2(1.8) = 0.85
})

test_that("results without a batch column have no batch information", {
  expect_null(demo_result()$batch)
})

test_that("batch diagnostics plots build", {
  sim <- batch_data_set()
  res <- run_tma_pipeline(sim$counts, sim$design, BPPARAM = serial)
  pca <- plot_batch_pca(res)
  rle <- plot_rle(res)
  expect_s3_class(pca, "ggplot")
  expect_equal(nrow(pca$data), 2 * nrow(sim$design))
  expect_setequal(levels(pca$data$stage), c("Before", "After"))
  expect_no_error(ggplot2::ggplot_build(rle))
  # RLE medians move towards 0 after the correction
  med <- tapply(abs(rle$data$middle), rle$data$stage, mean)
  expect_lt(med[["After"]], med[["Before"]])
})

test_that("validate_design rejects empty batch values", {
  design <- tma_demo_data()$design
  design$batch <- c(rep("b1", 11), "")
  expect_match(validate_design(design), "'batch' column")
  design$batch[12] <- "b2"
  expect_identical(validate_design(design), character(0))
})

test_that("batch_plots switches the PCA and heatmap to the corrected values", {
  sim <- batch_data_set()
  plain <- run_tma_pipeline(sim$counts, sim$design, BPPARAM = serial)
  corrected <- run_tma_pipeline(sim$counts, sim$design, batch_plots = TRUE, BPPARAM = serial)
  expect_equal(plain$meta_degs, corrected$meta_degs)

  s1 <- sim$design$sample[sim$design$study == "S1"]
  gene <- corrected$meta_degs[1]
  x <- corrected$batch$corrected_log[gene, s1]
  expect_equal(heatmap_matrix(corrected)[gene, s1], (x - mean(x)) / stats::sd(x))
  expect_false(isTRUE(all.equal(heatmap_matrix(plain), heatmap_matrix(corrected))))
  expect_false(isTRUE(all.equal(plot_pca(plain)$data$PC1, plot_pca(corrected)$data$PC1)))
})

test_that("the results module shows the batch status of each study", {
  sim <- batch_data_set()
  design <- sim$design
  design$batch[design$study == "S2"] <- "b1"
  res <- run_tma_pipeline(sim$counts, design, BPPARAM = serial)
  shiny::testServer(mod_results_server, args = list(result = shiny::reactive(res),
                                                    organism = list(go = 1, kegg = "x")), {
    status <- output$batch_status
    expect_match(status, "S1")
    expect_match(status, "Batch in the DE model")
    expect_match(status, "not used \\(single batch\\)")
    expect_s3_class(batch_pca(), "ggplot")
    expect_s3_class(batch_rle(), "ggplot")
  })
  # no batch column: no diagnostics
  shiny::testServer(mod_results_server, args = list(result = shiny::reactive(demo_result()),
                                                    organism = list(go = 1, kegg = "x")), {
    expect_error(batch_pca(), class = "shiny.silent.error")
  })
})
