serial <- BiocParallel::SerialParam()
array_cache <- new.env()
array_demo_analysis <- function() {
  if (is.null(array_cache$analysis)) {
    demo <- tma_microarray_demo_data()
    array_cache$analysis <- run_array_analysis(demo$matrices, demo$design, reference = "HC", BPPARAM = serial)
  }
  array_cache$analysis
}

test_that("run_limma_array matches a direct limma fit with contrasts and a trended eBayes", {
  prepared <- prepare_array_studies(tma_microarray_demo_data()$matrices, tma_microarray_demo_data()$design)
  study <- prepared$studies$GSE63352
  de <- run_limma_array(study$values, study$design, reference = "HC")
  expect_named(de$tables, c("LC_24h vs HC", "LC_3h vs HC"))

  condition <- factor(study$design$condition, levels = c("HC", "LC_24h", "LC_3h"))
  model <- stats::model.matrix(~ 0 + condition)
  cm <- limma::makeContrasts(conditionLC_3h - conditionHC, levels = model)
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(study$values, model), cm), trend = TRUE)
  direct <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none")
  expect_equal(de$tables[["LC_3h vs HC"]]$log2_fc, direct$logFC)
  expect_equal(de$tables[["LC_3h vs HC"]]$pvalue, direct$P.Value)
  expect_equal(de$tables[["LC_3h vs HC"]]$padj, direct$adj.P.Val)
  expect_equal(de$norm_counts, study$values)
})

test_that("run_limma_array adds batch to the model when asked", {
  prepared <- prepare_array_studies(tma_microarray_demo_data()$matrices, tma_microarray_demo_data()$design)
  study <- prepared$studies$GSE82120
  study$design$batch <- rep(c("b1", "b2"), 6)
  with_batch <- run_limma_array(study$values, study$design, use_batch = TRUE, reference = "HC")
  without <- run_limma_array(study$values, study$design, reference = "HC")
  expect_false(isTRUE(all.equal(with_batch$table$pvalue, without$table$pvalue)))
  model <- stats::model.matrix(~ 0 + factor(study$design$condition, levels = c("HC", "LC_24h")) + factor(study$design$batch))
  colnames(model) <- c("HC", "LC", "b2")
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(study$values, model),
                                            limma::makeContrasts(LC - HC, levels = model)), trend = TRUE)
  expect_equal(with_batch$table$pvalue, unname(fit$p.value[, 1]))
})

test_that("the microarray demo is meta-analysed per comparison with the replicate rule", {
  analysis <- array_demo_analysis()
  expect_s3_class(analysis, "tma_analysis")
  expect_equal(analysis$contrasts, c("LC_24h vs HC", "LC_3h vs HC"))
  expect_equal(analysis$results[["LC_24h vs HC"]]$studies, c("GSE82120", "GSE63352", "GSE68250"))
  expect_equal(analysis$results[["LC_3h vs HC"]]$studies, c("GSE63352", "GSE68250"))
  expect_match(analysis_messages(analysis), "Study GSE82120 has no samples of 'LC_3h'")
  res <- analysis$results[["LC_24h vs HC"]]
  expect_equal(res$scale, "log2")
  expect_equal(unname(res$de_methods), rep("limma", 3))
  expect_equal(nrow(res$norm_counts), 3358)
  expect_equal(attr(analysis$preparation$report, "common_genes"), 3358)
  # the bicarbonate transporters are induced at low CO2
  up <- c("slr1512", "slr0040", "sll1733")
  expect_true(all(up %in% res$meta_degs))
  expect_true(all(mean_fold_change(res$fold_changes[up, ]) > 1))
})

test_that("plots use log2 expression values for microarray results", {
  res <- array_demo_analysis()$results[["LC_24h vs HC"]]
  p <- plot_gene_counts(res, "slr1512")
  expect_equal(p$labels$y, "log2 expression")
  expect_equal(p$data$count, unname(res$norm_counts["slr1512", res$design$sample]))
  expect_false(any(vapply(p$scales$scales, function(s) inherits(s$trans, "transform") && s$trans$name == "log-10", logical(1))))
  expect_s3_class(plot_pca(res), "ggplot")
  heat <- heatmap_matrix(res)
  expect_equal(rownames(heat), res$meta_degs)
  expect_no_error(ggplot2::ggplot_build(plot_volcano(volcano_data(res))))
})

test_that("microarray batch handling: model term or none, not ComBat-seq", {
  demo <- tma_microarray_demo_data()
  expect_error(run_array_analysis(demo$matrices, demo$design, batch_method = "combat_seq"),
               "ComBat-seq corrects RNA-seq counts")
  design <- demo$design
  design$batch <- ave(design$sample, design$study, FUN = function(x) rep(c("b1", "b2"), length.out = length(x)))
  analysis <- run_array_analysis(demo$matrices, design, reference = "HC", batch_method = "model",
                                 batch_plots = TRUE, BPPARAM = serial)
  res <- analysis$results[["LC_24h vs HC"]]
  expect_true("applied" %in% res$batch$status)
  expect_equal(dim(res$batch$corrected_log), dim(res$norm_counts))
  expect_s3_class(plot_batch_pca(res), "ggplot")
})
