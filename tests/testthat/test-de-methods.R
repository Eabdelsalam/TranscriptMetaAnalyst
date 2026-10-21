demo_study <- function(id = "A") {
  demo <- tma_demo_data()
  split_studies(demo$counts, demo$design)[[id]]
}

test_that("every DE method returns the standard table", {
  study <- demo_study()
  for (m in de_methods()) {
    de <- run_de(study$counts, study$design, method = m)
    expect_named(de$table, c("gene", "base_mean", "log2_fc", "lfc_se", "stat", "pvalue", "padj"))
    expect_equal(de$table$gene, rownames(study$counts))
    expect_equal(dim(de$norm_counts), dim(study$counts))
    expect_true(all(de$table$pvalue >= 0 & de$table$pvalue <= 1, na.rm = TRUE))
  }
  expect_error(run_de(study$counts, study$design, method = "t-test"))
})

test_that("run_edger matches a direct edgeR quasi-likelihood fit", {
  study <- demo_study()
  de <- run_edger(study$counts, study$design)
  group <- factor(study$design$condition)
  model <- stats::model.matrix(~group)
  y <- edgeR::estimateDisp(edgeR::normLibSizes(edgeR::DGEList(study$counts, group = group)), model)
  direct <- edgeR::glmQLFTest(edgeR::glmQLFit(y, model), coef = 2)$table
  expect_equal(de$table$log2_fc, direct$logFC)
  expect_equal(de$table$pvalue, direct$PValue)
  expect_equal(de$table$padj, stats::p.adjust(direct$PValue, "BH"))
})

test_that("run_limma_voom matches a direct limma-voom fit", {
  study <- demo_study()
  de <- run_limma_voom(study$counts, study$design)
  group <- factor(study$design$condition)
  model <- stats::model.matrix(~group)
  y <- edgeR::normLibSizes(edgeR::DGEList(study$counts, group = group))
  fit <- limma::eBayes(limma::lmFit(limma::voom(y, model), model))
  expect_equal(de$table$log2_fc, unname(fit$coefficients[, 2]))
  expect_equal(de$table$pvalue, unname(fit$p.value[, 2]))
  expect_equal(de$table$stat, unname(fit$t[, 2]))
})

test_that("the DE methods agree on the direction of strongly changed genes", {
  study <- demo_study()
  fc <- sapply(de_methods(), function(m) run_de(study$counts, study$design, method = m)$table$log2_fc)
  strong <- abs(fc[, "DESeq2"]) > 2
  expect_gt(mean(sign(fc[strong, "edgeR"]) == sign(fc[strong, "DESeq2"])), 0.95)
  expect_gt(mean(sign(fc[strong, "limma-voom"]) == sign(fc[strong, "DESeq2"])), 0.95)
})

test_that("plots use the same normalised counts whatever the DE method", {
  study <- demo_study()
  deseq <- run_deseq2(study$counts, study$design)$norm_counts
  expect_equal(run_edger(study$counts, study$design)$norm_counts, deseq)
  expect_equal(run_limma_voom(study$counts, study$design)$norm_counts, deseq)
})

test_that("run_de_studies takes one method or one per study", {
  demo <- tma_demo_data()
  studies <- split_studies(demo$counts, demo$design)
  serial <- BiocParallel::SerialParam()
  mixed <- run_de_studies(studies, methods = c(A = "edgeR", B = "DESeq2", C = "limma-voom"),
                          BPPARAM = serial)
  expect_named(mixed, c("A", "B", "C"))
  expect_equal(mixed$A, run_edger(studies$A$counts, studies$A$design))
  expect_equal(mixed$C, run_limma_voom(studies$C$counts, studies$C$design))
  expect_error(run_de_studies(studies, methods = c(A = "edgeR"), BPPARAM = serial), "study: B, C")
  expect_error(run_de_studies(studies, methods = "nope", BPPARAM = serial), "Unknown DE method")
})
