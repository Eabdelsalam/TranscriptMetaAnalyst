serial <- BiocParallel::SerialParam()

test_that("default_reference prefers a control-like name, else the first alphabetically", {
  expect_equal(default_reference(c("Treated", "Control", "Heat")), "Control")
  expect_equal(default_reference(c("Heat", "WT", "Mutant")), "WT")
  expect_equal(default_reference(c("treated", "mock")), "mock")
  expect_equal(default_reference(c("Heat", "Cold")), "Cold")
  expect_equal(default_reference(tma_demo_data()$design$condition), "Cont")
})

test_that("every condition is compared with the reference, each in its own result", {
  sim <- multi_data_set()
  analysis <- run_tma_analysis(sim$counts, sim$design, reference = "Cont", BPPARAM = serial)
  expect_s3_class(analysis, "tma_analysis")
  expect_equal(analysis$contrasts, c("Cold vs Cont", "Heat vs Cont"))
  expect_named(analysis$results, analysis$contrasts)
  heat <- analysis$results[["Heat vs Cont"]]
  expect_equal(heat$contrast, "Heat vs Cont")
  expect_setequal(heat$design$condition, c("Cont", "Heat"))
  expect_equal(ncol(heat$norm_counts), 18L)
  # each contrast finds its own responding genes
  expect_gt(mean(sim$heat_genes %in% heat$meta_degs), 0.8)
  expect_lt(mean(sim$cold_genes %in% heat$meta_degs), 0.1)
  cold <- analysis$results[["Cold vs Cont"]]
  expect_gt(mean(sim$cold_genes %in% cold$meta_degs), 0.8)
  expect_true(all(analysis$report$used))
  expect_length(analysis_messages(analysis), 0)
  expect_output(print(analysis), "Heat vs Cont: 3 studies")
})

test_that("each study is fitted once and the contrasts come from that fit", {
  sim <- multi_data_set()
  fits <- 0
  original <- run_de
  local_mocked_bindings(run_de = function(...) {
    fits <<- fits + 1
    original(...)
  })
  analysis <- run_tma_analysis(sim$counts, sim$design, reference = "Cont", BPPARAM = serial)
  expect_equal(fits, 3)

  study <- filter_low_counts(split_studies(sim$counts, sim$design))$S2
  d <- study$design
  d$condition <- factor(d$condition, levels = c("Cont", "Cold", "Heat"))
  dds <- DESeq2::DESeq(suppressMessages(DESeq2::DESeqDataSetFromMatrix(study$counts, d, ~condition)),
                       quiet = TRUE)
  direct <- DESeq2::results(dds, contrast = c("condition", "Heat", "Cont"), alpha = 0.05)
  expect_equal(analysis$results[["Heat vs Cont"]]$de$S2$log2_fc, direct$log2FoldChange)
  expect_equal(analysis$results[["Heat vs Cont"]]$de$S2$padj, direct$padj)
})

test_that("the edgeR and limma-voom contrasts use the right coefficient", {
  sim <- multi_data_set()
  study <- split_studies(sim$counts, sim$design)$S1
  for (m in c("edgeR", "limma-voom")) {
    de <- run_de(study$counts, study$design, method = m, reference = "Cont")
    expect_named(de$tables, c("Cold vs Cont", "Heat vs Cont"))
    heat_fc <- de$tables[["Heat vs Cont"]]$log2_fc
    cold_fc <- de$tables[["Cold vs Cont"]]$log2_fc
    expect_gt(mean(heat_fc[1:80]), 1)
    expect_lt(abs(mean(heat_fc[81:160])), 0.3)
    expect_lt(mean(cold_fc[81:160]), -1)
  }
})

test_that("a condition with one sample in a study is excluded from that contrast only", {
  sim <- multi_data_set()
  design <- sim$design
  drop <- design$sample %in% c("S2_5", "S2_6") # S2 keeps one Heat sample
  counts <- sim$counts[, !names(sim$counts) %in% design$sample[drop]]
  design <- design[!drop, ]
  analysis <- run_tma_analysis(counts, design, reference = "Cont", BPPARAM = serial)

  expect_equal(analysis$results[["Heat vs Cont"]]$studies, c("S1", "S3"))
  expect_equal(analysis$results[["Cold vs Cont"]]$studies, c("S1", "S2", "S3"))
  expect_match(analysis_messages(analysis), "Study S2 has only 1 sample of 'Heat' \\(at least 2 are needed\\), so it is not used for 'Heat vs Cont'", all = FALSE)
  # the single Heat sample is not part of S2's fit
  expect_false("S2_4" %in% analysis$results[["Cold vs Cont"]]$design$sample)
})

test_that("a study without 2 reference samples is not used at all", {
  sim <- multi_data_set()
  design <- sim$design
  design$condition[design$sample %in% c("S3_1", "S3_2")] <- "Other"
  analysis <- run_tma_analysis(sim$counts, design, reference = "Cont", BPPARAM = serial)
  expect_false("S3" %in% analysis$results[["Heat vs Cont"]]$studies)
  expect_match(analysis_messages(analysis), "Study S3 has only 1 sample of the reference condition 'Cont'", all = FALSE)
  # "Other" has 2 samples only in S3, which cannot be used: not meta-analysed
  expect_false(analysis$summary$analysed[analysis$summary$contrast == "Other vs Cont"])
  expect_match(analysis_messages(analysis), "'Other vs Cont' was not meta-analysed: no study has both conditions", all = FALSE)
})

test_that("a contrast with fewer than 2 studies is not meta-analysed, and no contrast is an error", {
  sim <- multi_data_set()
  design <- sim$design
  # only S1 keeps two Cold samples
  design$condition[design$sample %in% c("S2_8", "S2_9", "S3_8", "S3_9")] <- "Heat"
  analysis <- run_tma_analysis(sim$counts, design, reference = "Cont", BPPARAM = serial)
  expect_null(analysis$results[["Cold vs Cont"]])
  expect_match(analysis_messages(analysis), "only 1 study \\(S1\\) has both conditions", all = FALSE)
  expect_s3_class(analysis$results[["Heat vs Cont"]], "tma_result")

  one_study <- sim$design$study == "S1"
  expect_error(run_tma_analysis(sim$counts[, c("Geneid", sim$design$sample[one_study])],
                                sim$design[one_study, ], BPPARAM = serial),
               "No contrast can be meta-analysed")
})

test_that("run_tma_pipeline is for two conditions", {
  sim <- multi_data_set()
  expect_error(run_tma_pipeline(sim$counts, sim$design, BPPARAM = serial), "use run_tma_analysis")
  res <- demo_result()
  expect_equal(res$contrast, "Treated vs Cont")
  expect_true(all(res$report$used))
})

test_that("the RNA-seq tab shows one results tab per contrast and the left-out studies", {
  sim <- multi_data_set()
  design <- sim$design
  drop <- design$sample %in% c("S2_5", "S2_6")
  counts <- sim$counts[, !names(sim$counts) %in% design$sample[drop]]
  design <- design[!drop, ]
  analysis <- run_tma_analysis(counts, design, reference = "Cont", BPPARAM = serial)
  local_mocked_bindings(run_tma_analysis = function(...) analysis)
  counts_file <- withr::local_tempfile(fileext = ".txt")
  design_file <- withr::local_tempfile(fileext = ".txt")
  utils::write.table(counts, counts_file, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(design, design_file, sep = "\t", quote = FALSE, row.names = FALSE)
  mock_organism_lists()
  app <- shiny::shinyAppDir(system.file("app", package = "TranscriptMetaAnalyst"))
  shiny::testServer(app, {
    session$setInputs(rnaseq_counts = list(datapath = counts_file),
                      rnaseq_design = list(datapath = design_file),
                      rnaseq_analyze = 1)
    for (i in 1:200) {
      later::run_now(0.05)
      session$flushReact()
      if (rnaseq_task$status() != "running") break
    }
    html <- as.character(output$rnaseq_results$html)
    expect_match(html, "Cold vs Cont")
    expect_match(html, "Heat vs Cont")
    expect_match(html, "Study S2 has only 1 sample of .Heat.")
    expect_match(output$rnaseq_summary, "Heat vs Cont")
    # each tab's module renders its own contrast
    expect_match(as.character(output$`rnaseq_c2-num_studies`$html), "Heat vs Cont")
    expect_match(as.character(output$`rnaseq_c2-num_studies`$html), "2 studies")
    expect_match(as.character(output$`rnaseq_c1-num_studies`$html), "3 studies")
  })
})
