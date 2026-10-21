test_that("read_array_matrix finds the ID, symbol and sample columns", {
  f <- withr::local_tempfile(fileext = ".txt")
  writeLines(c("probe\tGene Symbol\tS1\tS2", "p1\tgeneA\t5.1\t6.2", "p2\t\t7.3\t7.0"), f)
  m <- read_array_matrix(f)
  expect_named(m, c("ID", "Symbol", "S1", "S2"))
  expect_equal(m$Symbol, c("geneA", ""))
  g <- withr::local_tempfile(fileext = ".txt")
  writeLines(c("ID\tS1\tS2", "geneA\t5\t6"), g)
  expect_true(is.na(read_array_matrix(g)$Symbol))
})

test_that("log2 is detected like GEO2R and applied only when needed", {
  linear <- matrix(c(15, 230, 4100, 61000, 900, 12), 3)
  logged <- log2(linear)
  expect_true(needs_log2(linear))
  expect_false(needs_log2(logged))
  out <- log2_if_needed(cbind(linear, c(-1, 0, 8)))
  expect_true(out$transformed)
  expect_equal(out$values[, 1:2], log2(linear))
  expect_equal(out$values[, 3], c(NA, NA, 3))
  expect_false(log2_if_needed(logged)$transformed)
})

test_that("collapse_probes keeps the probe with the highest mean per gene", {
  m <- rbind(p1 = c(1, 1), p2 = c(5, 7), p3 = c(9, 9), p4 = c(3, 3), p5 = c(8, 2), p6 = c(4, 4))
  genes <- c("A", "A", "", "B", "B /// C", "B")
  out <- collapse_probes(m, genes)
  expect_equal(rownames(out$values), c("A", "B"))
  expect_equal(unname(out$values["A", ]), c(5, 7))
  expect_equal(unname(out$values["B", ]), c(4, 4))
  expect_equal(out$report, list(probes = 6, no_gene = 1, multiple_genes = 1, genes = 2))
  expect_equal(collapse_probes(m[1:2, ], c("x", "y"))$report$genes, 2)
})

test_that("prepare_array_studies matches files to studies and intersects genes", {
  demo <- tma_microarray_demo_data()
  prepared <- prepare_array_studies(demo$matrices, demo$design)
  expect_named(prepared$studies, unique(demo$design$study))
  common <- attr(prepared$report, "common_genes")
  expect_equal(common, 3358)
  for (s in prepared$studies) {
    expect_equal(nrow(s$values), common)
    expect_equal(colnames(s$values), s$design$sample)
    expect_equal(rownames(s$values), rownames(prepared$studies[[1]]$values))
  }
  expect_equal(prepared$report$no_gene, rep(20, 3))
  expect_match(prepared$notes, "3358 genes are measured in all 3 studies", all = FALSE)
  # the collapsed value is the gene's probe with the highest mean
  raw <- demo$matrices$GSE63352
  gene <- rownames(prepared$studies$GSE63352$values)[10]
  probes <- as.matrix(raw[raw$Symbol == gene, prepared$studies$GSE63352$design$sample])
  expect_equal(unname(prepared$studies$GSE63352$values[gene, ]),
               unname(probes[which.max(rowMeans(probes)), ]))
})

test_that("prepare_array_studies log2-transforms linear data and explains input problems", {
  demo <- tma_microarray_demo_data()
  linear <- demo$matrices
  samples <- demo$design$sample[demo$design$study == "GSE63352"]
  linear$GSE63352[, samples] <- 2^linear$GSE63352[, samples]
  prepared <- prepare_array_studies(linear, demo$design)
  expect_true(prepared$report$log2_applied[prepared$report$study == "GSE63352"])
  expect_match(prepared$notes, "Study GSE63352: the values looked like linear intensities", all = FALSE)
  expect_equal(prepared$studies$GSE63352$values,
               prepare_array_studies(demo$matrices, demo$design)$studies$GSE63352$values)

  without <- demo$matrices[names(demo$matrices) != "GSE68250"]
  expect_error(prepare_array_studies(without, demo$design), "No expression file was given for study GSE68250")
  merged <- demo$matrices$GSE63352
  merged <- cbind(merged, demo$matrices$GSE68250[match(merged$ID, demo$matrices$GSE68250$ID), -(1:2)])
  expect_error(prepare_array_studies(list(merged = merged), demo$design), "contains samples of several studies")
  expect_error(prepare_array_studies(list(x = data.frame(ID = "a", Symbol = "a", S9 = 1)), demo$design),
               "has no column that matches a sample")
})

test_that("genes are intersected across platforms with different gene sets", {
  demo <- tma_microarray_demo_data()
  # a second "platform" that lacks some genes and measures each gene once
  other <- demo$matrices$GSE68250
  other <- other[other$Symbol != "" & !duplicated(other$Symbol), ]
  removed <- other$Symbol[1:100]
  other <- other[-(1:100), ]
  other$ID <- paste0("other_", other$ID)
  matrices <- demo$matrices
  matrices$GSE68250 <- other
  prepared <- prepare_array_studies(matrices, demo$design)
  expect_equal(attr(prepared$report, "common_genes"), 3358 - 100)
  expect_false(any(removed %in% rownames(prepared$studies$GSE82120$values)))
})
