test_that("validate_design accepts the demo design and reports problems", {
  design <- tma_demo_data()$design
  expect_identical(validate_design(design), character(0))
  expect_match(validate_design(design[, c("sample", "condition")]), "'study'")
  expect_match(validate_design(design[, c("study", "condition")]), "'sample'")

  three <- design
  three$condition[1] <- "Other"
  expect_match(validate_design(three), "two unique conditions")

  dup <- design
  dup$sample[2] <- dup$sample[1]
  expect_match(validate_design(dup), "duplicated sample names")
})

test_that("validate_counts checks the Geneid column", {
  counts <- tma_demo_data()$counts
  expect_identical(validate_counts(counts), character(0))
  expect_match(validate_counts(counts[, -1]), "'Geneid'")
  dup <- counts
  dup$Geneid[2] <- dup$Geneid[1]
  expect_match(validate_counts(dup), "duplicated gene names")
})

test_that("validate_match needs a counts column for every design sample", {
  demo <- tma_demo_data()
  expect_identical(validate_match(demo$counts, demo$design), character(0))
  expect_match(validate_match(demo$counts[, -8], demo$design), "corresponding counts column")
})

test_that("the example files are bundled", {
  expect_true(file.exists(tma_example_file("counts")))
  demo <- tma_demo_data()
  expect_equal(dim(demo$design), c(12L, 3L))
  expect_true("Geneid" %in% names(demo$counts))
})
