# Two studies whose samples and gene order differ, so any value scrambling shows up
make_studies <- function() {
  a <- matrix(c(1, 2, 3,
                10, 20, 30,
                5, 5, 5,
                7, 8, 9),
              nrow = 4, byrow = TRUE,
              dimnames = list(c("g1", "g2", "g3", "g4"), c("a1", "a2", "a3")))
  b <- matrix(c(100, 300,
                40, 60,
                2, 4,
                11, 11),
              nrow = 4, byrow = TRUE,
              dimnames = list(c("g4", "g2", "g1", "g3"), c("b1", "b2")))
  list(A = a, B = b)
}

test_that("combine_study_matrices keeps every value with its gene and sample", {
  studies <- make_studies()
  combined <- combine_study_matrices(studies)

  expect_equal(dim(combined), c(4L, 5L))
  expect_equal(rownames(combined), c("g1", "g2", "g3", "g4"))
  expect_equal(colnames(combined), c("a1", "a2", "a3", "b1", "b2"))
  for (s in names(studies)) {
    m <- studies[[s]]
    expect_equal(combined[rownames(m), colnames(m)], m)
  }
})

test_that("the old matrix(unlist(...), byrow = TRUE) construction scrambles values", {
  studies <- make_studies()
  studies$B <- studies$B[rownames(studies$A), ]
  old <- matrix(unlist(studies), ncol = 5, byrow = TRUE)
  expect_false(isTRUE(all.equal(old, unname(combine_study_matrices(studies)))))
})

test_that("combine_study_matrices orders columns by the design samples", {
  combined <- combine_study_matrices(make_studies(), samples = c("b2", "a1", "b1", "a3", "a2"))
  expect_equal(colnames(combined), c("b2", "a1", "b1", "a3", "a2"))
  expect_equal(combined["g4", "b2"], 300)
})

test_that("combine_study_matrices rejects studies with different genes", {
  studies <- make_studies()
  studies$B <- studies$B[1:3, ]
  expect_error(combine_study_matrices(studies), "same genes")
})

test_that("scale_deg_counts z-scores each gene within each study", {
  heat <- scale_deg_counts(make_studies(), genes = c("g2", "g1"))

  expect_equal(rownames(heat), c("g2", "g1"))
  expect_equal(unname(heat["g1", c("a1", "a2", "a3")]), c(-1, 0, 1))
  expect_equal(unname(heat["g2", c("b1", "b2")]), c(-1, 1) / sqrt(2))
})

test_that("scale_deg_counts gives 0, not NaN, for constant genes", {
  heat <- scale_deg_counts(make_studies(), genes = c("g3", "g1"))
  expect_false(anyNA(heat))
  expect_equal(unname(heat["g3", c("a1", "a2", "a3")]), c(0, 0, 0))
})

test_that("scale_deg_counts keeps a matrix for a single gene", {
  heat <- scale_deg_counts(make_studies(), genes = "g1")
  expect_true(is.matrix(heat))
  expect_equal(dim(heat), c(1L, 5L))
})

test_that("heatmap_row_clusters never exceeds the number of genes", {
  expect_equal(heatmap_row_clusters(500), 2L)
  expect_equal(heatmap_row_clusters(1), 1L)
  expect_equal(heatmap_row_clusters(0), 1L)
  expect_lte(heatmap_row_clusters(3, k = 12), 3L)
})

test_that("heatmap_row_clusters uses the user's choice, capped by the number of genes", {
  expect_equal(heatmap_row_clusters(100, k = 5), 5L)
  expect_equal(heatmap_row_clusters(100, k = 1), 1L)
  expect_equal(heatmap_row_clusters(4, k = 10), 4L)
  expect_equal(heatmap_row_clusters(100, k = 3.7), 3L)
})

test_that("heatmap_row_clusters falls back to 2 for an empty or invalid choice", {
  expect_equal(heatmap_row_clusters(100, k = NA), 2L)
  expect_equal(heatmap_row_clusters(100, k = NULL), 2L)
  expect_equal(heatmap_row_clusters(100, k = 0), 2L)
  expect_equal(heatmap_row_clusters(100, k = "abc"), 2L)
})
