test_that("mean_fold_change averages each gene across studies", {
  fc <- list(A = c(1, -2, 0.5), B = c(3, -4, NA), C = c(2, 0, 1.5))
  expect_equal(mean_fold_change(fc), c(2, -2, 1))
})

test_that("the old matrix(unlist(...), byrow = TRUE) mean fold change is wrong", {
  fc <- list(A = c(1, -2, 0.5), B = c(3, -4, NA), C = c(2, 0, 1.5))
  old <- rowMeans(matrix(unlist(fc), ncol = 3, byrow = TRUE), na.rm = TRUE)
  expect_false(isTRUE(all.equal(old, mean_fold_change(fc))))
})
