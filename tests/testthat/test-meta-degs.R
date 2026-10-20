test_that("select_meta_degs keeps meta-significant genes with a consistent direction", {
  fc_select <- data.frame(
    A = c(1, 1, 0, 1, 0),
    DE.meta = c(1, 1, 0, NA, 1),
    signFC = c(1, 0, -1, 1, -1),
    row.names = c("up", "conflict", "not_meta", "na_meta", "down")
  )
  expect_equal(select_meta_degs(fc_select), c("up", "down"))
})

test_that("select_meta_degs drops genes with a missing direction", {
  fc_select <- data.frame(DE.meta = c(1, 1), signFC = c(NA, -1),
                          row.names = c("na_sign", "down"))
  expect_equal(select_meta_degs(fc_select), "down")
})

test_that("select_meta_degs returns an empty vector when nothing is significant", {
  fc_select <- data.frame(DE.meta = c(0, 0), signFC = c(1, -1),
                          row.names = c("g1", "g2"))
  expect_identical(select_meta_degs(fc_select), character(0))
})
