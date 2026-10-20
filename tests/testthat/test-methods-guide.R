test_that("the methods guide states the scope, the requirements and the ID types", {
  html <- as.character(methods_guide_ui())
  expect_match(html, "gene-level")
  expect_match(html, "isoform")
  expect_match(html, "At least 2 samples")
  expect_match(html, "At least 2 studies per comparison")
  for (label in names(id_types())) expect_match(html, label, fixed = TRUE)
})

test_that("the methods guide cites the method-choice references", {
  html <- as.character(methods_guide_ui())
  expect_match(html, "Conesa")
  expect_match(html, "Schurch")
  expect_match(html, "Rau, A.")
})

test_that("KEGG is marked as not accepting gene symbols", {
  html <- as.character(methods_guide_ui())
  row <- regmatches(html, regexpr("<td>Gene symbol</td>(\\s*<td>[a-z]+</td>){3}", html))
  expect_match(row, "<td>yes</td>\\s*<td>no</td>\\s*<td>yes</td>")
})
