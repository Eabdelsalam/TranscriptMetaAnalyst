# GEO is mocked: a small ExpressionSet stands in for a downloaded series
fake_eset <- function(two_colour = FALSE) {
  expr <- matrix(c(5, 6, 7, 8, 9, 10,
                   7, 7, 8, 8, 9, 9,
                   1, 2, 1, 2, 1, 2), nrow = 3, byrow = TRUE,
                 dimnames = list(c("p1", "p2", "p3"), paste0("GSM", 1:6)))
  pheno <- data.frame(title = paste0("sample", 1:6),
                      source_name_ch1 = "culture",
                      characteristics_ch1 = rep(c("treatment: control", "treatment: stress", "treatment: other"), each = 2),
                      label_ch1 = if (two_colour) rep(c("Cy3", "Cy5"), 3) else "biotin",
                      channel_count = "1", row.names = colnames(expr))
  feature <- data.frame(ID = rownames(expr), ORF = c("sll0001", "sll0002", ""),
                        Gene.symbol = c("abcA", "abcB", ""), row.names = rownames(expr))
  Biobase::ExpressionSet(expr, phenoData = Biobase::AnnotatedDataFrame(pheno),
                         featureData = Biobase::AnnotatedDataFrame(feature), annotation = "GPL123")
}

test_that("fetch_geo checks the accession and names the platforms", {
  seen <- NULL
  local_mocked_bindings(getGEO = function(GEO, ...) {
    seen <<- list(GEO = GEO, ...)
    list(fake_eset())
  }, .package = "GEOquery")
  sets <- fetch_geo(" gse123 ", destdir = withr::local_tempdir())
  expect_named(sets, "GPL123")
  expect_equal(seen$GEO, "GSE123")
  expect_true(seen$GSEMatrix)
  expect_error(fetch_geo("GDS99"), "not a GEO series accession", class = "tma_enrichment_error")

  local_mocked_bindings(getGEO = function(...) stop("HTTP error 404"), .package = "GEOquery")
  expect_error(fetch_geo("GSE999999999"), "could not be downloaded from GEO \\(HTTP error 404\\)",
               class = "tma_enrichment_error")
})

test_that("the GEO helpers suggest the condition and gene columns", {
  es <- fake_eset()
  expect_equal(geo_condition_columns(es)[1], "characteristics_ch1")
  expect_false("source_name_ch1" %in% geo_condition_columns(es)) # one value only
  expect_equal(geo_gene_columns(es)[1], "Gene.symbol")
  expect_false(geo_is_two_colour(es))
  expect_true(geo_is_two_colour(fake_eset(two_colour = TRUE)))
})

test_that("geo_study keeps the labelled samples and builds the study inputs", {
  study <- geo_study(fake_eset(), "GSE123", "characteristics_ch1",
                     labels = c("treatment: control" = "Control", "treatment: stress" = "Stress",
                                "treatment: other" = ""),
                     gene_column = "ORF")
  expect_equal(study$design, data.frame(sample = paste0("GSM", 1:4), study = "GSE123",
                                        condition = c("Control", "Control", "Stress", "Stress")))
  expect_named(study$matrix, c("ID", "Symbol", paste0("GSM", 1:4)))
  expect_equal(study$matrix$Symbol, c("sll0001", "sll0002", ""))

  # the result feeds straight into the microarray preparation
  prepared <- prepare_array_studies(list(GSE123 = study$matrix), study$design)
  expect_equal(rownames(prepared$studies$GSE123$values), c("sll0001", "sll0002"))

  expect_error(geo_study(fake_eset(), "GSE123", "characteristics_ch1", labels = c(x = "y")),
               "No sample of GSE123 has a condition label")
  expect_error(geo_study(fake_eset(), "GSE123", "nope", labels = c()), "no column 'nope'")
  no_gene <- geo_study(fake_eset(), "GSE123", "characteristics_ch1",
                       labels = c("treatment: control" = "C"), gene_column = NULL)
  expect_true(all(is.na(no_gene$matrix$Symbol)))
})
