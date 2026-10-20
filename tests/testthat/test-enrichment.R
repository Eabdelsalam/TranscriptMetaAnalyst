# The PANTHER and KEGG web services are mocked so the tests run offline

fake_panther <- function(genes, organism, annot_dataset, cutoff) {
  list(result = data.frame(
    term.id = c(annot_dataset, paste0(annot_dataset, "1"), NA, paste0(annot_dataset, "2")),
    term.label = c("root", "term one", "unclassified", "term two"),
    number_in_list = c(10, 5, 1, 3),
    fold_enrichment = c(1, 2.5, 1, 0),
    fdr = c(0.001, 0.02, 0.5, 0.01),
    plus_minus = "+"
  ))
}

test_that("run_go_enrichment queries the three GO aspects and cleans the result", {
  calls <- list()
  local_mocked_bindings(
    rba_panther_enrich = function(...) {
      args <- list(...)
      calls[[length(calls) + 1]] <<- args
      do.call(fake_panther, args)
    },
    .package = "rbioapi"
  )
  go <- run_go_enrichment(c("g1", "g2"), organism = "1111708")

  expect_length(calls, 3)
  expect_equal(vapply(calls, `[[`, "", "annot_dataset"), c("GO:0008150", "GO:0005575", "GO:0003674"))
  expect_equal(calls[[1]]$genes, c("g1", "g2"))
  expect_equal(calls[[1]]$organism, 1111708)
  # roots, unclassified and zero-enrichment terms removed; sorted by FDR
  expect_equal(nrow(go), 3)
  expect_false(any(go$term.id %in% c("GO:0008150", "GO:0005575", "GO:0003674")))
  expect_setequal(go$category, c("BP", "CC", "MF"))
  expect_false(is.unsorted(go$fdr))
  expect_s3_class(plot_go_lollipop(go), "ggplot")
})

test_that("run_kegg_enrichment passes the genes and organism to enrichKEGG", {
  seen <- NULL
  local_mocked_bindings(
    enrichKEGG = function(gene, organism, qvalueCutoff, pAdjustMethod, ...) {
      seen <<- list(gene = gene, organism = organism, q = qvalueCutoff, adj = pAdjustMethod)
      data.frame(Description = c("Photosynthesis", "Ribosome"),
                 subcategory = c("Energy metabolism", "Translation"),
                 FoldEnrichment = c(3, 2), qvalue = c(0.001, 0.03), Count = c(12, 8))
    },
    .package = "clusterProfiler"
  )
  kegg <- run_kegg_enrichment(c("sll1214", "slr0612"), organism = "syn")

  expect_equal(seen, list(gene = c("sll1214", "slr0612"), organism = "syn", q = 0.05, adj = "BH"))
  expect_equal(nrow(kegg), 2)
  expect_no_error(ggplot2::ggplot_build(plot_kegg_lollipop(kegg)))
})
