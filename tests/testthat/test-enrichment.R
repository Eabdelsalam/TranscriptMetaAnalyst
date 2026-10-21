# The PANTHER, KEGG and STRING web services are mocked so the tests run offline

fake_panther <- function(genes, organism, annot_dataset, cutoff, ref_genes = NULL, ref_organism = NULL, ...) {
  mapped <- setdiff(genes, c("unknown1", "unknown2"))
  list(
    result = data.frame(
      term.id = c(annot_dataset, paste0(annot_dataset, "1"), NA, paste0(annot_dataset, "2")),
      term.label = c("root", "term one", "unclassified", "term two"),
      number_in_list = c(10, 5, 1, 3),
      fold_enrichment = c(1, 2.5, 1, 0),
      fdr = c(0.001, 0.02, 0.5, 0.01),
      plus_minus = "+"
    ),
    input_list = list(organism = "Synechocystis", mapped_ids = paste(mapped, collapse = ","),
                      unmapped_ids = paste(setdiff(genes, mapped), collapse = ",")),
    reference = list(mapped_count = list(content = if (is.null(ref_genes)) NULL else length(ref_genes) - 1)),
    tool_release_date = 20240807,
    annot_version_release_date = "GO Ontology database DOI:   Released 2026-04-28"
  )
}

test_that("run_go_enrichment uses the tested genes as reference and reports the mapping", {
  calls <- list()
  local_mocked_bindings(
    rba_panther_enrich = function(...) {
      args <- list(...)
      calls[[length(calls) + 1]] <<- args
      do.call(fake_panther, args)
    },
    .package = "rbioapi"
  )
  genes <- c("g1", "g2", "g3", "unknown1")
  background <- paste0("g", 1:20)
  go <- run_go_enrichment(genes, organism = "1111708", background = background)

  expect_s3_class(go, "tma_enrichment")
  expect_length(calls, 3)
  expect_equal(vapply(calls, `[[`, "", "annot_dataset"), c("GO:0008150", "GO:0005575", "GO:0003674"))
  expect_equal(calls[[1]]$ref_genes, background)
  expect_equal(calls[[1]]$ref_organism, 1111708)
  # the mapping report and labels
  expect_equal(go$mapping$submitted, 4)
  expect_equal(go$mapping$mapped, 3)
  expect_equal(go$mapping$unmapped_ids, "unknown1")
  expect_null(go$mapping$warning)
  expect_match(go$version, "PANTHER release 20240807; GO Ontology database DOI: Released 2026-04-28")
  expect_equal(go$background, "all 20 tested genes (19 recognised by PANTHER)")
  expect_equal(go$organism, "Synechocystis (1111708)")
  # roots, unclassified and zero-enrichment terms removed; sorted by FDR
  expect_equal(nrow(go$table), 3)
  expect_false(any(go$table$term.id %in% c("GO:0008150", "GO:0005575", "GO:0003674")))
  expect_setequal(go$table$category, c("BP", "CC", "MF"))
  expect_false(is.unsorted(go$table$fdr))
  expect_s3_class(plot_go_lollipop(go$table), "ggplot")
  expect_output(print(go), "3 of 4 IDs mapped")

  # without a background, PANTHER's whole-genome default is used
  calls <- list()
  default <- run_go_enrichment(genes, organism = 1111708)
  expect_null(calls[[1]]$ref_genes)
  expect_equal(default$background, "whole genome (PANTHER default)")
})

test_that("run_go_enrichment explains when no ID is recognised", {
  local_mocked_bindings(
    rba_panther_enrich = function(genes, ...) fake_panther(rep("unknown1", length(genes)), ...),
    .package = "rbioapi"
  )
  expect_error(run_go_enrichment(c("unknown1", "unknown2"), organism = 1), class = "tma_enrichment_error")
})

fake_kegg_result <- function(genes, universe) {
  methods::new("enrichResult",
               result = data.frame(ID = c("syn00195", "syn03010"),
                                   Description = c("Photosynthesis", "Ribosome"),
                                   subcategory = c("Energy metabolism", "Translation"),
                                   FoldEnrichment = c(3, 2), pvalue = c(1e-4, 0.01),
                                   p.adjust = c(0.001, 0.03), qvalue = c(0.001, 0.03), Count = c(12, 8)),
               gene = genes, universe = universe,
               geneSets = list(syn00195 = c("sll1214", "slr0612"), syn03010 = c("sll1399")),
               keytype = "kegg", organism = "syn", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
               pAdjustMethod = "BH", ontology = "KEGG", readable = FALSE)
}

test_that("run_kegg_enrichment passes the ID type, universe and organism to enrichKEGG", {
  seen <- NULL
  local_mocked_bindings(
    enrichKEGG = function(gene, organism, keyType, universe, qvalueCutoff, pAdjustMethod, ...) {
      seen <<- list(gene = gene, organism = organism, keyType = keyType, universe = universe,
                    q = qvalueCutoff, adj = pAdjustMethod)
      fake_kegg_result(gene, if (is.null(universe)) character(0) else universe[1:3])
    },
    .package = "clusterProfiler"
  )
  local_mocked_bindings(kegg_release = function() "pathway database of 2026/10/15")
  genes <- c("sll1214", "slr0612", "sll9999")
  background <- c(genes, "sll1399", "slr0001", "slr0002")
  kegg <- run_kegg_enrichment(genes, organism = "syn", background = background)

  expect_equal(seen$keyType, "kegg")
  expect_equal(seen$universe, background)
  expect_equal(seen$organism, "syn")
  expect_equal(c(seen$q, seen$adj), c(0.05, "BH"))
  expect_equal(kegg$mapping$mapped, 2)
  expect_equal(kegg$mapping$unmapped_ids, "sll9999")
  expect_equal(kegg$organism, "Synechocystis sp. PCC 6803 (syn)")
  expect_equal(kegg$version, "KEGG pathway database of 2026/10/15")
  expect_equal(kegg$background, "all 6 tested genes (3 recognised by KEGG)")
  expect_equal(nrow(kegg$table), 2)
  expect_no_error(ggplot2::ggplot_build(plot_kegg_lollipop(kegg$table)))

  run_kegg_enrichment(c("1234", "5678"), organism = "syn", id_type = "entrez")
  expect_equal(seen$keyType, "ncbi-geneid")
  run_kegg_enrichment(c("P12345"), organism = "syn", id_type = "uniprot")
  expect_equal(seen$keyType, "uniprot")
})

test_that("run_kegg_enrichment refuses unsupported ID types and unknown organisms before querying", {
  queried <- FALSE
  local_mocked_bindings(enrichKEGG = function(...) {
    queried <<- TRUE
    NULL
  }, .package = "clusterProfiler")
  expect_error(run_kegg_enrichment("TP53", organism = "hsa", id_type = "symbol"),
               "KEGG does not accept gene symbols", class = "tma_enrichment_error")
  expect_error(run_kegg_enrichment("ENSG0001", organism = "hsa", id_type = "ensembl"),
               "KEGG does not accept Ensembl gene IDs", class = "tma_enrichment_error")
  expect_error(run_kegg_enrichment("sll1214", organism = "xyz123"),
               "'xyz123' is not a KEGG organism code", class = "tma_enrichment_error")
  expect_false(queried)
  # nothing mapped: enrichKEGG returns NULL
  expect_error(run_kegg_enrichment("sll1214", organism = "syn"), "None of the submitted IDs",
               class = "tma_enrichment_error")
  expect_true(queried)
})

test_that("kegg_organism_name looks codes up offline", {
  expect_equal(kegg_organism_name("syn"), "Synechocystis sp. PCC 6803")
  expect_equal(kegg_organism_name("hsa"), "Homo sapiens")
  expect_true(is.na(kegg_organism_name("xyz123")))
  expect_true(is.na(kegg_organism_name("")))
})

fake_string_map <- function(ids, species, ...) {
  known <- ids[!grepl("^unknown", ids)]
  data.frame(queryItem = known, stringId = paste0(species, ".", known),
             taxonName = rep("Synechocystis sp. PCC 6803", length(known)))
}

test_that("run_string_enrichment maps the input and the background to STRING identifiers", {
  enrich_args <- NULL
  local_mocked_bindings(
    rba_string_map_ids = fake_string_map,
    rba_string_enrichment = function(ids, species, background = NULL, ...) {
      enrich_args <<- list(ids = ids, species = species, background = background)
      list(
        Process = data.frame(category = "Process", term = "GO:0015979", description = "Photosynthesis",
                             number_of_genes = 12, number_of_genes_in_background = 40,
                             p_value = 1e-6, fdr = 1e-4, inputGenes = "a,b"),
        PMID = data.frame(category = "PMID", term = "PMID:1", description = "A paper",
                          number_of_genes = 5, number_of_genes_in_background = 9,
                          p_value = 1e-3, fdr = 0.01, inputGenes = "a"),
        KEGG = data.frame(category = "KEGG", term = "syn00195", description = "Photosynthesis",
                          number_of_genes = 3, number_of_genes_in_background = 30,
                          p_value = 0.01, fdr = 0.2, inputGenes = "a")
      )
    },
    rba_string_version = function(...) list(string_version = "12.0"),
    .package = "rbioapi"
  )
  genes <- c("g1", "g2", "unknown1")
  background <- c(genes, "g3", "unknown2")
  string <- run_string_enrichment(genes, species = "1148", background = background)

  expect_equal(enrich_args$ids, c("1148.g1", "1148.g2"))
  expect_equal(enrich_args$species, 1148)
  expect_equal(enrich_args$background, c("1148.g1", "1148.g2", "1148.g3"))
  expect_equal(string$mapping$mapped, 2)
  expect_equal(string$mapping$unmapped_ids, "unknown1")
  expect_equal(string$version, "STRING version 12.0")
  expect_equal(string$organism, "Synechocystis sp. PCC 6803 (1148)")
  expect_equal(string$background, "all 5 tested genes (3 recognised by STRING)")
  # FDR cutoff applied, sorted by FDR
  expect_equal(string$table$term, c("GO:0015979", "PMID:1"))
  p <- plot_string_lollipop(string$table)
  expect_equal(p$data$term, "GO:0015979") # publications are not plotted
})

test_that("STRING errors and failed connections give readable messages", {
  local_mocked_bindings(
    rba_string_map_ids = function(...) data.frame(error = "unknown_error", message = "Something went wrong!"),
    .package = "rbioapi"
  )
  err <- tryCatch(run_string_enrichment("g1", species = 1148), error = function(e) e)
  expect_match(conditionMessage(err), "STRING could not run the request: Something went wrong!")
  msg <- enrichment_message(err, "STRING")
  expect_match(msg, "The STRING enrichment could not be completed")
  expect_match(msg, "check your internet connection")

  local_mocked_bindings(
    rba_string_map_ids = function(...) stop("Timeout was reached: Operation timed out after 30001 milliseconds"),
    .package = "rbioapi"
  )
  err <- tryCatch(run_string_enrichment("g1", species = 1148), error = function(e) e)
  expect_match(enrichment_message(err, "STRING"), "Timeout was reached")
  # TMA's own explanations are shown as they are
  own <- enrichment_error("KEGG does not accept gene symbols.")
  expect_equal(enrichment_message(own, "KEGG"), "KEGG does not accept gene symbols.")
})

test_that("mapping_report warns about low mapping rates unless the background is similar", {
  low <- mapping_report(letters[1:10], letters[1:3], "KEGG")
  expect_equal(low$rate, 0.3)
  expect_match(low$warning, "Only 30% of the IDs were recognised by KEGG")
  similar <- mapping_report(letters[1:10], letters[1:3], "KEGG", background_rate = 0.32,
                            what = "in a KEGG pathway")
  expect_null(similar$warning)
  worse <- mapping_report(letters[1:10], "a", "KEGG", background_rate = 0.32, what = "in a KEGG pathway")
  expect_match(worse$warning, "Only 10% of the IDs were in a KEGG pathway \\(32% of all tested genes are\\)")
  expect_null(mapping_report(letters[1:4], letters[1:3], "STRING")$warning)
})

test_that("id_types lists the ID types and KEGG's key types", {
  expect_setequal(id_types(), c("symbol", "kegg", "entrez", "uniprot", "ensembl"))
  expect_equal(kegg_key_type(c("kegg", "entrez", "uniprot")), c("kegg", "ncbi-geneid", "uniprot"))
  expect_true(all(is.na(kegg_key_type(c("symbol", "ensembl")))))
})
