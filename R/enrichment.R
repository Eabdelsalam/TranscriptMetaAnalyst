#' Organisms supported by PANTHER GO enrichment
#'
#' Queries the PANTHER web service for its list of organisms.
#'
#' @return A data frame of organisms, including `long_name` and `taxon_id`.
#' @export
#' @examples
#' \dontrun{
#' head(organisms_go()$long_name)
#' }
organisms_go <- function() {
  suppressMessages(rbioapi::rba_panther_info(what = "organisms"))
}

#' GO enrichment with PANTHER
#'
#' Runs PANTHER overrepresentation tests for the three GO aspects (biological
#' process, cellular component, molecular function) through `rbioapi`, and
#' returns the significant terms sorted by FDR. Needs an internet connection.
#'
#' @param genes Character vector of gene IDs (e.g. the meta-analysis DEGs).
#' @param organism NCBI taxon ID of the organism, e.g. `1111708` for
#'   *Synechocystis* sp. PCC 6803.
#' @param cutoff FDR cutoff passed to PANTHER.
#' @return A data frame of enriched GO terms with a `category` column
#'   (`"BP"`, `"CC"` or `"MF"`), sorted by FDR.
#' @export
#' @examples
#' \dontrun{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' head(run_go_enrichment(res$meta_degs, organism = 1111708))
#' }
run_go_enrichment <- function(genes, organism, cutoff = 0.05) {
  aspects <- c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
  parts <- lapply(names(aspects), function(category) {
    res <- suppressMessages(rbioapi::rba_panther_enrich(
      genes = genes, organism = as.numeric(organism),
      annot_dataset = aspects[[category]], cutoff = cutoff
    ))$result
    res$category <- category
    res
  })
  all_go <- do.call(rbind, parts)
  all_go <- all_go[!is.na(all_go$term.id), ] # unknown categories
  all_go <- all_go[all_go$fold_enrichment > 0, ] # categories with 0 genes
  all_go <- all_go[!all_go$term.id %in% aspects, ] # the three root terms
  all_go[order(all_go$fdr), ]
}

#' KEGG pathway enrichment
#'
#' Runs [clusterProfiler::enrichKEGG()], which queries the KEGG REST service.
#' Needs an internet connection.
#'
#' @param genes Character vector of gene IDs in a KEGG-supported ID type
#'   (e.g. locus tags such as `"sll1214"`).
#' @param organism KEGG organism code, e.g. `"syn"`.
#' @param cutoff q-value cutoff.
#' @return A data frame of enriched pathways.
#' @export
#' @examples
#' \dontrun{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' head(run_kegg_enrichment(res$meta_degs, organism = "syn"))
#' }
run_kegg_enrichment <- function(genes, organism, cutoff = 0.05) {
  res <- suppressMessages(clusterProfiler::enrichKEGG(gene = genes,
                                                     organism = organism,
                                                     qvalueCutoff = cutoff,
                                                     pAdjustMethod = "BH"))
  as.data.frame(res)
}
