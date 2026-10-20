#' Gene ID types
#'
#' The ID types that the `Geneid` column of the counts file can contain. The
#' enrichment services map them to their own identifiers: PANTHER and STRING
#' accept all of them; KEGG accepts locus tags / KEGG gene IDs, NCBI Entrez
#' Gene IDs and UniProt accessions.
#'
#' @return A named character vector of ID-type codes, named by their labels.
#' @export
#' @examples
#' id_types()
id_types <- function() {
  c("Gene symbol" = "symbol",
    "Locus tag / KEGG gene ID" = "kegg",
    "NCBI Entrez Gene ID" = "entrez",
    "UniProt accession" = "uniprot",
    "Ensembl gene ID" = "ensembl")
}

# KEGG keyType of an ID type, NA when KEGG does not accept it
kegg_key_type <- function(id_type) {
  unname(c(kegg = "kegg", entrez = "ncbi-geneid", uniprot = "uniprot",
           symbol = NA, ensembl = NA)[id_type])
}

# An error with a message written for the user
enrichment_error <- function(message) {
  structure(class = c("tma_enrichment_error", "error", "condition"),
            list(message = message, call = NULL))
}

#' Readable message for a failed enrichment
#'
#' @param error The condition raised by an enrichment function.
#' @param service Name of the service, e.g. `"PANTHER"`.
#' @return A message for the user: the error's own message for problems TMA
#'   explains itself (unsupported ID type, unknown organism, nothing mapped),
#'   otherwise a note that the service could not be reached.
#' @export
#' @examples
#' enrichment_message(simpleError("Timeout was reached"), "STRING")
enrichment_message <- function(error, service) {
  if (inherits(error, "tma_enrichment_error")) return(conditionMessage(error))
  sprintf("The %s enrichment could not be completed (%s). This is usually a network problem or a temporary problem of the %s service; please check your internet connection and try again later.",
          service, conditionMessage(error), service)
}

#' ID-mapping report
#'
#' Summarises how many of the submitted IDs an enrichment service recognised.
#' A warning is given when fewer than half of them were recognised, unless
#' the background (all tested genes) is recognised at a similar rate: for
#' example, KEGG only counts genes that are in a pathway, which is often a
#' minority of all genes, so a low rate is then expected.
#'
#' @param submitted Character vector of submitted IDs.
#' @param mapped Character vector of the submitted IDs that were recognised.
#' @param service Name of the service.
#' @param background_rate Fraction of the background that the service
#'   recognised, or `NA` if unknown.
#' @param what How the service recognises IDs, used in the warning.
#' @return A list with `submitted`, `mapped` and `unmapped` (counts),
#'   `unmapped_ids`, `rate` (fraction mapped), `background_rate` and
#'   `warning` (a message, or `NULL`).
#' @export
#' @examples
#' mapping_report(c("a", "b", "c"), "a", "KEGG")$warning
#' mapping_report(c("a", "b", "c"), "a", "KEGG", background_rate = 0.3)$warning
mapping_report <- function(submitted, mapped, service, background_rate = NA,
                           what = paste("recognised by", service)) {
  submitted <- unique(submitted)
  mapped <- intersect(submitted, mapped)
  rate <- if (length(submitted) > 0) length(mapped) / length(submitted) else 0
  similar_to_background <- !is.na(background_rate) && rate >= 0.8 * background_rate
  warning <- if (rate < 0.5 && !similar_to_background) {
    context <- if (is.na(background_rate)) "" else
      sprintf(" (%.0f%% of all tested genes are)", 100 * background_rate)
    sprintf("Only %.0f%% of the IDs were %s%s. Check that the selected ID type matches the counts file, or try another service (STRING recognises most ID types and organisms).",
            100 * rate, what, context)
  }
  list(submitted = length(submitted), mapped = length(mapped),
       unmapped = length(submitted) - length(mapped),
       unmapped_ids = setdiff(submitted, mapped), rate = rate,
       background_rate = background_rate, warning = warning)
}

# The result of one enrichment run
new_enrichment <- function(service, table, mapping, version, background, organism, id_type) {
  structure(list(service = service, table = table, mapping = mapping, version = version,
                 background = background, organism = organism, id_type = id_type),
            class = "tma_enrichment")
}

#' @export
print.tma_enrichment <- function(x, ...) {
  cat(enrichment_label(x), "\n")
  cat(" ", x$mapping$mapped, "of", x$mapping$submitted, "IDs mapped;",
      nrow(x$table), "enriched terms\n")
  invisible(x)
}

#' Label of an enrichment result
#'
#' The service, its database version and the background used, for reporting
#' the result.
#'
#' @param enrichment A `tma_enrichment` object.
#' @return A character string.
#' @export
#' @examples
#' \dontrun{
#' enrichment_label(run_string_enrichment(genes, species = 1148))
#' }
enrichment_label <- function(enrichment) {
  sprintf("%s (%s); organism: %s; background: %s.", enrichment$service, enrichment$version,
          enrichment$organism, enrichment$background)
}

# "all N tested genes (M recognised by <service>)" or the service default
background_label <- function(background, recognised, service, default) {
  if (is.null(background)) return(default)
  sprintf("all %d tested genes (%d recognised by %s)", length(unique(background)), recognised, service)
}

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

#' Organisms supported by STRING
#'
#' Downloads the species list of the current STRING version (about 12,500
#' organisms).
#'
#' @return A data frame with `taxon_id`, `name` and `domain`.
#' @export
#' @examples
#' \dontrun{
#' head(organisms_string())
#' }
organisms_string <- function() {
  version <- suppressMessages(rbioapi::rba_string_version())$string_version
  url <- sprintf("https://stringdb-downloads.org/download/species.v%s.txt", version)
  species <- utils::read.delim(url, quote = "", comment.char = "", check.names = FALSE)
  data.frame(taxon_id = species[[1]], name = species[[4]], domain = species[[5]])
}

#' Check a KEGG organism code
#'
#' Looks the code up in the KEGG organism table bundled with clusterProfiler,
#' without a network call.
#'
#' @param code A KEGG organism code such as `"syn"` or `"hsa"`.
#' @return The organism's scientific name, or `NA` if the code is unknown.
#' @export
#' @examples
#' kegg_organism_name("syn")
#' kegg_organism_name("xyz123")
kegg_organism_name <- function(code) {
  if (length(code) != 1 || is.na(code) || !nzchar(code)) return(NA_character_)
  hits <- clusterProfiler::search_kegg_organism(paste0("^", code, "$"), by = "kegg_code")
  if (nrow(hits) == 0) NA_character_ else hits$scientific_name[1]
}

# Date of the KEGG pathway database, e.g. "pathway database of 2026/10/15";
# NA if KEGG cannot be reached
kegg_release <- function() {
  info <- tryCatch(readLines("https://rest.kegg.jp/info/kegg", warn = FALSE),
                   error = function(e) character(0))
  pathway <- grep("^\\s*pathway\\s", info, value = TRUE)
  if (length(pathway) == 0) return(NA_character_)
  paste("pathway database of", utils::tail(strsplit(trimws(pathway[1]), "\\s+")[[1]], 1))
}

#' GO enrichment with PANTHER
#'
#' Runs PANTHER overrepresentation tests (Fisher's exact test, FDR) for the
#' three GO aspects (biological process, cellular component, molecular
#' function) through `rbioapi`. PANTHER maps the submitted IDs (gene symbols,
#' UniProt, Ensembl, NCBI and organism-specific IDs such as locus tags) on its
#' server. Needs an internet connection.
#'
#' @param genes Character vector of gene IDs (e.g. the meta-analysis DEGs).
#' @param organism NCBI taxon ID of the organism, e.g. `1111708` for
#'   *Synechocystis* sp. PCC 6803.
#' @param background The reference list: all genes tested in the analysis.
#'   `NULL` uses PANTHER's whole-genome default.
#' @param id_type ID type of `genes`, one of [id_types()] (for the record;
#'   PANTHER recognises the type itself).
#' @param cutoff FDR cutoff passed to PANTHER.
#' @return A `tma_enrichment` object: a list with `service`, `table` (enriched
#'   GO terms with a `category` column, `"BP"`, `"CC"` or `"MF"`, sorted by
#'   FDR), `mapping` (see [mapping_report()]), `version`, `background`,
#'   `organism` and `id_type`.
#' @export
#' @examples
#' \dontrun{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' go <- run_go_enrichment(res$meta_degs, organism = 1111708, background = res$genes)
#' go$mapping$rate
#' }
run_go_enrichment <- function(genes, organism, background = NULL, id_type = "kegg", cutoff = 0.05) {
  aspects <- c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
  responses <- lapply(aspects, function(dataset) {
    args <- list(genes = genes, organism = as.numeric(organism), annot_dataset = dataset,
                 cutoff = cutoff)
    if (!is.null(background)) {
      args$ref_genes <- background
      args$ref_organism <- as.numeric(organism)
    }
    suppressMessages(do.call(rbioapi::rba_panther_enrich, args))
  })
  first <- responses[[1]]
  if (is.null(first$result) && is.null(first$input_list)) {
    stop("PANTHER returned no result")
  }
  ids <- function(x) if (is.null(x) || all(is.na(x))) character(0) else trimws(unlist(strsplit(paste(x, collapse = ","), ",")))
  ref_mapped <- first$reference$mapped_count
  ref_mapped <- if (is.list(ref_mapped)) ref_mapped$content else ref_mapped
  bg_rate <- if (!is.null(background) && !is.null(ref_mapped)) ref_mapped / length(unique(background)) else NA
  mapping <- mapping_report(genes, ids(first$input_list$mapped_ids), "PANTHER", background_rate = bg_rate)
  if (mapping$mapped == 0) {
    stop(enrichment_error("PANTHER did not recognise any of the submitted IDs for this organism. Check the organism and the ID type, or try STRING."))
  }
  parts <- lapply(names(aspects), function(category) {
    res <- responses[[category]]$result
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res$category <- category
    res
  })
  all_go <- do.call(rbind, parts)
  if (is.null(all_go)) all_go <- data.frame(term.id = character(0), term.label = character(0),
                                            fold_enrichment = numeric(0), fdr = numeric(0),
                                            number_in_list = integer(0), category = character(0))
  all_go <- all_go[!is.na(all_go$term.id), ] # unknown categories
  all_go <- all_go[all_go$fold_enrichment > 0, ] # categories with 0 genes
  all_go <- all_go[!all_go$term.id %in% aspects, ] # the three root terms
  all_go <- all_go[order(all_go$fdr), ]
  rownames(all_go) <- NULL

  new_enrichment(
    service = "PANTHER GO enrichment",
    table = all_go,
    mapping = mapping,
    version = sprintf("PANTHER release %s; %s", first$tool_release_date,
                      gsub("\\s+", " ", first$annot_version_release_date)),
    background = background_label(background, if (is.null(ref_mapped)) NA else ref_mapped,
                                  "PANTHER", "whole genome (PANTHER default)"),
    organism = as.character(organism),
    id_type = id_type
  )
}

#' KEGG pathway enrichment
#'
#' Runs [clusterProfiler::enrichKEGG()], which queries the KEGG REST service.
#' KEGG accepts locus tags / KEGG gene IDs, NCBI Entrez Gene IDs and UniProt
#' accessions; for gene symbols or Ensembl IDs use [run_string_enrichment()]
#' or [run_go_enrichment()] instead. Needs an internet connection.
#'
#' @inheritParams run_go_enrichment
#' @param organism KEGG organism code, e.g. `"syn"`; checked with
#'   [kegg_organism_name()] before the query.
#' @param background The universe: all genes tested in the analysis. `NULL`
#'   uses all genes of the organism in KEGG.
#' @param cutoff q-value cutoff.
#' @return A `tma_enrichment` object (see [run_go_enrichment()]) with the
#'   enriched pathways in `table`.
#' @export
#' @examples
#' \dontrun{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' kegg <- run_kegg_enrichment(res$meta_degs, organism = "syn", background = res$genes)
#' }
run_kegg_enrichment <- function(genes, organism, background = NULL, id_type = "kegg", cutoff = 0.05) {
  key_type <- kegg_key_type(id_type)
  if (is.na(key_type)) {
    label <- c(symbol = "gene symbols", ensembl = "Ensembl gene IDs")[[id_type]]
    stop(enrichment_error(sprintf("KEGG does not accept %s. Use the STRING tab (which maps most ID types) or the Gene Ontology tab instead.",
                                  label)))
  }
  organism_name <- kegg_organism_name(organism)
  if (is.na(organism_name)) {
    stop(enrichment_error(sprintf("'%s' is not a KEGG organism code. Look up the code (e.g. 'hsa' for human, 'syn' for Synechocystis sp. PCC 6803) in the KEGG organism list.",
                                  organism)))
  }
  res <- suppressMessages(clusterProfiler::enrichKEGG(gene = genes, organism = organism,
                                                      keyType = key_type, universe = background,
                                                      qvalueCutoff = cutoff, pAdjustMethod = "BH"))
  if (is.null(res)) {
    stop(enrichment_error(sprintf("None of the submitted IDs was found in the KEGG pathways of %s. Check the organism code and the ID type, or try STRING.",
                                  organism_name)))
  }
  in_pathways <- unique(unlist(res@geneSets))
  mapped <- genes[genes %in% in_pathways]
  table <- as.data.frame(res)
  rownames(table) <- NULL
  release <- kegg_release()
  bg_rate <- if (!is.null(background)) length(res@universe) / length(unique(background)) else NA
  new_enrichment(
    service = "KEGG pathway enrichment (clusterProfiler)",
    table = table,
    mapping = mapping_report(genes, mapped, "KEGG", background_rate = bg_rate,
                             what = "in a KEGG pathway"),
    version = if (is.na(release)) "KEGG, release date not available" else paste("KEGG", release),
    background = background_label(background, length(res@universe), "KEGG",
                                  "all genes of the organism in KEGG"),
    organism = sprintf("%s (%s)", organism_name, organism),
    id_type = id_type
  )
}

#' Functional enrichment with STRING
#'
#' Maps the submitted IDs to STRING identifiers ([rbioapi::rba_string_map_ids()];
#' STRING recognises gene symbols, locus tags, UniProt and most other IDs for
#' more than 12,000 organisms) and runs STRING's functional enrichment
#' ([rbioapi::rba_string_enrichment()]), which covers GO, KEGG, Reactome, Pfam,
#' InterPro and other annotation sources. The background (all tested genes) is
#' mapped to STRING identifiers as well. Needs an internet connection.
#'
#' @inheritParams run_go_enrichment
#' @param species NCBI taxon ID of the organism as used by STRING, e.g. `1148`
#'   for *Synechocystis* sp. PCC 6803 (see [organisms_string()]).
#' @param background All genes tested in the analysis. `NULL` uses STRING's
#'   default background (the whole genome).
#' @param cutoff FDR cutoff.
#' @return A `tma_enrichment` object (see [run_go_enrichment()]) with the
#'   enriched terms in `table`: `category`, `term`, `description`,
#'   `number_of_genes`, `number_of_genes_in_background`, `p_value`, `fdr` and
#'   `inputGenes`, sorted by FDR.
#' @export
#' @examples
#' \dontrun{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' string <- run_string_enrichment(res$meta_degs, species = 1148, background = res$genes)
#' table(string$table$category)
#' }
run_string_enrichment <- function(genes, species, background = NULL, id_type = "kegg", cutoff = 0.05) {
  species <- as.numeric(species)
  string_ids <- function(ids) {
    mapped <- suppressMessages(rbioapi::rba_string_map_ids(ids, species = species))
    check_string_response(mapped)
    if (is.null(mapped) || nrow(mapped) == 0) return(data.frame(queryItem = character(0), stringId = character(0)))
    mapped
  }
  input <- string_ids(genes)
  if (length(intersect(genes, input$queryItem)) == 0) {
    stop(enrichment_error("STRING did not recognise any of the submitted IDs for this organism. Check the organism and the ID type."))
  }
  bg <- if (!is.null(background)) string_ids(background)
  bg_rate <- if (!is.null(bg)) length(intersect(background, bg$queryItem)) / length(unique(background)) else NA
  mapping <- mapping_report(genes, input$queryItem, "STRING", background_rate = bg_rate)
  args <- list(ids = unique(input$stringId), species = species)
  if (!is.null(bg)) args$background <- unique(bg$stringId)
  res <- suppressMessages(do.call(rbioapi::rba_string_enrichment, args))
  check_string_response(res)
  table <- if (is.data.frame(res)) res else do.call(rbind, unname(res))
  cols <- c("category", "term", "description", "number_of_genes", "number_of_genes_in_background",
            "p_value", "fdr", "inputGenes", "preferredNames")
  if (is.null(table) || nrow(table) == 0) {
    table <- stats::setNames(data.frame(matrix(nrow = 0, ncol = length(cols))), cols)
  } else {
    table <- table[table$fdr <= cutoff, intersect(cols, names(table)), drop = FALSE]
    table <- table[order(table$fdr), , drop = FALSE]
  }
  rownames(table) <- NULL
  version <- suppressMessages(rbioapi::rba_string_version())
  new_enrichment(
    service = "STRING functional enrichment",
    table = table,
    mapping = mapping,
    version = paste("STRING version", version$string_version),
    background = background_label(background, if (is.null(bg)) NA else length(unique(bg$queryItem)),
                                  "STRING", "whole genome (STRING default)"),
    organism = as.character(species),
    id_type = id_type
  )
}

# STRING reports problems as a data frame with an "error" column
check_string_response <- function(x) {
  if (is.data.frame(x) && "error" %in% names(x)) {
    stop(sprintf("STRING returned an error: %s", paste(x$message, collapse = " ")))
  }
  invisible(x)
}
