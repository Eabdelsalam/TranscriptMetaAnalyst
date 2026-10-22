# S3. Microarray case study (Supplementary Section S3).
#
# The bundled microarray demo (wild-type Synechocystis sp. PCC 6803, high CO2
# vs 3 h and 24 h after the shift to low CO2; GSE82120, GSE63352, GSE68250;
# Fisher's method): genes kept after probe collapsing and intersection,
# meta-DEGs per comparison, new and side-tracked DEGs, recovery of known
# low-CO2-inducible genes (locus tags checked against KEGG), and the STRING and
# KEGG enrichment of the up-regulated DEGs. Needs an internet connection (KEGG,
# STRING).
#
#   Rscript dev/supplementary/02_microarray_case_study.R

source("dev/supplementary/00_setup.R")

arrays <- tma_microarray_demo_data()
analysis <- run_array_analysis(arrays$matrices, arrays$design, reference = "HC", meta_method = "fisher")
genes_kept <- attr(analysis$preparation$report, "common_genes")
prep <- analysis$preparation$report
save_table(prep, "S3_microarray_preparation.csv")

up_degs <- function(res) {
  fc <- mean_fold_change(res$fold_changes)
  res$meta_degs[fc[match(res$meta_degs, res$genes)] > 0]
}
summary_table <- do.call(rbind, lapply(names(analysis$results), function(con) {
  res <- analysis$results[[con]]
  d <- res$idd_irr
  up <- up_degs(res)
  data.frame(comparison = con, studies = length(res$studies), genes = length(res$genes),
             meta_degs = length(res$meta_degs), up = length(up), down = length(res$meta_degs) - length(up),
             new = unname(d[["IDD"]]), new_pct = unname(d[["IDR"]]),
             sidetracked = unname(d[["Loss"]]), sidetracked_pct = unname(d[["IRR"]]))
}))
save_table(summary_table, "S3_microarray_summary.csv")
print(summary_table)

# known low-CO2-inducible genes: expected symbols, verified against KEGG ------
kegg_genes <- local({
  lines <- readLines(url("https://rest.kegg.jp/list/syn"), warn = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  data.frame(locus = sub("^syn:", "", vapply(parts, `[`, "", 1)),
             description = vapply(parts, function(p) p[length(p)], ""))
})
kegg_symbol <- function(desc) ifelse(grepl(";", desc), trimws(sub(";.*", "", desc)), NA)
known <- data.frame(
  locus = c("slr1512", "slr1513", "slr0040", "slr0041", "slr0043", "slr0044",
            "sll1732", "sll1733", "sll1734", "sll0217", "sll0218", "sll0219"),
  expected = c("sbtA", "sbtB", "cmpA", "cmpB", "cmpC", "cmpD",
               "ndhF3", "ndhD3", "cupA", "flv4", "sll0218", "flv2"),
  system = c("Na+/HCO3- symporter SbtA", "SbtA-associated protein SbtB", "BCT1 HCO3- transporter",
             "BCT1 HCO3- transporter", "BCT1 HCO3- transporter", "BCT1 HCO3- transporter",
             "NDH-1(3) CO2 uptake complex", "NDH-1(3) CO2 uptake complex", "NDH-1(3) CO2 uptake complex",
             "flv4-2 operon", "flv4-2 operon", "flv4-2 operon"))
known$kegg <- kegg_genes$description[match(known$locus, kegg_genes$locus)]
known$kegg_symbol <- kegg_symbol(known$kegg)
# second source: UniProt gene and protein names (Synechocystis sp. PCC 6803
# substr. Kazusa, taxon 1111708). Some literature symbols (sbtA, cupA, flv4/2)
# are not used by KEGG or UniProt, which is why both names are reported.
uniprot <- local({
  query <- utils::URLencode(sprintf("(%s) AND (organism_id:1111708)",
                                    paste0("gene:", known$locus, collapse = " OR ")), reserved = TRUE)
  tsv <- utils::read.delim(url(paste0("https://rest.uniprot.org/uniprotkb/search?query=", query,
                                      "&fields=accession,gene_names,protein_name&format=tsv")),
                           check.names = FALSE)
  tsv$locus <- vapply(strsplit(tsv$`Gene Names`, " "), function(g) {
    hit <- g[g %in% known$locus]
    if (length(hit)) hit[1] else NA_character_
  }, "")
  tsv
})
# entries that do not list their locus tag as a gene name: free-text search
for (locus in setdiff(known$locus, uniprot$locus)) {
  query <- utils::URLencode(sprintf("%s AND (organism_id:1111708)", locus), reserved = TRUE)
  hit <- utils::read.delim(url(paste0("https://rest.uniprot.org/uniprotkb/search?query=", query,
                                      "&fields=accession,gene_names,protein_name&format=tsv")),
                           check.names = FALSE)
  if (nrow(hit) > 0) uniprot <- rbind(uniprot, cbind(hit[1, ], locus = locus))
}
i <- match(known$locus, uniprot$locus)
known$uniprot <- uniprot$Entry[i]
known$uniprot_genes <- uniprot$`Gene Names`[i]
known$uniprot_protein <- uniprot$`Protein names`[i]
known$verified_locus <- !is.na(known$kegg) & !is.na(known$uniprot)
in_comparison <- function(con) {
  res <- analysis$results[[con]]
  fc <- mean_fold_change(res$fold_changes)
  i <- match(known$locus, res$genes)
  data.frame(tested = !is.na(i), meta_deg = known$locus %in% res$meta_degs, mean_log2fc = round(fc[i], 2))
}
for (con in names(analysis$results)) {
  x <- in_comparison(con)
  tag <- gsub("[^A-Za-z0-9]+", "_", con)
  known[[paste0(tag, "_meta_deg")]] <- x$meta_deg
  known[[paste0(tag, "_log2fc")]] <- x$mean_log2fc
}
save_table(known, "S3_microarray_known_genes.csv")
print(known)

# other strongly induced genes at 24 h, with their KEGG description
res24 <- analysis$results[["LC_24h vs HC"]]
fc24 <- mean_fold_change(res24$fold_changes)
top_up <- utils::head(res24$meta_degs[order(-fc24[match(res24$meta_degs, res24$genes)])], 25)
top_up <- data.frame(locus = top_up, mean_log2fc = round(fc24[match(top_up, res24$genes)], 2),
                     kegg = kegg_genes$description[match(top_up, kegg_genes$locus)])
save_table(top_up, "S3_microarray_top_induced_24h.csv")

# enrichment of the up-regulated DEGs against all tested genes ------------------
enrich <- list()
for (con in names(analysis$results)) {
  res <- analysis$results[[con]]
  up <- up_degs(res)
  string <- tryCatch(run_string_enrichment(up, species = 1148, background = res$genes, id_type = "kegg"),
                     error = function(e) { message("STRING: ", conditionMessage(e)); NULL })
  kegg <- tryCatch(run_kegg_enrichment(up, organism = "syn", background = res$genes, id_type = "kegg"),
                   error = function(e) { message("KEGG: ", conditionMessage(e)); NULL })
  tag <- gsub("[^A-Za-z0-9]+", "_", con)
  if (!is.null(string)) write_enrichment(string, out(paste0("S3_", tag, "_up_STRING.txt")))
  if (!is.null(kegg)) write_enrichment(kegg, out(paste0("S3_", tag, "_up_KEGG.txt")))
  enrich[[con]] <- list(
    up = length(up),
    string_label = if (!is.null(string)) enrichment_label(string) else NA,
    string_top = if (!is.null(string)) utils::head(string$table[string$table$category != "PMID",
                                                                c("category", "description", "fdr")], 10),
    kegg_label = if (!is.null(kegg)) enrichment_label(kegg) else NA,
    kegg_top = if (!is.null(kegg)) utils::head(kegg$table[, c("Description", "qvalue")], 10))
}

saveRDS(list(genes_kept = genes_kept, preparation = prep, summary = summary_table, known = known,
             top_up = top_up, enrichment = enrich, session = session_line(),
             kegg_list_date = format(Sys.Date())),
        out("S3_microarray.rds"))
message("S3 done")
