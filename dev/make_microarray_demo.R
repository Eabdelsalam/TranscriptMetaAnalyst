# Builds the microarray demo data in inst/extdata/microarray_demo/ (D13).
#
# Wild-type Synechocystis sp. PCC 6803 grown at high CO2 (HC) and shifted to
# low CO2 (air) for 3 h (LC_3h) or 24 h (LC_24h), from three GEO series
# (platform GPL15867):
#   GSE82120: HC 6, LC_24h 6
#   GSE63352: HC 2, LC_3h 2, LC_24h 2
#   GSE68250: HC 2, LC_3h 2, LC_24h 2
# GSE31672 (GPL6546) has the same design but is a two-colour experiment whose
# channels are not comparable as single-channel samples, so it is not used.
# The series-matrix values (normalised, log2) are kept as published. To keep
# the package small, each file holds at most 3 probes per gene (locus tag) plus
# a few control probes without a gene, so the demo still shows probe
# collapsing and the removal of probes without a gene.
#
# Run from the package root (needs an internet connection):
#   Rscript dev/make_microarray_demo.R

suppressPackageStartupMessages(library(GEOquery))

out_dir <- file.path("inst", "extdata", "microarray_demo")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cache <- file.path(tempdir(), "geo")
dir.create(cache, showWarnings = FALSE)

# condition of each sample from its title and characteristics; NA = not used
classify <- function(p) {
  text <- tolower(apply(p[, grep("^(title|characteristics)", names(p))], 1, paste, collapse = " "))
  wild_type <- grepl("wt|wild ?type", text) & !grepl("mutant|delta|∆", text)
  condition <- ifelse(grepl("high c", text) | grepl("(^|[ _])hc([ _]|$)", text), "HC",
                      ifelse(grepl("3 ?h", text) & grepl("lc|low c", text), "LC_3h",
                             ifelse(grepl("24 ?h", text) & grepl("lc|low c", text), "LC_24h", NA)))
  ifelse(wild_type, condition, NA)
}

series <- c("GSE82120", "GSE63352", "GSE68250")
designs <- list()
for (gse in series) {
  es <- suppressMessages(getGEO(gse, GSEMatrix = TRUE, destdir = cache, getGPL = TRUE))[[1]]
  p <- Biobase::pData(es)
  condition <- classify(p)
  keep <- !is.na(condition)
  expr <- Biobase::exprs(es)[, keep, drop = FALSE]
  gene <- as.character(Biobase::fData(es)$ORF)
  is_gene <- grepl("^(s[lm][lr]|ss[lr]|sg[lr])[0-9]{4}", gene)
  # at most 3 probes per gene, in platform order, plus 20 probes without a gene
  rank_in_gene <- stats::ave(seq_along(gene), gene, FUN = seq_along)
  rows <- which((is_gene & rank_in_gene <= 3))
  rows <- sort(c(rows, utils::head(which(!is_gene), 20)))
  out <- data.frame(ID = rownames(expr)[rows], Symbol = ifelse(is_gene[rows], gene[rows], ""),
                    round(expr[rows, , drop = FALSE], 4), check.names = FALSE)
  utils::write.table(out, gzfile(file.path(out_dir, paste0(gse, ".txt.gz"))),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  designs[[gse]] <- data.frame(sample = colnames(expr), study = gse, condition = condition[keep])
  message(gse, ": ", ncol(expr), " samples (", paste(names(table(condition[keep])), table(condition[keep]), collapse = ", "),
          "), ", nrow(out), " probes written")
}
design <- do.call(rbind, unname(designs))
utils::write.table(design, file.path(out_dir, "design.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table(design$study, design$condition))
