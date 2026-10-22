# S2. RNA-Seq demo case study (Supplementary Section S2).
#
# The bundled RNA-Seq demo (3 studies x 4 samples, Treated vs Cont): meta-DEGs,
# new and side-tracked DEGs with DESeq2 + Fisher, and the same with edgeR and
# limma-voom and with the inverse-normal method, with their overlap.
#
#   Rscript dev/supplementary/03_rnaseq_demo_case_study.R

source("dev/supplementary/00_setup.R")

demo <- tma_demo_data()
contrast <- "Treated vs Cont"
de_methods <- c("DESeq2", "edgeR", "limma-voom")
metas <- c("Fisher" = "fisher", "Inverse normal" = "invnorm")

runs <- list()
rows <- list()
for (de in de_methods) {
  for (m in names(metas)) {
    res <- run_tma_analysis(demo$counts, demo$design, reference = "Cont", de_method = de,
                            meta_method = metas[[m]], batch_method = "none")$results[[contrast]]
    key <- paste(de, m, sep = " + ")
    runs[[key]] <- res
    d <- res$idd_irr
    rows[[key]] <- data.frame(de_method = de, meta_method = m, genes_tested = length(res$genes),
                              meta_degs = length(res$meta_degs),
                              new = unname(d[["IDD"]]), new_pct = unname(d[["IDR"]]),
                              sidetracked = unname(d[["Loss"]]), sidetracked_pct = unname(d[["IRR"]]),
                              study_degs = paste(lengths(res$study_degs), collapse = " / "))
  }
}
summary_table <- do.call(rbind, unname(rows))
save_table(summary_table, "S2_rnaseq_demo_methods.csv")
print(summary_table)

# overlap: share (%) of the row's meta-DEGs also found by the column's analysis
keys <- names(runs)
overlap <- outer(keys, keys, Vectorize(function(a, b) {
  da <- runs[[a]]$meta_degs
  100 * length(intersect(da, runs[[b]]$meta_degs)) / length(da)
}))
dimnames(overlap) <- list(keys, keys)
overlap_table <- data.frame(analysis = keys, round(overlap, 1), check.names = FALSE)
save_table(overlap_table, "S2_rnaseq_demo_overlap.csv")

reference <- runs[["DESeq2 + Fisher"]]$meta_degs
all_three <- Reduce(intersect, lapply(paste(de_methods, "Fisher", sep = " + "), function(k) runs[[k]]$meta_degs))
method_overlap <- 100 * length(intersect(reference, all_three)) / length(reference)
cat(sprintf("%.1f%% of the %d DESeq2 + Fisher meta-DEGs are also found with edgeR and limma-voom\n",
            method_overlap, length(reference)))

saveRDS(list(summary = summary_table, overlap = overlap_table, method_overlap = method_overlap,
             reference_degs = length(reference), shared_degs = length(intersect(reference, all_three)),
             session = session_line()),
        out("S2_rnaseq_demo.rds"))
message("S2 done")
