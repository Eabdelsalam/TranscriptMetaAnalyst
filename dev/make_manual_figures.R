# Regenerates the result figures of manual/TMA_Manual.md from the demo data,
# drawn the same way as the app (inst/TMA/server.R). Run from the package root:
#
#   Rscript -e 'devtools::load_all(); source("dev/make_manual_figures.R")'
#
# The GO figure queries the PANTHER web service and needs an internet connection.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

testthat::source_test_helpers("tests/testthat", env = environment())
demo <- demo_pipeline()
fig_dir <- file.path("manual", "figures")
design <- demo$design
rownames(design) <- design$sample

# PCA (Figure 10)
full_counts <- combine_study_matrices(demo$norm_counts, design$sample)
se <- SummarizedExperiment(log2(full_counts + 1), colData = design)
pca <- suppressMessages(plotPCA(DESeqTransform(se)))
ggsave(file.path(fig_dir, "fig10_pca.png"), pca, width = 15, height = 12, units = "cm", dpi = 150)

# Heatmap with the default 2 gene clusters (Figure 11)
heat <- scale_deg_counts(demo$norm_counts, demo$meta_degs, design$sample)
png(file.path(fig_dir, "fig11_heatmap.png"), width = 7, height = 10, units = "in", res = 120)
pheatmap::pheatmap(heat,
                   cutree_rows = heatmap_row_clusters(nrow(heat)),
                   cutree_cols = length(demo$studies),
                   show_rownames = FALSE,
                   color = colorRampPalette(c("yellow", "black", "purple"))(150),
                   annotation_col = design[, c("condition", "study")])
dev.off()

# Volcano plot (Figure 13)
for_volcano <- data.frame(genes = demo$genenames,
                          fc = mean_fold_change(demo$fold_changes),
                          pval = demo$fishcomb$adjpval)
volcano <- EnhancedVolcano::EnhancedVolcano(for_volcano, lab = for_volcano$genes,
                                            x = "fc", y = "pval", title = NULL,
                                            subtitle = NULL, pCutoff = 0.05,
                                            border = "full", caption = NULL)
ggsave(file.path(fig_dir, "fig13_volcano.png"), suppressWarnings(volcano),
       width = 22, height = 15, units = "cm", dpi = 150)

# GO enrichment lollipop on the meta-analysis DEGs (Figure 15)
go_part <- function(dataset, category) {
  res <- suppressMessages(rbioapi::rba_panther_enrich(
    genes = demo$meta_degs, organism = 1111708, annot_dataset = dataset, cutoff = 0.05
  ))$result
  res$category <- category
  res
}
all_go <- rbind(go_part("GO:0008150", "BP"), go_part("GO:0005575", "CC"), go_part("GO:0003674", "MF"))
all_go <- all_go[!is.na(all_go$term.id) & all_go$fold_enrichment > 0, ]
all_go <- all_go[!all_go$term.id %in% c("GO:0008150", "GO:0003674", "GO:0005575"), ]
all_go <- all_go[order(all_go$fdr), ]
for_plot <- data.frame(full_term = paste(all_go$term.id, all_go$term.label),
                       fold_enrichment = all_go$fold_enrichment,
                       fdr = all_go$fdr,
                       genes = all_go$number_in_list,
                       category = all_go$category)
for_plot <- utils::head(for_plot, 30)
go_plot <- ggplot(for_plot, aes(x = reorder(full_term, fold_enrichment), y = fold_enrichment)) +
  geom_segment(aes(xend = reorder(full_term, fold_enrichment), y = 0,
                   yend = fold_enrichment, colour = fdr), linewidth = 1.5) +
  scale_color_gradient2(low = "blue", high = "red", mid = "green",
                        midpoint = mean(for_plot$fdr), name = bquote(-log[10](FDR))) +
  geom_point(aes(colour = fdr, size = genes)) +
  labs(x = "GO terms", y = "Fold enrichment", size = "No. of genes") +
  facet_grid(rows = vars(category), scales = "free", space = "free") +
  coord_flip() +
  theme(axis.text.x = element_text(color = "grey20", size = 12),
        axis.text.y = element_text(color = "grey20", size = 12))
ggsave(file.path(fig_dir, "fig15_go_lollipop.png"), go_plot,
       width = 30, height = 20, units = "cm", dpi = 120)

message("Figures written to ", fig_dir, " (", nrow(all_go), " enriched GO terms)")
