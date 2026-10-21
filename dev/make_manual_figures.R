# Regenerates the result figures of manual/TMA_Manual.md from the demo data,
# drawn with the same package functions the app uses. Run from the package root:
#
#   Rscript -e 'devtools::load_all(); source("dev/make_manual_figures.R")'
#
# The enrichment figure queries the STRING web service and needs an internet connection.

demo <- tma_demo_data()
res <- run_tma_pipeline(demo$counts, demo$design)
fig_dir <- file.path("manual", "figures")

# Meta-analysis p-value distribution (Figure 12)
grDevices::png(file.path(fig_dir, "fig12_pvalue_dist.png"), width = 6, height = 4.5, units = "in", res = 150)
plot_pvalue_hist(res)
grDevices::dev.off()

# PCA (Figure 13)
ggplot2::ggsave(file.path(fig_dir, "fig13_pca.png"), plot_pca(res),
                width = 15, height = 12, units = "cm", dpi = 150)

# Heatmap with the default 2 gene clusters (Figure 14)
heat <- heatmap_matrix(res)
grDevices::png(file.path(fig_dir, "fig14_heatmap.png"), width = 7, height = 10, units = "in", res = 120)
plot_heatmap_static(heat, res$design, heatmap_row_clusters(nrow(heat)), length(res$studies))
grDevices::dev.off()

# Venn diagram and UpSet plot of the DEGs of the studies and the meta-analysis (Figure 15)
sets <- deg_sets(res)
ggplot2::ggsave(file.path(fig_dir, "fig15a_venn.png"), plot_venn(sets),
                width = 15, height = 15, units = "cm", dpi = 150)
grDevices::png(file.path(fig_dir, "fig15b_upset.png"), width = 6, height = 3, units = "in", res = 150)
ComplexHeatmap::draw(plot_upset(sets))
grDevices::dev.off()

# Volcano plot (Figure 16)
ggplot2::ggsave(file.path(fig_dir, "fig16_volcano.png"), plot_volcano(volcano_data(res)),
                width = 22, height = 15, units = "cm", dpi = 150)

# Counts of the most significant meta-analysis DEG in every study (Figure 17b)
top_gene <- res$meta_degs[order(res$meta[res$meta_degs, "adjpval"])][1]
ggplot2::ggsave(file.path(fig_dir, "fig17b_counts_plot.png"), plot_gene_counts(res, top_gene),
                width = 20, height = 10, units = "cm", dpi = 150)

# STRING enrichment of the meta-analysis DEGs against all tested genes (Figure 18)
string <- run_string_enrichment(res$meta_degs, species = 1148, background = res$genes)
ggplot2::ggsave(file.path(fig_dir, "fig18_string_lollipop.png"), plot_string_lollipop(string$table),
                width = 30, height = 22, units = "cm", dpi = 120)

message("Figures written to ", fig_dir, " (", nrow(string$table), " enriched STRING terms)")
