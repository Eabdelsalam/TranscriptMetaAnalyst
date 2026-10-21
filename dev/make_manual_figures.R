# Regenerates the result figures of manual/TMA_Manual.md from the demo data,
# drawn with the same package functions the app uses. Run from the package root:
#
#   Rscript -e 'devtools::load_all(); source("dev/make_manual_figures.R")'
#
# The enrichment figure queries the STRING web service and needs an internet connection.

demo <- tma_demo_data()
res <- run_tma_pipeline(demo$counts, demo$design)
fig_dir <- file.path("manual", "figures")

# PCA (Figure 10)
ggplot2::ggsave(file.path(fig_dir, "fig10_pca.png"), plot_pca(res),
                width = 15, height = 12, units = "cm", dpi = 150)

# Heatmap with the default 2 gene clusters (Figure 11)
heat <- heatmap_matrix(res)
grDevices::png(file.path(fig_dir, "fig11_heatmap.png"), width = 7, height = 10, units = "in", res = 120)
plot_heatmap_static(heat, res$design, heatmap_row_clusters(nrow(heat)), length(res$studies))
grDevices::dev.off()

# Volcano plot (Figure 13)
ggplot2::ggsave(file.path(fig_dir, "fig13_volcano.png"), plot_volcano(volcano_data(res)),
                width = 22, height = 15, units = "cm", dpi = 150)

# STRING enrichment of the meta-analysis DEGs against all tested genes (Figure 15)
string <- run_string_enrichment(res$meta_degs, species = 1148, background = res$genes)
ggplot2::ggsave(file.path(fig_dir, "fig15_string_lollipop.png"), plot_string_lollipop(string$table),
                width = 30, height = 22, units = "cm", dpi = 120)

message("Figures written to ", fig_dir, " (", nrow(string$table), " enriched STRING terms)")
