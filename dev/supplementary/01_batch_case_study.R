# S1. Batch-effect case study (Supplementary Section S1, Figure 3).
#
# Simulated data of dev/simulate_studies.R (6 studies x 6 samples, 3 vs 3,
# 30,000 genes, 10% truly DE, seed 2026), so the true DEGs are known. In every
# study the samples are split into two batches balanced with condition
# (batch 1: 2 control + 1 treated; batch 2: 1 control + 2 treated), and the
# counts of batch 2 are multiplied by gene-specific log-normal factors
# (meanlog 0, sdlog 0.5; sdlog 1 as a strong-batch scenario).
#
#   Rscript dev/supplementary/01_batch_case_study.R

source("dev/supplementary/00_setup.R")

sim <- simulate_studies(seed = seed)
truth <- sim$truth
contrast <- "Treated vs Cont"

# balanced batches and the batch-2 factors: with_batches() and
# add_batch_effect() in 00_setup.R
batch_design <- with_batches(sim$design)

analyse <- function(counts, design, batch_method, meta_method) {
  run_tma_analysis(counts, design, reference = "Cont", de_method = "DESeq2", meta_method = meta_method,
                   batch_method = batch_method, batch_plots = batch_method != "none")$results[[contrast]]
}

strategies <- c("Batch ignored" = "none", "Batch in the DE model" = "model", "ComBat-seq" = "combat_seq")
meta_methods <- c("Fisher" = "fisher", "Inverse normal" = "invnorm")

# (i) reference: original data, no batch effect
message("Reference analyses")
reference <- lapply(meta_methods, function(m) analyse(sim$counts, sim$design, "none", m))

rows <- list()
checks <- list()
examples <- list()
for (sdlog in c(0.5, 1)) {
  batch_counts <- add_batch_effect(sim$counts, batch_design, sdlog, seed = seed + 1)
  check <- batch_check(batch_counts, batch_design)
  checks[[as.character(sdlog)]] <- cbind(sdlog = sdlog, check$table)
  for (m in names(meta_methods)) {
    ref <- reference[[m]]
    rows[[length(rows) + 1]] <- cbind(sdlog = sdlog, strategy = "No batch effect (reference)", meta = m,
                                      evaluate_degs(ref$meta_degs, truth, ref$meta_degs))
    for (st in names(strategies)) {
      message(sprintf("sdlog %s, %s, %s", sdlog, st, m))
      res <- analyse(batch_counts, batch_design, strategies[[st]], meta_methods[[m]])
      rows[[length(rows) + 1]] <- cbind(sdlog = sdlog, strategy = st, meta = m,
                                        evaluate_degs(res$meta_degs, truth, ref$meta_degs))
      if (sdlog == 0.5 && m == "Fisher" && st == "ComBat-seq") examples$combat <- res
    }
  }
}
performance <- do.call(rbind, rows)
batch_checks <- do.call(rbind, checks)
batch_checks$flagged <- grepl("as much variance as condition or more|confounded", batch_checks$message)
save_table(performance, "S1_batch_performance.csv")
save_table(batch_checks, "S1_batch_check.csv")
print(performance)
print(batch_checks[, c("sdlog", "study", "condition_var", "batch_var", "flagged")])

# A pure between-study shift does not change the per-study p-values -----------
# The same simulation drawn twice from the same random numbers (inverse CDF of
# the negative binomial), once with and once without the gene-wise baseline
# shift of each study, so the only difference is the shift.
simulate_coupled <- function(shift_sd, n_studies = 6, n_genes = 30000, n_per_group = 3,
                             de_fraction = 0.1, seed = 2026) {
  set.seed(seed)
  base_mean <- exp(stats::rnorm(n_genes, mean = 5, sd = 1.5))
  dispersion <- 0.05 + 2 / sqrt(base_mean)
  n_de <- round(de_fraction * n_genes)
  lfc <- c(sample(c(-1.5, 1.5), n_de, replace = TRUE), rep(0, n_genes - n_de))
  condition <- rep(c("Cont", "Treated"), each = n_per_group)
  counts <- list()
  design <- list()
  for (s in LETTERS[seq_len(n_studies)]) {
    shift <- exp(stats::rnorm(n_genes, 0, 0.5))
    u <- matrix(stats::runif(n_genes * length(condition)), nrow = n_genes)
    study_mean <- base_mean * if (shift_sd > 0) shift else 1
    mu <- study_mean * 2^(outer(lfc, as.numeric(condition == "Treated")))
    samples <- paste0(s, "_", seq_along(condition))
    counts[[s]] <- matrix(stats::qnbinom(u, mu = mu, size = 1 / dispersion), nrow = n_genes,
                          dimnames = list(NULL, samples))
    design[[s]] <- data.frame(sample = samples, study = s, condition = condition)
  }
  list(counts = data.frame(Geneid = sprintf("gene%05d", seq_len(n_genes)), do.call(cbind, unname(counts)),
                           check.names = FALSE),
       design = do.call(rbind, unname(design)))
}
message("Between-study shift check")
with_shift <- simulate_coupled(0.5)
without_shift <- simulate_coupled(0)
res_with <- analyse(with_shift$counts, with_shift$design, "none", "fisher")
res_without <- analyse(without_shift$counts, without_shift$design, "none", "fisher")
shift_rows <- lapply(res_with$studies, function(s) {
  a <- res_with$de[[s]]
  b <- res_without$de[[s]]
  genes <- intersect(a$gene[!is.na(a$pvalue)], b$gene[!is.na(b$pvalue)])
  pa <- a$pvalue[match(genes, a$gene)]
  pb <- b$pvalue[match(genes, b$gene)]
  sa <- a$padj[match(genes, a$gene)] < 0.05
  sb <- b$padj[match(genes, b$gene)] < 0.05
  data.frame(study = s, genes = length(genes),
             spearman = stats::cor(pa, pb, method = "spearman"),
             median_abs_diff_log10p = stats::median(abs(log10(pa) - log10(pb))),
             same_call = 100 * mean(sa == sb, na.rm = TRUE))
})
shift_check <- do.call(rbind, shift_rows)
shift_summary <- data.frame(
  degs_with_shift = length(res_with$meta_degs), degs_without_shift = length(res_without$meta_degs),
  jaccard = length(intersect(res_with$meta_degs, res_without$meta_degs)) /
    length(union(res_with$meta_degs, res_without$meta_degs)))
save_table(shift_check, "S1_study_shift_pvalues.csv")
save_table(shift_summary, "S1_study_shift_summary.csv")
print(shift_check)
print(shift_summary)

# Figure 3: PCA and RLE before/after ComBat-seq for two studies, and TPR/FDP --
two_studies <- examples$combat
two_studies$studies <- two_studies$studies[1:2]
theme_fig <- ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(plot.title = ggplot2::element_blank(), legend.key.size = ggplot2::unit(0.35, "cm"))
pca <- plot_batch_pca(two_studies) + theme_fig + ggplot2::guides(shape = ggplot2::guide_legend(order = 2))
pca$layers[[1]]$aes_params$size <- 1.8
rle <- plot_rle(two_studies) + theme_fig +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
bars_data <- subset(performance, meta == "Fisher")
bars_data <- rbind(transform(bars_data, measure = "TPR (%)", value = tpr),
                   transform(bars_data, measure = "FDP (%)", value = fdp))
bars_data$strategy <- factor(bars_data$strategy,
                             levels = c("No batch effect (reference)", names(strategies)))
bars_data$scenario <- paste0("sdlog = ", bars_data$sdlog)
bars <- ggplot2::ggplot(bars_data, ggplot2::aes(x = strategy, y = value, fill = strategy)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", value)), vjust = -0.3, size = 2.2) +
  ggplot2::facet_grid(measure ~ scenario, scales = "free_y") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
  ggplot2::scale_fill_manual(values = c("grey60", "#D55E00", "#0072B2", "#009E73")) +
  ggplot2::labs(x = NULL, y = NULL, fill = NULL) + theme_fig +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                 legend.position = "bottom")
fig3 <- patchwork::wrap_plots(pca, rle, bars, ncol = 1, heights = c(1, 1, 1.1)) +
  patchwork::plot_annotation(tag_levels = "a")
save_figure(fig3, "fig3_batch", width = width_double, height = 21)

saveRDS(list(performance = performance, batch_checks = batch_checks, shift_check = shift_check,
             shift_summary = shift_summary, session = session_line()),
        out("S1_batch.rds"))
message("S1 done")
