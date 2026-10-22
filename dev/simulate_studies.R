# Simulated multi-study RNA-seq data, used by dev/benchmark.R and the
# batch-effect case study (dev/supplementary/01_batch_case_study.R).
#
# Simulated multi-study RNA-seq counts --------------------------------------
# Negative-binomial counts. 10% of genes are DE with the same direction in all
# studies (|log2FC| = 1.5); every study also has its own gene-wise baseline
# shift, as real independent experiments do.
simulate_studies <- function(n_studies = 6, n_genes = 30000, n_per_group = 3,
                             de_fraction = 0.1, seed = 2026) {
  set.seed(seed)
  base_mean <- exp(stats::rnorm(n_genes, mean = 5, sd = 1.5))
  dispersion <- 0.05 + 2 / sqrt(base_mean)
  n_de <- round(de_fraction * n_genes)
  lfc <- c(sample(c(-1.5, 1.5), n_de, replace = TRUE), rep(0, n_genes - n_de))

  study_ids <- LETTERS[seq_len(n_studies)]
  counts <- list()
  design <- list()
  for (s in study_ids) {
    study_mean <- base_mean * exp(stats::rnorm(n_genes, 0, 0.5))
    condition <- rep(c("Cont", "Treated"), each = n_per_group)
    samples <- paste0(s, "_", seq_along(condition))
    mu <- outer(study_mean, ifelse(condition == "Treated", 1, 0), function(m, t) m) *
      2^(outer(lfc, ifelse(condition == "Treated", 1, 0)))
    m <- matrix(stats::rnbinom(length(mu), mu = mu, size = 1 / dispersion),
                nrow = n_genes, dimnames = list(NULL, samples))
    counts[[s]] <- m
    design[[s]] <- data.frame(sample = samples, study = s, condition = condition)
  }
  genes <- sprintf("gene%05d", seq_len(n_genes))
  list(
    counts = data.frame(Geneid = genes, do.call(cbind, unname(counts)), check.names = FALSE),
    design = do.call(rbind, unname(design)),
    truth = genes[lfc != 0]
  )
}

