# Three studies with three conditions (Cont, Heat, Cold), 3 samples each.
# Genes 1-80 respond to Heat (x3), genes 81-160 to Cold (x0.3).
simulate_multi_data <- function(n_genes = 800, seed = 21) {
  set.seed(seed)
  base <- exp(stats::rnorm(n_genes, 5, 1))
  condition <- rep(c("Cont", "Heat", "Cold"), each = 3)
  counts <- list()
  design <- list()
  for (s in c("S1", "S2", "S3")) {
    mu <- outer(base * exp(stats::rnorm(n_genes, 0, 0.3)), rep(1, 9))
    mu[1:80, condition == "Heat"] <- mu[1:80, condition == "Heat"] * 3
    mu[81:160, condition == "Cold"] <- mu[81:160, condition == "Cold"] * 0.3
    samples <- paste0(s, "_", 1:9)
    counts[[s]] <- matrix(stats::rnbinom(length(mu), mu = mu, size = 20), n_genes,
                          dimnames = list(NULL, samples))
    design[[s]] <- data.frame(sample = samples, study = s, condition = condition)
  }
  list(counts = data.frame(Geneid = sprintf("g%04d", seq_len(n_genes)),
                           do.call(cbind, unname(counts)), check.names = FALSE),
       design = do.call(rbind, unname(design)),
       heat_genes = sprintf("g%04d", 1:80), cold_genes = sprintf("g%04d", 81:160))
}

multi_cache <- new.env()
multi_data_set <- function() {
  if (is.null(multi_cache$data)) multi_cache$data <- simulate_multi_data()
  multi_cache$data
}

# A tma_analysis wrapping one result, for app tests that mock the analysis
fake_analysis <- function(result) {
  contrast <- if (is.null(result$contrast)) "Treated vs Cont" else result$contrast
  structure(list(reference = "Cont", contrasts = contrast,
                 results = stats::setNames(list(result), contrast),
                 report = data.frame(contrast = contrast, study = result$studies, used = TRUE, reason = ""),
                 summary = data.frame(contrast = contrast, studies = length(result$studies),
                                      analysed = TRUE, reason = "")),
            class = "tma_analysis")
}
