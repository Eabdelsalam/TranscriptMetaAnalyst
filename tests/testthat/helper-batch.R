# Two studies of 8 samples (4 vs 4) in two balanced batches. 10% of the genes
# are DE (x1.8); 30% have a strong batch effect (x4 in batch b2).
simulate_batch_data <- function(n_genes = 1000, seed = 7) {
  set.seed(seed)
  base <- exp(stats::rnorm(n_genes, 5, 1))
  de <- seq_len(n_genes) <= 0.1 * n_genes
  batchy <- seq_len(n_genes) %in% sample(n_genes, 0.3 * n_genes)
  condition <- rep(c("Cont", "Treated"), each = 4)
  batch <- rep(c("b1", "b2"), 4)
  counts <- list()
  design <- list()
  for (s in c("S1", "S2")) {
    mu <- outer(base, rep(1, 8))
    mu[de, condition == "Treated"] <- mu[de, condition == "Treated"] * 1.8
    mu[batchy, batch == "b2"] <- mu[batchy, batch == "b2"] * 4
    samples <- paste0(s, "_", 1:8)
    counts[[s]] <- matrix(stats::rnbinom(length(mu), mu = mu, size = 10), n_genes,
                          dimnames = list(NULL, samples))
    design[[s]] <- data.frame(sample = samples, study = s, condition = condition, batch = batch)
  }
  list(counts = data.frame(Geneid = sprintf("g%04d", seq_len(n_genes)),
                           do.call(cbind, unname(counts)), check.names = FALSE),
       design = do.call(rbind, unname(design)),
       de_genes = sprintf("g%04d", which(de)),
       batch_genes = sprintf("g%04d", which(batchy)))
}

batch_cache <- new.env()
batch_data_set <- function() {
  if (is.null(batch_cache$data)) batch_cache$data <- simulate_batch_data()
  batch_cache$data
}
