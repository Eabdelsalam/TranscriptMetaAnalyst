# Benchmark of the meta-analysis pipeline (task 2.7).
#
# Compares the baseline pipeline (dev/bench_baseline_pipeline.R, the code
# before the Phase 2 refactor) with the current package pipeline, on
#   - the bundled demo data (3 studies x 4 samples x 3405 genes), and
#   - a simulated dataset of 6 studies x 6 samples (3 vs 3) x 30,000 genes.
# For each it reports the run time, the memory allocated during the run (main
# R process only), and the size of the result that the app keeps in memory.
#
# Run from the package root:
#   Rscript dev/benchmark.R
# The numbers are summarised in dev/BENCHMARKS.md.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(bench)
})
source("dev/bench_baseline_pipeline.R")

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
  list(
    counts = data.frame(Geneid = sprintf("gene%05d", seq_len(n_genes)),
                        do.call(cbind, unname(counts)), check.names = FALSE),
    design = do.call(rbind, unname(design))
  )
}

read_demo <- function() {
  list(
    counts = utils::read.delim(demo_file("counts")),
    design = utils::read.delim(demo_file("design"))
  )
}

# Locates the demo files before (data/) and after (inst/extdata) task 2.6
demo_file <- function(which) {
  name <- paste0("example_rnaseq_", which, ".txt")
  path <- system.file("extdata", name, package = "TranscriptMetaAnalyst")
  if (path == "") path <- system.file("data", name, package = "TranscriptMetaAnalyst")
  path
}

# Pipelines to compare ------------------------------------------------------
# "current, serial" is the fair memory comparison: bench::mark() only counts
# memory allocated in the main R process, not in parallel workers.
pipelines <- list(baseline = baseline_pipeline)
if (exists("run_tma_pipeline", mode = "function")) {
  pipelines$`current, serial` <- function(counts, design) {
    run_tma_pipeline(counts, design, BPPARAM = BiocParallel::SerialParam())
  }
  pipelines$`current, parallel` <- function(counts, design) run_tma_pipeline(counts, design)
}

run_benchmark <- function(data, label, iterations) {
  rows <- lapply(names(pipelines), function(p) {
    f <- pipelines[[p]]
    result <- NULL
    in_parallel <- grepl("parallel", p)
    bm <- bench::mark(result <- f(data$counts, data$design),
                      iterations = iterations, check = FALSE, filter_gc = FALSE,
                      memory = !in_parallel)
    data.frame(
      dataset = label,
      pipeline = p,
      median_time = format(bm$median),
      mem_alloc = if (in_parallel) "n/a" else format(bm$mem_alloc),
      result_size = format(structure(utils::object.size(result), class = "object_size"),
                           units = "auto"),
      meta_degs = length(result$meta_degs)
    )
  })
  do.call(rbind, rows)
}

message("R ", getRversion(), ", DESeq2 ", utils::packageVersion("DESeq2"),
        ", ", parallel::detectCores(), " cores")
results <- rbind(
  run_benchmark(read_demo(), "demo (3 x 4 samples, 3405 genes)", iterations = 3),
  run_benchmark(simulate_studies(), "simulated (6 x 6 samples, 30k genes)", iterations = 2)
)
print(results, row.names = FALSE)
