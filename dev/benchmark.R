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

# Simulated multi-study RNA-seq counts (dev/simulate_studies.R)
source("dev/simulate_studies.R")

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
