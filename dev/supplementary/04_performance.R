# S4. Performance (Supplementary Section S4).
#
# Machine details, the pipeline benchmarks of dev/BENCHMARKS.md, and the time
# from clicking "Load demo data" to the results in the app, for the RNA-Seq and
# the microarray demo (3 repeats each, a fresh app session each time, headless
# Chrome via shinytest2). Needs Chrome or Chromium.
#
#   NOT_CRAN=true Rscript dev/supplementary/04_performance.R

source("dev/supplementary/00_setup.R")
Sys.setenv(NOT_CRAN = "true")
options(chromote.timeout = 60)

# machine ---------------------------------------------------------------------
cpu <- tryCatch(sub(".*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1]),
                error = function(e) NA)
ram_gb <- tryCatch(as.numeric(sub("\\D*(\\d+).*", "\\1", grep("^MemTotal", readLines("/proc/meminfo"),
                                                                value = TRUE))) / 1024^2,
                   error = function(e) NA)
machine <- data.frame(
  item = c("Operating system", "CPU", "Logical cores", "RAM (GB)", "R", "TranscriptMetaAnalyst",
           "DESeq2", "edgeR", "limma", "metaRNASeq", "shiny"),
  value = c(utils::sessionInfo()$running, cpu, parallel::detectCores(), sprintf("%.1f", ram_gb),
            as.character(getRversion()), as.character(utils::packageVersion("TranscriptMetaAnalyst")),
            as.character(utils::packageVersion("DESeq2")), as.character(utils::packageVersion("edgeR")),
            as.character(utils::packageVersion("limma")), as.character(utils::packageVersion("metaRNASeq")),
            as.character(utils::packageVersion("shiny"))))
save_table(machine, "S4_machine.csv")

# pipeline benchmarks of dev/BENCHMARKS.md -------------------------------------
bench_lines <- readLines("dev/BENCHMARKS.md")
rows <- grep("^\\| (Demo|Simulated) \\|", bench_lines, value = TRUE)
cells <- lapply(strsplit(rows, "\\|"), function(x) trimws(x[-1]))
benchmarks <- as.data.frame(do.call(rbind, cells))
names(benchmarks) <- c("dataset", "pipeline", "median_time", "memory", "result_size", "meta_degs")
save_table(benchmarks, "S4_benchmarks.csv")

# time from click to results in the app ---------------------------------------
pkg <- normalizePath(".")
app_fun <- function() {
  pkgload::load_all(PKG, quiet = TRUE)
  shiny::shinyAppDir(TranscriptMetaAnalyst:::tma_app_dir())
}
body(app_fun) <- do.call(substitute, list(body(app_fun), list(PKG = pkg)))

time_demo <- function(tab, button) {
  app <- shinytest2::AppDriver$new(app_fun, width = 1400, height = 1000,
                                   load_timeout = 120 * 1000, timeout = 300 * 1000)
  on.exit(app$stop())
  app$run_js(sprintf("document.querySelector('.sidebar-menu a[data-value=\\'%s\\']').click();", tab))
  app$wait_for_idle()
  start <- Sys.time()
  app$click(button, wait_ = FALSE)
  app$wait_for_js("document.body.innerText.includes('genes were identified as DEGs')",
                  timeout = 600 * 1000, interval = 100)
  as.numeric(difftime(Sys.time(), start, units = "secs"))
}
app_times <- do.call(rbind, lapply(1:3, function(i) {
  data.frame(repeat_no = i,
             rnaseq_demo_s = time_demo("rnaseq", "rnaseq_demo"),
             microarray_demo_s = time_demo("microarray", "array_demo"))
}))
save_table(app_times, "S4_app_times.csv")
print(app_times)

saveRDS(list(machine = machine, benchmarks = benchmarks, app_times = app_times,
             rnaseq_median = stats::median(app_times$rnaseq_demo_s),
             microarray_median = stats::median(app_times$microarray_demo_s),
             session = session_line()),
        out("S4_performance.rds"))
message("S4 done")
