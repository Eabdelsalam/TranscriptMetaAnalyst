# S5. Contents of the ready-to-run analysis bundle (Supplementary Section S5).
#
# Builds the bundle of the RNA-Seq demo analysis (DESeq2 + Fisher), as the
# app's "Download reproducible analysis (ZIP)" does, keeps the ZIP in
# output/, lists its files and re-runs its analysis.R with
# Rscript to show that the re-run reproduces the tables of the run.
#
#   Rscript dev/supplementary/06_bundle_example.R

source("dev/supplementary/00_setup.R")

demo <- tma_demo_data()
settings <- list(de_method = "DESeq2", meta_method = "fisher", batch_method = "none", batch_plots = FALSE)
start <- Sys.time()
analysis <- run_tma_analysis(demo$counts, demo$design, reference = "Cont", de_method = "DESeq2",
                             meta_method = "fisher", batch_method = "none")
run_time <- as.numeric(difftime(Sys.time(), start, units = "secs"))
zip_file <- out("S5_example_bundle.zip")
unlink(zip_file)
write_tma_bundle(zip_file, analysis, demo, settings, source = "demo",
                 enrichment = list(`Treated vs Cont` = list(params = list(id_type = "kegg", go = 1111708,
                                                                           kegg = "syn", string = 1148,
                                                                           run = character(0)))),
                 run_time = run_time)

unpacked <- out("S5_example_bundle")
unlink(unpacked, recursive = TRUE)
zip::unzip(zip_file, exdir = unpacked)
bundle <- file.path(unpacked, "tma_analysis")
files <- list.files(bundle, recursive = TRUE)
listing <- data.frame(file = files, bytes = file.size(file.path(bundle, files)))
save_table(listing, "S5_bundle_files.csv")
file.copy(file.path(bundle, "README.md"), out("S5_example_README.md"), overwrite = TRUE)

# re-run analysis.R as a user would (the setup step installs the release from
# GitHub, so it is skipped here; this source tree is loaded instead)
lines <- readLines(file.path(bundle, "analysis.R"))
lines[lines == "library(TranscriptMetaAnalyst)"] <- sprintf("pkgload::load_all(%s, quiet = TRUE)",
                                                           deparse(normalizePath(".")))
writeLines(lines, file.path(bundle, "analysis_local.R"))
rerun_start <- Sys.time()
log <- withr::with_dir(bundle, system2(file.path(R.home("bin"), "Rscript"), "analysis_local.R",
                                       stdout = TRUE, stderr = TRUE))
rerun_time <- as.numeric(difftime(Sys.time(), rerun_start, units = "secs"))
unlink(file.path(bundle, "analysis_local.R"))
tables <- grep("\\.csv$", list.files(file.path(bundle, "results")), value = TRUE)
identical_tables <- vapply(tables, function(t) {
  identical(readLines(file.path(bundle, "rerun", t)), readLines(file.path(bundle, "results", t)))
}, logical(1))
comparison <- data.frame(table = tables, identical = identical_tables)
save_table(comparison, "S5_rerun_comparison.csv")
print(comparison)

# the ZIP, its README and the file list are kept; the unpacked copy is not
unlink(unpacked, recursive = TRUE)

saveRDS(list(listing = listing, comparison = comparison, run_time = run_time, rerun_time = rerun_time,
             zip_bytes = file.size(zip_file), session = session_line()),
        out("S5_bundle.rds"))
message("S5 done")
