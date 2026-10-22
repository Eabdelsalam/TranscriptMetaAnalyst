# Runs all supplementary analyses in order and renders the supplementary file
# (output/TMA_Supplementary.docx). Needs an internet connection (GEO, KEGG,
# STRING, UniProt) and Chrome or Chromium (app timings and screenshots).
# Takes about an hour. From the package root:
#
#   NOT_CRAN=true Rscript dev/supplementary/run_all.R

scripts <- c("01_batch_case_study.R", "02_microarray_case_study.R", "03_rnaseq_demo_case_study.R",
             "04_performance.R", "05_figures.R", "06_bundle_example.R", "07_placeholders.R")
Sys.setenv(NOT_CRAN = "true")
for (s in scripts) {
  message("== ", s)
  status <- system2(file.path(R.home("bin"), "Rscript"), file.path("dev", "supplementary", s))
  if (status != 0) stop(s, " failed.", call. = FALSE)
}
rmarkdown::render(file.path("dev", "supplementary", "TMA_Supplementary.Rmd"),
                  output_dir = file.path("dev", "supplementary", "output"), quiet = TRUE)
message("Supplementary file written to dev/supplementary/output/TMA_Supplementary.docx")
