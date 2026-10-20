# Regenerates fixtures/demo-reference.txt, the recorded demo meta-analysis result
# used by test-demo-regression.R. Run from the package root after an intended
# change to the statistical pipeline or after upgrading DESeq2/metaRNASeq:
#
#   Rscript -e 'devtools::load_all(); source("tests/testthat/fixtures/make-demo-reference.R")'
#
# Review the diff of demo-reference.txt before committing it.

testthat::source_test_helpers("tests/testthat", env = environment())
demo <- demo_pipeline()

out_file <- file.path("tests", "testthat", "fixtures", "demo-reference.txt")
header <- c(
  "# Demo meta-analysis reference (meta-significant, direction-consistent DEGs)",
  "# Regenerate with tests/testthat/fixtures/make-demo-reference.R",
  paste0("# created: ", format(Sys.Date())),
  paste0("# DESeq2: ", utils::packageVersion("DESeq2")),
  paste0("# metaRNASeq: ", utils::packageVersion("metaRNASeq")),
  paste0("# consistent_genes: ", sum(abs(demo$fc_select$signFC) == 1, na.rm = TRUE))
)
writeLines(c(header, sort(demo$meta_degs)), out_file)
message("Wrote ", length(demo$meta_degs), " DEGs to ", out_file)
