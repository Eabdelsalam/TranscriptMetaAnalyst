# Shared setup of the supplementary analyses (dev/supplementary/). Every script
# sources this file; run the scripts from the package root, e.g.
#   Rscript dev/supplementary/01_batch_case_study.R
# or all of them with  Rscript dev/supplementary/run_all.R

suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))
source("dev/simulate_studies.R")

supp_dir <- file.path("dev", "supplementary")
out_dir <- file.path(supp_dir, "output")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
seed <- 2026

# Bioinformatics Advances figure widths (cm)
width_single <- 8.4
width_double <- 17.4

out <- function(name) file.path(out_dir, name)
save_table <- function(x, name) utils::write.csv(x, out(name), row.names = FALSE)

# Saves a ggplot or patchwork as 300 dpi PNG and as PDF
save_figure <- function(plot, name, width, height) {
  ggplot2::ggsave(file.path(fig_dir, paste0(name, ".png")), plot, width = width, height = height,
                  units = "cm", dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot, width = width, height = height,
                  units = "cm", bg = "white")
}

# Meta-analysis DEGs against a known truth or a reference set
evaluate_degs <- function(degs, truth, reference = NULL) {
  tp <- length(intersect(degs, truth))
  data.frame(degs = length(degs),
             tpr = 100 * tp / length(truth),
             fdp = if (length(degs) == 0) 0 else 100 * (length(degs) - tp) / length(degs),
             jaccard = if (is.null(reference)) NA_real_ else
               length(intersect(degs, reference)) / length(union(degs, reference)))
}

session_line <- function() {
  sprintf("R %s; TranscriptMetaAnalyst %s; DESeq2 %s; limma %s; edgeR %s",
          getRversion(), utils::packageVersion("TranscriptMetaAnalyst"), utils::packageVersion("DESeq2"),
          utils::packageVersion("limma"), utils::packageVersion("edgeR"))
}

# Within-study batches balanced with condition (per study the samples are
# Cont, Cont, Cont, Treated, Treated, Treated): batch 1 = 2 control + 1
# treated, batch 2 = 1 control + 2 treated
with_batches <- function(design) {
  design$batch <- rep(c("b1", "b1", "b2", "b1", "b2", "b2"), length(unique(design$study)))
  design
}

# Multiplies the counts of batch 2 by gene-specific log-normal factors (one set
# of factors per study), keeping integer counts
add_batch_effect <- function(counts, design, sdlog, seed) {
  set.seed(seed)
  out <- counts
  for (s in unique(design$study)) {
    factor <- stats::rlnorm(nrow(counts), meanlog = 0, sdlog = sdlog)
    cols <- design$sample[design$study == s & design$batch == "b2"]
    out[cols] <- round(as.matrix(counts[cols]) * factor)
  }
  out
}
