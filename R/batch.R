#' Batch-effect handling methods
#'
#' How an optional `batch` column of the experimental design is used in each
#' study:
#' \describe{
#'   \item{`"model"`}{batch is added to the DE model, `~ batch + condition`.}
#'   \item{`"combat_seq"`}{the raw counts are corrected with
#'     [sva::ComBat_seq()] (keeping the condition effect), then analysed with
#'     `~ condition`.}
#'   \item{`"none"`}{batch is ignored.}
#' }
#'
#' @return A named character vector of method codes, named by their labels.
#' @export
#' @examples
#' batch_methods()
batch_methods <- function() {
  c("Batch in the DE model (~ batch + condition)" = "model",
    "ComBat-seq on the raw counts" = "combat_seq",
    "Ignore batch" = "none")
}

#' Can a study's batch be used?
#'
#' Checks, for one study, whether the requested batch method can be applied:
#' the study needs at least two batches, batch must not be confounded with
#' condition, the model `~ batch + condition` must leave residual degrees of
#' freedom, and ComBat-seq needs at least two samples per batch.
#'
#' @param design The study's rows of the experimental design.
#' @param method One of [batch_methods()].
#' @return `"applied"` if the batch method can be used, otherwise the reason
#'   why not (`"no batch column"`, `"ignored"`, `"single batch"`,
#'   `"confounded with condition"`, `"no residual degrees of freedom"` or
#'   `"a batch has only one sample"`).
#' @export
#' @examples
#' d <- data.frame(condition = c("c", "c", "t", "t", "c", "t"),
#'                 batch = c(1, 2, 1, 2, 1, 2))
#' batch_status(d, "model")
#' batch_status(d[1:4, ], "model")
batch_status <- function(design, method) {
  method <- match.arg(method, batch_methods())
  if (is.null(design$batch)) return("no batch column")
  if (method == "none") return("ignored")
  batch <- factor(design$batch)
  if (nlevels(batch) < 2) return("single batch")
  model <- stats::model.matrix(~ batch + factor(design$condition))
  rank <- qr(model)$rank
  if (rank < ncol(model)) return("confounded with condition")
  if (method == "model" && nrow(model) - rank < 1) return("no residual degrees of freedom")
  if (method == "combat_seq" && any(table(batch) < 2)) return("a batch has only one sample")
  "applied"
}

#' Correct raw counts with ComBat-seq
#'
#' @param counts Count matrix (genes x samples).
#' @param design Data frame with `batch` and `condition` columns, one row per
#'   sample.
#' @return The batch-corrected count matrix.
#' @export
#' @examples
#' set.seed(1)
#' counts <- matrix(rpois(600, 50), 100,
#'                  dimnames = list(paste0("g", 1:100), paste0("s", 1:6)))
#' design <- data.frame(batch = c(1, 1, 1, 2, 2, 2), condition = c("c", "t", "c", "t", "c", "t"))
#' corrected <- combat_seq_counts(counts, design)
combat_seq_counts <- function(counts, design) {
  corrected <- NULL
  utils::capture.output(
    corrected <- sva::ComBat_seq(counts, batch = as.character(design$batch),
                                 group = as.character(design$condition))
  )
  dimnames(corrected) <- dimnames(counts)
  corrected
}

# log2 values of one study with the batch effect removed for plotting:
# ComBat-seq counts when they were used, otherwise limma::removeBatchEffect()
# (keeping the condition effect). Unchanged when batch was not applied.
batch_corrected_log <- function(counts, design, status, method) {
  if (status == "applied" && method == "combat_seq") {
    return(log2(size_factor_normalise(combat_seq_counts(counts, design)) + 1))
  }
  log_counts <- log2(size_factor_normalise(counts) + 1)
  if (is.null(design$batch) || nlevels(factor(design$batch)) < 2 ||
      batch_status(design, "model") == "confounded with condition") {
    return(log_counts)
  }
  limma::removeBatchEffect(log_counts, batch = as.character(design$batch),
                           design = stats::model.matrix(~ factor(design$condition)))
}

#' Batch-effect diagnostics
#'
#' PCA and relative log expression (RLE) plots of each study. For a
#' `tma_result` they show the data before and after the batch correction, to
#' judge whether the correction removed the batch effect: "before" is the log2
#' normalised counts; "after" is the ComBat-seq corrected counts if ComBat-seq
#' was used, and otherwise the values from [limma::removeBatchEffect()]. For a
#' [batch_check()] they show the uploaded data before any analysis. Samples are
#' coloured by batch, or by condition when the design has no `batch` column.
#'
#' The PCA is computed separately in each study on its 500 most variable
#' genes, so that the differences between studies do not hide the batch
#' effect. The RLE plot shows, for each sample, the distribution of each
#' gene's log expression minus the gene's median in the study; boxes centred
#' on 0 with similar widths indicate comparable samples.
#'
#' @param result A `tma_result` from a design with a `batch` column, or a
#'   `tma_batch_check` from [batch_check()].
#' @return A ggplot.
#' @name batch_diagnostics
#' @examples
#' set.seed(1)
#' demo <- tma_demo_data()
#' demo$design$batch <- rep(c("b1", "b2"), 6)
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' plot_batch_pca(res)
#' plot_rle(res)
NULL

batch_data <- function(result) {
  if (inherits(result, "tma_batch_check")) return(list(Before = result$log_values))
  list(Before = log_values(result), After = result$batch$corrected_log)
}

# Colour of the samples in the diagnostic plots: batch, or condition without it
sample_groups <- function(design, samples) {
  rows <- match(samples, design$sample)
  if (is.null(design$batch)) design$condition[rows] else as.character(design$batch[rows])
}

#' @rdname batch_diagnostics
#' @export
plot_batch_pca <- function(result) {
  design <- result$design
  scores <- lapply(names(batch_data(result)), function(stage) {
    mat <- batch_data(result)[[stage]]
    do.call(rbind, lapply(result$studies, function(s) {
      samples <- design$sample[design$study == s]
      x <- mat[, samples, drop = FALSE]
      variances <- apply(x, 1, stats::var)
      top <- utils::head(order(variances, decreasing = TRUE), 500)
      pc <- stats::prcomp(t(x[top, , drop = FALSE]))
      data.frame(stage = stage, study = s, sample = samples,
                 PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                 batch = sample_groups(design, samples),
                 condition = design$condition[match(samples, design$sample)])
    }))
  })
  scores <- do.call(rbind, scores)
  scores$stage <- factor(scores$stage, levels = c("Before", "After"))
  ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, colour = batch, shape = condition)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_grid(stage ~ study, scales = "free") +
    ggplot2::labs(title = diagnostic_title("PCA of each study", result),
                  colour = if (is.null(design$batch)) "condition" else "batch")
}

#' @rdname batch_diagnostics
#' @export
plot_rle <- function(result) {
  design <- result$design
  boxes <- lapply(names(batch_data(result)), function(stage) {
    mat <- batch_data(result)[[stage]]
    do.call(rbind, lapply(result$studies, function(s) {
      samples <- design$sample[design$study == s]
      x <- mat[, samples, drop = FALSE]
      rle <- x - apply(x, 1, stats::median)
      q <- apply(rle, 2, stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
      data.frame(stage = stage, study = s, sample = samples,
                 batch = sample_groups(design, samples),
                 ymin = q[1, ], lower = q[2, ], middle = q[3, ], upper = q[4, ], ymax = q[5, ])
    }))
  })
  boxes <- do.call(rbind, boxes)
  boxes$stage <- factor(boxes$stage, levels = c("Before", "After"))
  fill_label <- if (is.null(design$batch)) "condition" else "batch"
  ggplot2::ggplot(boxes, ggplot2::aes(x = sample, fill = batch)) +
    ggplot2::geom_boxplot(ggplot2::aes(ymin = ymin, lower = lower, middle = middle,
                                       upper = upper, ymax = ymax), stat = "identity") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::facet_grid(stage ~ study, scales = "free_x") +
    ggplot2::labs(title = diagnostic_title("Relative log expression", result),
                  x = "Sample", y = "Relative log expression", fill = fill_label) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}

diagnostic_title <- function(what, result) {
  if (inherits(result, "tma_batch_check")) paste(what, "before the analysis")
  else paste(what, "before and after batch correction")
}

#' Check for batch effects before the analysis
#'
#' Looks for batch effects in each study before any differential expression
#' analysis. On the log2 normalised counts of each study's 500 most variable
#' genes, it runs a PCA and measures how much of the variance of the first
#' principal components is explained by condition and, if the design has a
#' `batch` column, by batch: for each component, the R-squared of a linear
#' model of the component on the factor, weighted by the component's share of
#' the variance, summed over the first `n_pcs` components (at most one fewer
#' than the number of samples).
#'
#' A study is flagged when batch explains as much variance as condition or
#' more, or when batch and condition are confounded. Without a `batch` column,
#' the check warns when the first principal component carries a large share of
#' the variance but is not explained by the condition, which can point to an
#' unrecorded batch effect. These flags are guides for looking at the PCA and
#' RLE plots ([plot_batch_pca()], [plot_rle()]), not statistical tests.
#'
#' @param counts Counts data frame with a `Geneid` column and one column per
#'   sample.
#' @param design Experimental design data frame with `sample`, `study` and
#'   `condition` columns, and optionally a `batch` column.
#' @param n_pcs Number of principal components to use.
#' @return An object of class `tma_batch_check`: a list with `table` (one row
#'   per study: `study`, `samples`, `condition_var` and `batch_var` (percent of
#'   variance explained), `pc1_var` (percent of variance on PC1),
#'   `pc1_condition_r2`, and `message`), `log_values` (log2 normalised counts of
#'   all samples), `design` and `studies`.
#' @export
#' @examples
#' demo <- tma_demo_data()
#' check <- batch_check(demo$counts, demo$design)
#' check$table
#' plot_batch_pca(check)
batch_check <- function(counts, design, n_pcs = 5) {
  studies <- filter_low_counts(split_studies(counts, design))
  has_batch <- "batch" %in% colnames(design)
  rows <- lapply(names(studies), function(id) {
    st <- studies[[id]]
    x <- log2(size_factor_normalise(st$counts) + 1)
    top <- utils::head(order(apply(x, 1, stats::var), decreasing = TRUE), 500)
    pc <- stats::prcomp(t(x[top, , drop = FALSE]))
    k <- min(n_pcs, ncol(x) - 1, ncol(pc$x))
    share <- pc$sdev^2 / sum(pc$sdev^2)
    explained <- function(f) {
      f <- factor(f)
      if (nlevels(f) < 2) return(0)
      r2 <- vapply(seq_len(k), function(i) summary(stats::lm(pc$x[, i] ~ f))$r.squared, numeric(1))
      100 * sum(share[seq_len(k)] * r2)
    }
    condition_var <- explained(st$design$condition)
    batch_var <- if (has_batch) explained(st$design$batch) else NA_real_
    pc1_condition_r2 <- summary(stats::lm(pc$x[, 1] ~ factor(st$design$condition)))$r.squared
    data.frame(study = id, samples = ncol(x),
               condition_var = condition_var, batch_var = batch_var,
               pc1_var = 100 * share[1], pc1_condition_r2 = pc1_condition_r2,
               message = batch_check_message(st$design, has_batch, condition_var, batch_var,
                                             100 * share[1], pc1_condition_r2))
  })
  log_values <- lapply(studies, function(st) log2(size_factor_normalise(st$counts) + 1))
  structure(list(table = do.call(rbind, rows),
                 log_values = combine_study_matrices(log_values, design$sample),
                 design = design, studies = names(studies)),
            class = "tma_batch_check")
}

# The advice shown for one study of a batch check
batch_check_message <- function(design, has_batch, condition_var, batch_var,
                                pc1_var, pc1_condition_r2) {
  if (has_batch) {
    status <- batch_status(design, "model")
    if (status == "single batch") return("Only one batch in this study.")
    if (status == "confounded with condition") {
      return("Batch is confounded with condition: the two effects cannot be separated.")
    }
    if (batch_var >= condition_var) {
      return(sprintf("Batch explains as much variance as condition or more (%.0f%% vs %.0f%%): consider 'Batch in the DE model' or ComBat-seq.",
                     batch_var, condition_var))
    }
    return(sprintf("Batch explains less variance than condition (%.0f%% vs %.0f%%).",
                   batch_var, condition_var))
  }
  if (pc1_var >= 30 && pc1_condition_r2 < 0.5) {
    return(sprintf("PC1 (%.0f%% of the variance) is not explained by the condition: there may be an unrecorded batch effect. Check the PCA, and add a 'batch' column if the samples were processed in batches.",
                   pc1_var))
  }
  "No sign of a batch effect: the main variation follows the condition."
}

#' @export
print.tma_batch_check <- function(x, ...) {
  cat("TMA batch-effect check (before the analysis)\n")
  for (i in seq_len(nrow(x$table))) {
    cat("  ", x$table$study[i], ": ", x$table$message[i], "\n", sep = "")
  }
  invisible(x)
}
