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
#' \donttest{
#' set.seed(1)
#' counts <- matrix(rpois(600, 50), 100,
#'                  dimnames = list(paste0("g", 1:100), paste0("s", 1:6)))
#' design <- data.frame(batch = c(1, 1, 1, 2, 2, 2), condition = c("c", "t", "c", "t", "c", "t"))
#' corrected <- combat_seq_counts(counts, design)
#' }
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
#' PCA and relative log expression (RLE) plots of each study before and after
#' the batch correction, to judge whether the correction removed the batch
#' effect. "Before" is the log2 normalised counts; "after" is the ComBat-seq
#' corrected counts if ComBat-seq was used, and otherwise the values from
#' [limma::removeBatchEffect()].
#'
#' The PCA is computed separately in each study on its 500 most variable
#' genes, so that the differences between studies do not hide the batch
#' effect. The RLE plot shows, for each sample, the distribution of each
#' gene's log expression minus the gene's median in the study; boxes centred
#' on 0 with similar widths indicate comparable samples.
#'
#' @param result A `tma_result` from a design with a `batch` column.
#' @return A ggplot.
#' @name batch_diagnostics
#' @examples
#' \donttest{
#' set.seed(1)
#' demo <- tma_demo_data()
#' demo$design$batch <- rep(c("b1", "b2"), 6)
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' plot_batch_pca(res)
#' plot_rle(res)
#' }
NULL

batch_data <- function(result) {
  list(Before = log2(result$norm_counts + 1), After = result$batch$corrected_log)
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
                 batch = as.character(design$batch[match(samples, design$sample)]),
                 condition = design$condition[match(samples, design$sample)])
    }))
  })
  scores <- do.call(rbind, scores)
  scores$stage <- factor(scores$stage, levels = c("Before", "After"))
  ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, colour = batch, shape = condition)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_grid(stage ~ study, scales = "free") +
    ggplot2::labs(title = "PCA of each study before and after batch correction")
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
                 batch = as.character(design$batch[match(samples, design$sample)]),
                 ymin = q[1, ], lower = q[2, ], middle = q[3, ], upper = q[4, ], ymax = q[5, ])
    }))
  })
  boxes <- do.call(rbind, boxes)
  boxes$stage <- factor(boxes$stage, levels = c("Before", "After"))
  ggplot2::ggplot(boxes, ggplot2::aes(x = sample, fill = batch)) +
    ggplot2::geom_boxplot(ggplot2::aes(ymin = ymin, lower = lower, middle = middle,
                                       upper = upper, ymax = ymax), stat = "identity") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::facet_grid(stage ~ study, scales = "free_x") +
    ggplot2::labs(title = "Relative log expression before and after batch correction",
                  x = "Sample", y = "Relative log expression") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}
