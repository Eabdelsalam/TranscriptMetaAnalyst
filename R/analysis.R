#' Run the TMA RNA-seq meta-analysis
#'
#' Compares every condition with a reference condition ("all vs reference")
#' across the studies. Each study is fitted once with all its conditions, and
#' each contrast is extracted from that fit. For each contrast, the p-values of
#' the studies are combined with the chosen meta-analysis method, genes whose
#' direction conflicts between studies are removed, and the result is
#' collected in a `tma_result` object.
#'
#' Replicate rule: a condition with fewer than 2 samples in a study is removed
#' from that study (its samples are not used in the fit), so the study is
#' excluded from the contrasts that involve that condition. A study without at
#' least 2 samples of the reference condition is not used at all. A contrast is
#' meta-analysed only if at least 2 studies remain for it. If no contrast can
#' be meta-analysed, the function stops with an error that explains why. The
#' `report` element lists, for every contrast and study, whether the study was
#' used and why not.
#'
#' Before the analysis, genes with low counts in every study are removed
#' ([filter_low_counts()]).
#'
#' @param counts Counts data frame with a `Geneid` column and one column per
#'   sample.
#' @param design Experimental design data frame with `sample`, `study` and
#'   `condition` columns, and optionally a `batch` column.
#' @param reference The reference condition. `NULL` picks one with
#'   [default_reference()].
#' @param alpha Significance level for the per-study and meta-analysis tests.
#' @param de_method DE method for all studies (one of [de_methods()]), or a
#'   named character vector with one method per study.
#' @param meta_method Meta-analysis method, `"fisher"` or `"invnorm"` (see
#'   [run_meta()]).
#' @param batch_method How to use an optional `batch` column of `design`, one
#'   of [batch_methods()]. By default batch is added to the DE model in every
#'   study where that is possible (see [batch_status()]).
#' @param batch_plots If `TRUE`, the PCA and heatmap of the results use the
#'   batch-corrected values (see [batch_diagnostics]).
#' @param filter If `TRUE`, remove genes with low counts in every study
#'   before the analysis (see [filter_low_counts()]).
#' @param BPPARAM BiocParallel back-end for the per-study analyses, see
#'   [tma_bpparam()].
#' @return `run_tma_analysis()` returns an object of class `tma_analysis`, a
#'   list with `reference`, `contrasts` (all "condition vs reference"
#'   contrasts of the design), `results` (a named list with one `tma_result`
#'   per contrast that could be meta-analysed), `report` (one row per contrast
#'   and study: `contrast`, `study`, `used`, `reason`) and `summary` (one row
#'   per contrast: `contrast`, `studies`, `analysed`, `reason`).
#'
#'   `run_tma_pipeline()` is the same analysis for a design with exactly two
#'   conditions and returns the single `tma_result`, a list with
#'   \describe{
#'     \item{contrast, reference}{the contrast ("B vs A") and its reference.}
#'     \item{design}{the design rows of the samples used in this contrast.}
#'     \item{studies, genes}{study names and the names of the analysed genes.}
#'     \item{de}{per-study DE tables of this contrast (see [run_deseq2()]).}
#'     \item{norm_counts}{normalised counts of the samples used (genes x samples).}
#'     \item{fold_changes}{log2 fold changes (genes x studies).}
#'     \item{meta}{meta-analysis table: `TestStatistic`, `rawpval`, `adjpval`
#'       and `significant` ("yes"/"no"), one row per gene.}
#'     \item{fc_select}{per-study DE flags, the meta-analysis flag `DE.meta`,
#'       fold changes and the direction `signFC` of each gene.}
#'     \item{meta_degs}{meta-significant, direction-consistent DEGs.}
#'     \item{study_degs}{direction-consistent DEGs of each study.}
#'     \item{idd_irr}{integration-driven discoveries and lost DEGs
#'       ([metaRNASeq::IDD.IRR()]).}
#'     \item{de_methods}{the DE method used for each study.}
#'     \item{batch}{`NULL` without a `batch` column; otherwise a list with the
#'       `method` requested, its `status` in each study (see
#'       [batch_status()]), the batch-corrected log2 values
#'       `corrected_log`, and `plots` (use them in the PCA and heatmap).}
#'     \item{meta_method}{the meta-analysis method used.}
#'     \item{report}{the rows of the analysis report for this contrast.}
#'     \item{alpha}{the significance level used.}
#'   }
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' res <- run_tma_pipeline(demo$counts, demo$design)
#' length(res$meta_degs)
#' res2 <- run_tma_pipeline(demo$counts, demo$design, de_method = "edgeR",
#'                          meta_method = "invnorm")
#'
#' # three conditions: study C's treated samples become "High". "Treated vs
#' # Cont" is meta-analysed with studies A and B; "High vs Cont" exists only in
#' # study C, so it is not meta-analysed (the replicate rule needs 2 studies)
#' design3 <- demo$design
#' design3$condition[design3$study == "C" & design3$condition == "Treated"] <- "High"
#' analysis <- run_tma_analysis(demo$counts, design3, reference = "Cont")
#' analysis
#' analysis_messages(analysis)
#' }
run_tma_analysis <- function(counts, design, reference = NULL, alpha = 0.05,
                             de_method = "DESeq2", meta_method = "fisher",
                             batch_method = "model", batch_plots = FALSE, filter = TRUE,
                             BPPARAM = tma_bpparam(length(unique(design$study)))) {
  if (is.null(reference)) reference <- default_reference(design$condition)
  conditions <- sort(unique(as.character(design$condition)))
  if (!reference %in% conditions) {
    stop("The reference condition '", reference, "' is not in the design.", call. = FALSE)
  }
  others <- setdiff(conditions, reference)
  contrasts <- paste(others, "vs", reference)
  meta_method <- unname(match.arg(meta_method, meta_methods()))
  batch_method <- match.arg(batch_method, batch_methods())

  # replicate rule: keep the conditions with at least 2 samples in each study
  studies <- split_studies(counts, design)
  report <- replicate_report(design, reference, others)
  fit_studies <- list()
  for (id in names(studies)) {
    sizes <- table(as.character(studies[[id]]$design$condition))
    kept <- names(sizes)[sizes >= 2]
    if (!reference %in% kept || length(kept) < 2) next
    st <- studies[[id]]
    keep <- st$design$condition %in% kept
    fit_studies[[id]] <- list(counts = st$counts[, keep, drop = FALSE],
                              design = st$design[keep, , drop = FALSE])
  }
  summary <- contrast_summary(report, contrasts)
  if (!any(summary$analysed)) {
    stop("No contrast can be meta-analysed. ", paste(summary$reason, collapse = " "), call. = FALSE)
  }

  if (filter) fit_studies <- filter_low_counts(fit_studies)
  de_methods_used <- study_methods(de_method, names(fit_studies))
  de <- run_de_studies(fit_studies, methods = de_methods_used, alpha = alpha,
                       batch_method = batch_method, reference = reference, BPPARAM = BPPARAM)
  batch <- NULL
  if ("batch" %in% colnames(design)) {
    status <- vapply(fit_studies, function(st) batch_status(st$design, batch_method), character(1))
    corrected <- lapply(names(fit_studies), function(id) {
      batch_corrected_log(fit_studies[[id]]$counts, fit_studies[[id]]$design, status[[id]], batch_method)
    })
    batch <- list(method = batch_method, status = status,
                  corrected_log = do.call(cbind, unname(corrected)), plots = isTRUE(batch_plots))
  }

  results <- list()
  for (i in seq_along(others)) {
    if (!summary$analysed[i]) next
    used <- report$study[report$contrast == contrasts[i] & report$used]
    results[[contrasts[i]]] <- build_tma_result(
      fit_studies[used], lapply(de[used], function(d) d$tables[[contrasts[i]]]),
      design, conditions = c(reference, others[i]), contrast = contrasts[i],
      reference = reference, alpha = alpha, meta_method = meta_method,
      de_methods = de_methods_used[used], batch = batch,
      report = report[report$contrast == contrasts[i], , drop = FALSE]
    )
  }

  structure(list(reference = reference, contrasts = contrasts, results = results,
                 report = report, summary = summary),
            class = "tma_analysis")
}

#' @rdname run_tma_analysis
#' @export
run_tma_pipeline <- function(counts, design, alpha = 0.05, de_method = "DESeq2",
                             meta_method = "fisher", batch_method = "model",
                             batch_plots = FALSE, filter = TRUE, reference = NULL,
                             BPPARAM = tma_bpparam(length(unique(design$study)))) {
  n_conditions <- length(unique(design$condition))
  if (n_conditions != 2) {
    stop("The design has ", n_conditions, " conditions; use run_tma_analysis() to compare them with a reference.",
         call. = FALSE)
  }
  analysis <- run_tma_analysis(counts, design, reference = reference, alpha = alpha,
                               de_method = de_method, meta_method = meta_method,
                               batch_method = batch_method, batch_plots = batch_plots,
                               filter = filter, BPPARAM = BPPARAM)
  analysis$results[[1]]
}

#' Default reference condition
#'
#' Picks a condition whose name looks like a control (control, ctrl, cont,
#' wt, wild type, mock, untreated; case-insensitive, at the start of the
#' name), or else the first condition in alphabetical order.
#'
#' @param conditions The `condition` column of the design.
#' @return One condition name.
#' @export
#' @examples
#' default_reference(c("Treated", "Control", "Heat"))
#' default_reference(c("Heat", "Cold"))
default_reference <- function(conditions) {
  lvls <- sort(unique(as.character(conditions)))
  looks_like_control <- grepl("^(control|ctrl|cont|wt|wild[ _.-]?type|mock|untreated|untr)",
                              lvls, ignore.case = TRUE)
  if (any(looks_like_control)) lvls[looks_like_control][1] else lvls[1]
}

# One row per contrast and study: is the study used, and if not, why
replicate_report <- function(design, reference, others) {
  rows <- list()
  for (lvl in others) {
    for (id in unique(design$study)) {
      cond <- as.character(design$condition[design$study == id])
      n_ref <- sum(cond == reference)
      n_lvl <- sum(cond == lvl)
      reason <- if (n_ref == 0) {
        sprintf("has no samples of the reference condition '%s'", reference)
      } else if (n_ref < 2) {
        sprintf("has only 1 sample of the reference condition '%s' (at least 2 are needed)", reference)
      } else if (n_lvl == 0) {
        sprintf("has no samples of '%s'", lvl)
      } else if (n_lvl < 2) {
        sprintf("has only 1 sample of '%s' (at least 2 are needed)", lvl)
      } else {
        ""
      }
      rows[[length(rows) + 1]] <- data.frame(contrast = paste(lvl, "vs", reference), study = id,
                                             used = reason == "", reason = reason)
    }
  }
  do.call(rbind, rows)
}

# One row per contrast: number of studies used, and why a contrast is not analysed
contrast_summary <- function(report, contrasts) {
  rows <- lapply(contrasts, function(con) {
    used <- report$study[report$contrast == con & report$used]
    analysed <- length(used) >= 2
    reason <- if (analysed) "" else {
      sprintf("'%s' was not meta-analysed: %s at least 2 studies are needed.", con,
              if (length(used) == 0) "no study has both conditions with at least 2 samples;"
              else sprintf("only 1 study (%s) has both conditions with at least 2 samples;", used))
    }
    data.frame(contrast = con, studies = length(used), analysed = analysed, reason = reason)
  })
  do.call(rbind, rows)
}

# The tma_result of one contrast
build_tma_result <- function(studies, tables, design, conditions, contrast, reference,
                             alpha, meta_method, de_methods, batch, report) {
  # samples of the two compared conditions in the studies used
  studies <- lapply(studies, function(st) {
    keep <- st$design$condition %in% conditions
    list(counts = st$counts[, keep, drop = FALSE], design = st$design[keep, , drop = FALSE])
  })
  samples <- unlist(lapply(studies, function(st) st$design$sample), use.names = FALSE)
  design <- design[design$sample %in% samples, , drop = FALSE]
  genes <- rownames(studies[[1]]$counts)
  study_ids <- names(studies)
  de_column <- function(col) {
    m <- do.call(cbind, lapply(tables, `[[`, col))
    dimnames(m) <- list(genes, study_ids)
    m
  }
  pvalues <- de_column("pvalue")
  padj <- de_column("padj")
  fold_changes <- de_column("log2_fc")

  nrep <- vapply(studies, function(st) ncol(st$counts), integer(1))
  meta <- run_meta(pvalues, padj, method = meta_method, nrep = nrep, alpha = alpha)
  sign_fc <- consistent_direction(fold_changes)
  consistent <- abs(sign_fc) == 1 & !is.na(sign_fc)
  study_flags <- ifelse(padj <= alpha, 1, 0)

  fc_select <- data.frame(study_flags,
                          DE.meta = ifelse(meta$adjpval <= alpha, 1, 0),
                          fold_changes,
                          signFC = sign_fc,
                          row.names = genes)
  meta_degs <- select_meta_degs(fc_select)
  study_degs <- lapply(stats::setNames(study_ids, study_ids), function(s) {
    genes[which(consistent & study_flags[, s] == 1)]
  })
  meta$significant <- ifelse(meta$adjpval <= alpha & consistent, "yes", "no")
  rownames(meta) <- genes
  if (!is.null(batch)) {
    batch$status <- batch$status[study_ids]
    batch$corrected_log <- batch$corrected_log[, design$sample, drop = FALSE]
  }

  structure(list(
    contrast = contrast,
    reference = reference,
    design = design,
    studies = study_ids,
    genes = genes,
    de = tables,
    norm_counts = combine_study_matrices(
      lapply(studies, function(st) size_factor_normalise(st$counts)), design$sample
    ),
    fold_changes = fold_changes,
    meta = meta,
    fc_select = fc_select,
    meta_degs = meta_degs,
    study_degs = study_degs,
    idd_irr = metaRNASeq::IDD.IRR(meta_degs, unname(study_degs)),
    de_methods = de_methods,
    batch = batch,
    meta_method = meta_method,
    report = report,
    alpha = alpha
  ), class = "tma_result")
}

#' Messages about excluded studies and contrasts
#'
#' Turns the report of a `tma_analysis` into sentences for the user: which
#' studies were left out of which contrast and why, and which contrasts could
#' not be meta-analysed.
#'
#' @param analysis A `tma_analysis` from [run_tma_analysis()].
#' @return A character vector (empty when every study was used in every
#'   contrast).
#' @export
#' @examples
#' \donttest{
#' demo <- tma_demo_data()
#' design <- demo$design
#' design$condition[design$sample == "Sample2"] <- "High"
#' # study A keeps only one "Treated" sample, so it is left out
#' analysis_messages(run_tma_analysis(demo$counts, design[design$condition != "High", ]))
#' }
analysis_messages <- function(analysis) {
  skipped <- analysis$report[!analysis$report$used, , drop = FALSE]
  c(sprintf("Study %s %s, so it is not used for '%s'.", skipped$study, skipped$reason, skipped$contrast),
    analysis$summary$reason[!analysis$summary$analysed])
}

#' @export
print.tma_analysis <- function(x, ...) {
  cat("TMA meta-analysis: all conditions vs", x$reference, "\n")
  for (con in x$contrasts) {
    res <- x$results[[con]]
    if (is.null(res)) {
      cat("  ", con, ": not meta-analysed\n", sep = "")
    } else {
      cat("  ", con, ": ", length(res$studies), " studies, ", length(res$meta_degs),
          " meta-analysis DEGs\n", sep = "")
    }
  }
  msgs <- analysis_messages(x)
  if (length(msgs) > 0) cat(paste0("  - ", msgs, collapse = "\n"), "\n")
  invisible(x)
}
