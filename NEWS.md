# TranscriptMetaAnalyst 2.1.0

## New features

* **Download reproducible analysis (ZIP)** on the RNA-Seq and Microarray tabs:
  a ready-to-run bundle of an analysis with the input data as used, `setup.R`
  (installs the same TMA release and reports package versions that differ from
  the original run), `analysis.R` (every setting of the run, including the plot
  display options and enrichment choices of each comparison; relative paths),
  `settings.json`, the results of the run and the R session information. GEO
  imports include the downloaded series files and the sample-to-condition
  mapping, so the bundle runs offline.
* New functions `write_tma_bundle()`, `save_tma_results()`,
  `tma_plot_options()` and `run_tma_enrichment()`. The app's downloads, the
  R script and the bundle use the same functions, so a re-run gives the same
  files.
* **Download R script** now also records the plot options of every
  comparison and re-runs the enrichment services used in the app.

# TranscriptMetaAnalyst 2.0.0

Revision in response to the reviewers' comments.

## Bug fixes

* The heatmap and PCA used values attached to the wrong genes or samples; the
  volcano plot used the wrong fold changes. All three now use correctly aligned
  matrices.
* Meta-analysis p-values below 1e-16 are kept exact instead of being rounded to
  0, and the inverse-normal z-scores no longer become infinite.
* GO and KEGG enrichment use the meta-analysis DEGs only.
* Closing a browser window stops the app only when it was started with
  `runTMA()`, not on a shared Shiny Server.

## New features

* Differential expression with DESeq2, edgeR or limma-voom, for all studies or
  per study, and meta-analysis with Fisher's or the inverse-normal method.
* Batch effects: an optional `batch` column, handled in the DE model or with
  ComBat-seq, a check for batch effects before the analysis, and PCA and RLE
  diagnostics before and after the correction.
* More than two conditions: every condition is compared with a reference
  condition, with one results tab per comparison. Conditions with fewer than 2
  samples in a study and comparisons with fewer than 2 studies are left out,
  and the app says why.
* Plot options (DEG selection, labels, colours, clustering, cutoffs) and
  downloads as PDF, PNG or TIFF.
* Enrichment against all tested genes, a choice of gene ID type, STRING
  enrichment and an ID-mapping report for every service.
* A reproducible R script of every analysis.
* A Microarray Meta-Analysis tab for normalised single-colour arrays, uploaded
  as files or imported from GEO, analysed with limma, with its own demo data.
* A Methods Guide tab and a vignette.
* The RNA-Seq Meta-Analysis tab is laid out like the Microarray tab and has a
  Load demo data button; the separate Run Demo Analysis tab was removed.

## Internal changes

* The app is in `inst/app`, its logic is in package functions and Shiny
  modules, and nothing is installed at run time.
* Studies are analysed in parallel and uploaded analyses run in the
  background; genes that are low in every study are filtered out first.
* Unit tests for the statistical functions and shinytest2 smoke tests for every
  tab.
