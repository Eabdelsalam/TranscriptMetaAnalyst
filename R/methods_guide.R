#' Methods Guide
#'
#' Content of the *Methods Guide* tab of the app: what TMA analyses, what the
#' studies must fulfil, and when to use which differential expression,
#' meta-analysis, batch and enrichment method, with references.
#'
#' @return A [shiny::tagList()] of `shinydashboard` boxes.
#' @export
#' @examples
#' guide <- methods_guide_ui()
#' class(guide)
methods_guide_ui <- function() {
  b <- shiny::strong
  li <- shiny::tags$li
  guide_box <- function(title, ..., status = "primary") {
    shinydashboard::box(title = title, status = status, solidHeader = TRUE, width = 12,
                        collapsible = TRUE, ...)
  }
  shiny::tagList(
    shiny::fluidRow(guide_box(
      "Scope",
      shiny::p("TMA does a ", b("gene-level"), " meta-analysis of ", b("bulk RNA-seq"),
               " (raw read counts per gene) and of ", b("single-colour microarray"),
               " data (normalised expression per probe or gene). In each study, every condition is compared with a reference condition; the p-values of the studies are then combined for each gene."),
      shiny::tags$ul(
        li("Transcript-level, isoform and alternative-splicing analyses are out of scope: the input must be one value per gene (RNA-seq) or per probe that is collapsed to its gene (microarrays)."),
        li("RNA-seq and microarray studies are analysed in separate tabs and cannot be combined in one meta-analysis."),
        li("Two-colour microarrays are not supported, because their channels are not comparable as independent samples.")
      )
    )),
    shiny::fluidRow(guide_box(
      "Study requirements",
      shiny::tags$ul(
        li(b("At least 2 samples (replicates) per condition in a study."),
           " A condition with fewer samples is left out of that study, and the study is not used for the comparisons that involve it."),
        li(b("At least 2 studies per comparison."),
           " A comparison with fewer studies is not meta-analysed. If no comparison can be meta-analysed, the analysis stops and says why."),
        li(b("Gene-level IDs."), " The same type of gene ID in every study, so that genes can be matched across studies. For microarrays, the probes are collapsed to genes (the probe with the highest mean expression is kept) and only the genes measured in every study are used."),
        li(b("Supported ID types for the enrichment:"), " see the table below. The app reports, for every enrichment, how many IDs the service recognised.")
      ),
      shiny::tags$table(
        class = "table table-condensed", style = "width: auto;",
        shiny::tags$thead(shiny::tags$tr(shiny::tags$th("ID type"), shiny::tags$th("PANTHER (GO)"),
                                         shiny::tags$th("KEGG"), shiny::tags$th("STRING"))),
        shiny::tags$tbody(
          lapply(names(id_types()), function(label) {
            kegg <- if (is.na(kegg_key_type(id_types()[[label]]))) "no" else "yes"
            shiny::tags$tr(shiny::tags$td(label), shiny::tags$td("yes"), shiny::tags$td(kegg), shiny::tags$td("yes"))
          })
        )
      )
    )),
    shiny::fluidRow(guide_box(
      "Differential expression (in each study)",
      shiny::p(b("RNA-seq."), " Genes with too few reads in every study are removed first (edgeR's filterByExpr [5]). Then:"),
      shiny::tags$ul(
        li(b("DESeq2"), " (default) [4]: negative-binomial model with shrunken dispersion estimates. Works well with few replicates per group."),
        li(b("edgeR"), " (quasi-likelihood F-test) [5]: also negative-binomial; the quasi-likelihood test keeps the false discovery rate under control with few replicates."),
        li(b("limma-voom"), " [6, 7]: linear models on log-counts with precision weights. Fast, and a good choice for larger studies or when library sizes differ a lot.")
      ),
      shiny::p("In a benchmark of RNA-seq tools with few replicates, DESeq2 and edgeR were among the best performers [2]; a review of best practices gives the same advice and discusses the number of replicates needed [1]. Choose one method for all studies unless a study clearly differs from the others (e.g. many more samples), and report the method in the methods section."),
      shiny::p(b("Microarrays."), " limma [7]: linear models on the log2 expression values (log2-transformed first when needed), every condition compared with the reference through contrasts, and empirical Bayes moderation with an intensity trend.")
    )),
    shiny::fluidRow(guide_box(
      "Meta-analysis (combining the studies)",
      shiny::tags$ul(
        li(b("Fisher"), " (default) [8]: sums -log(p) over the studies. A very strong signal in one study can be enough to make a gene significant."),
        li(b("Inverse normal"), " [3]: combines the studies' z-scores, weighted by the square root of each study's number of samples. It favours genes with consistent evidence across studies and gives larger studies more weight.")
      ),
      shiny::p("Both combine the raw p-values of the studies (a study whose adjusted p-value of a gene is missing, e.g. after DESeq2's independent filtering, does not contribute for that gene), and the combined p-values are adjusted with the Benjamini-Hochberg method. A gene is a meta-analysis DEG when its adjusted p-value is at most the significance level (0.05) and its fold change has the same sign in every study [3]. Use Fisher to find genes that respond strongly in at least some studies, and the inverse-normal method when you want genes that respond consistently across studies.")
    )),
    shiny::fluidRow(guide_box(
      "Batch effects (needs a 'batch' column in the design file)",
      shiny::tags$ul(
        li(b("Batch in the DE model"), " (recommended when there is a batch column): fits ~ batch + condition in each study, so the p-values account for the batch."),
        li(b("ComBat-seq"), " [9] (RNA-seq only): corrects the raw counts for batch, then fits ~ condition. The corrected counts are then treated as raw data, which can make p-values too small in small studies [10]."),
        li(b("Ignore batch"), " (default).")
      ),
      shiny::p("Batch can only be used in a study with at least two batches that are not confounded with the condition. The ", b("Check for batch effects"),
               " box shows, before the analysis, how much of the variation in each study is explained by batch and by condition. The PCA and heatmap can be drawn from batch-corrected values (ComBat-seq counts or limma's removeBatchEffect [7]); this changes only the plots.")
    )),
    shiny::fluidRow(guide_box(
      "Functional enrichment",
      shiny::tags$ul(
        li(b("PANTHER"), " [11]: GO biological process, molecular function and cellular component, for the organisms in PANTHER."),
        li(b("KEGG"), " [12]: KEGG pathways, for the organisms in KEGG; needs locus tags / KEGG gene IDs, NCBI Gene IDs or UniProt accessions."),
        li(b("STRING"), " [13]: GO, KEGG, Reactome, Pfam, InterPro and other annotation sources, for more than 12,000 organisms and most ID types.")
      ),
      shiny::p("The meta-analysis DEGs are compared with all genes tested in the analysis (the background), not with the whole genome, so that terms are not enriched only because the tested genes differ from the genome. All three services are queried online; no local annotation database is needed.")
    )),
    shiny::fluidRow(guide_box(
      "Reproducibility",
      shiny::p("After an analysis, the Comparisons box offers two downloads:"),
      shiny::tags$ul(
        li(b("Download reproducible analysis (ZIP)"), ": a ready-to-run bundle with the input data as used, a setup script that installs the same TMA release and checks the package versions, the analysis script with every setting of the run (analysis methods, plot display options, enrichment choices), the settings as JSON, the results of the run and the R session information. Run ", shiny::code("Rscript setup.R"), " and then ", shiny::code("Rscript analysis.R"), "; the re-run writes the same files to rerun/, for comparison with results/."),
        li(b("Download R script"), ": the analysis script alone, for the methods section; for uploaded data its input paths are placeholders.")
      ),
      shiny::p("Report the TMA version, the methods and the settings (all in the bundle's README.md and settings.json) in the methods section.")
    )),
    shiny::fluidRow(guide_box(
      "References", status = "info",
      shiny::tags$ol(
        li("Conesa, A. et al. A survey of best practices for RNA-seq data analysis. Genome Biology 17, 13 (2016)."),
        li("Schurch, N. J. et al. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? RNA 22, 839-851 (2016)."),
        li("Rau, A., Marot, G. & Jaffrezic, F. Differential meta-analysis of RNA-seq data from multiple studies. BMC Bioinformatics 15, 91 (2014)."),
        li("Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550 (2014)."),
        li("Chen, Y., Lun, A. T. L. & Smyth, G. K. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. F1000Research 5, 1438 (2016)."),
        li("Law, C. W., Chen, Y., Shi, W. & Smyth, G. K. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29 (2014)."),
        li("Ritchie, M. E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43, e47 (2015)."),
        li("Fisher, R. A. Statistical Methods for Research Workers. Oliver and Boyd, Edinburgh (1925)."),
        li("Zhang, Y., Parmigiani, G. & Johnson, W. E. ComBat-seq: batch effect adjustment for RNA-seq count data. NAR Genomics and Bioinformatics 2, lqaa078 (2020)."),
        li("Nygaard, V., Rodland, E. A. & Hovig, E. Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. Biostatistics 17, 29-39 (2016)."),
        li("Mi, H. et al. Protocol update for large-scale genome and gene function analysis with the PANTHER classification system (v.14.0). Nature Protocols 14, 703-721 (2019)."),
        li("Kanehisa, M. & Goto, S. KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Research 28, 27-30 (2000)."),
        li("Szklarczyk, D. et al. The STRING database in 2023. Nucleic Acids Research 51, D638-D646 (2023).")
      )
    ))
  )
}
