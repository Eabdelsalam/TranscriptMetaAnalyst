#package version
VERSION <- paste0("v", utils::packageVersion("TranscriptMetaAnalyst"))

#attach the packages used by ui.R and server.R; they are all package imports,
#so they are installed together with TranscriptMetaAnalyst
required_packages <- c("shiny", "shinydashboard", "shinyjs", "shinyalert", "shinycssloaders", "DESeq2", "metaRNASeq", "dplyr", "ggplot2", "ggrepel", "ggvenn", "ComplexHeatmap", "pheatmap", "EnhancedVolcano", "DT", "rbioapi", "clusterProfiler", "heatmaply")
invisible(lapply(required_packages, function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE, warn.conflicts = FALSE))
}))
