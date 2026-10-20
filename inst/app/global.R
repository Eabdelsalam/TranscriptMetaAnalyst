#package version
VERSION <- paste0("v", utils::packageVersion("TranscriptMetaAnalyst"))

#attach the packages used directly in ui.R and server.R; the analysis code
#lives in the TranscriptMetaAnalyst package (already attached under
#devtools::load_all())
if (!"package:TranscriptMetaAnalyst" %in% search()) library(TranscriptMetaAnalyst)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyalert)
