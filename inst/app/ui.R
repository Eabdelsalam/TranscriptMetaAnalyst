# Define UI for MetaAnalyst application
ui <- dashboardPage(
  #define the skin
  skin = "green",
  
  dashboardHeader(title = "TMA"),
  #close the header
  
  ###########################
  #####-----SIDEBAR-----#####
  ###########################
  dashboardSidebar(sidebarMenu(
    id = "main_tabs",
    #create a home tab to explain the app
    menuItem("Home",
             tabName = "home",
             icon = icon("home"),
             selected = T), #close home tab
    
    #create a RNA-Seq tab
    convertMenuItem(menuItem("RNA-Seq Meta Analysis",
                             tabName = "rnaseq",
                             icon = icon("list-alt", lib = "glyphicon"),
                             # upload button for counts
                             fileInput(inputId = "rnaseq_counts",
                                       label = "Upload the raw counts file (tab-delimited)",
                                       accept = c("text/csv", "text/comma-separated-values,text/plain", "txt")),
                             #upload button for experimental design
                             fileInput(inputId = "rnaseq_design",
                                       label = "Upload the experimental design file (tab-delimited)",
                                       accept = c("text/csv", "text/comma-separated-values,text/plain", ".txt")),
                             #example files
                             downloadLink("example_counts", label = "Download example counts file \n"),
                             p(),
                             downloadLink("example_design", label = "Download example design file"),
                             #add button to start analysis
                             actionButton("rnaseq_analyze", "Start Analysis!", icon = icon("laptop"))
                             ),
                    tabName = "rnaseq"), #close rna-seq tab

    #create a microarray tab
    menuItem("Microarray Meta-Analysis",
             tabName = "microarray",
             icon = icon("table-cells")),
    
    #create a tab for demo analysis
    menuItem("Run Demo Analysis",
             tabName = "demo",
             icon = icon("person-running")),
    
    menuItem("Methods Guide",
             tabName = "methods",
             icon = icon("compass")),

    #create a tabe for user manual
    menuItem("User Manual",
             tabName = "manual",
             icon = icon("book"))
  )),  #close the sidebar
  
  ########################
  #####-----BODY-----#####
  ########################
  dashboardBody(
    useShinyjs(),
    
    tabItems(
      ##### HOME TAB #####
      ####################
      # add box explaining the application
      tabItem(tabName = "home", class = "active",
        fluidRow(
          box(title = "Quick overview",
            h3(strong("TranscriptMetaAnalyst"), ": A graphical user interface to meta-analyze RNA-Seq and microarray transcriptomics data"),
            p("This app is designed to help biologists with less programming skills to perform a meta-analysis of several transcriptomic datasets (RNA-Seq or microarray), generated in their own studies or available online, to get robust results and new insights into the data."),
            p("TMA does a gene-level meta-analysis of bulk RNA-Seq (raw counts per gene) and single-colour microarray data. Transcript-level, isoform and splicing analyses are out of scope. Each study needs at least 2 samples per condition, and each comparison at least 2 studies; see the Methods Guide tab."),
            h4("If you would like to check the app, please click on", strong("Run Demo Analysis"), "tab..."),
            hr(),
            p("As shown on the side panel, the app consists of several tabs as follows:"),
            tags$ul(
              tags$li(strong("Home:"), "contains the preliminary information about the app and information about how to use it."),
              tags$li(strong("RNA-Seq Meta Analysis:"), "is designed to perform meta-analysis for a group of RNA-Seq-generated transcriptomic datasets."),
              tags$li(strong("Microarray Meta-Analysis:"), "performs the meta-analysis for normalised single-colour microarray datasets, uploaded as files or imported from GEO, and has its own demo data."),
              tags$li(strong("Run Demo Analysis:"), "could be used to perform a sample RNA-Seq meta analysis on a dataset provided with the app. When the user clicks this tab, the analysis will start automatically and the results will appear on the screen after a few seconds. The provided figures and tables are examples of what the application will provide under RNA-Seq Meta Analysis tab using the user's uploaded files"),
              tags$li(strong("Methods Guide:"), "explains when to use which differential expression, meta-analysis, batch and enrichment method, lists the study requirements, and gives the references."),
              tags$li(strong("User Manual:"), "contains a detailed user manual to guide the user during the app usage and explain all the files required to perform the analysis. It also contains information on how and/or where to get help, and how to cite the application.")
            ),
            solidHeader = T,
            width = 12,
            status = "primary")
        )
      ), #close home tab item
      
      ##### RNA-SEQ ANALYSIS TAB #####
      ################################
      tabItem(tabName = "rnaseq",
        fluidRow(mod_settings_ui("rnaseq_settings", per_study = TRUE)),
        fluidRow(mod_batch_check_ui("rnaseq_batch_check")),
        div(id = "rnaseq_steps",
            fluidPage(
              box(title = "RNA-Seq Quick Start",
                  p("To start the RNA-Seq meta analysis please perform the following steps:",
                    tags$ul(
                      tags$li("After cleaning the raw RNA-Seq read files (either downloaded from public databases of generated via your own experiments) of each sample, align them individually to the studied organism's reference genome."),
                      tags$li("Count the mapped reads to the reference genomic features of interest and combine all the counts in one file. Keeping genomic features as rows and the samples that will be included in the analysis as columns. To produce the final raw (unnormalized) counts comma-separated file.", strong("Note that the counts file should have a 'geneid' column.")),
                      tags$li("Upload the ", strong("counts"),"file produced by any feature counting tool, e.g., featureCounts or htseq-count, as explained above to the raw counts placeholder.", strong(" Please refer to example counts file via the provided download link.")),
                      tags$li("Upload the experimental design file.", strong(" Please refer to the example file via the provided download link."), "This file should contain at least 3 columns, namely 'sample' refers to the sample names used as columns in the counts file, 'study' refers to the study to which the sample belong, 'condition' refers to the treatment or control conditions. With more than two conditions, each one is compared with the reference condition chosen in the Analysis settings. Every condition needs at least 2 samples in a study, and every comparison at least 2 studies (an optional 'batch' column can give the batch of each sample)."),
                      tags$li("Hit",strong('Start Analysis!'), "button and wait for the results to appear.")
                    )),
                  solidHeader = T, width = 12, status = "info")
            )), #close steps box
        #results of the analysis, shown once it has run
        uiOutput("rnaseq_results")
      ), #close rnaseq tab item
      ##### MICROARRAY ANALYSIS TAB #####
      ###################################
      tabItem(tabName = "microarray",
        fluidRow(mod_settings_ui("array_settings", per_study = TRUE, platform = "array")),
        fluidRow(box(
          title = "Microarray data", status = "primary", solidHeader = TRUE, width = 12,
          p("Meta-analysis of normalised microarray expression data. Load the demo data, upload one expression file per study with a design file, or import series from GEO."),
          actionButton("array_demo", "Load demo data", icon = icon("person-running")),
          hr(),
          radioButtons("array_source", "Data source", inline = TRUE,
                       choices = c("Upload expression files" = "upload", "Import from GEO" = "geo")),
          conditionalPanel(
            "input.array_source == 'upload'",
            fluidRow(
              column(6, fileInput("array_files", "Expression files (tab-delimited, one per study)", multiple = TRUE,
                                  accept = c("text/plain", ".txt", ".tsv", ".gz"))),
              column(6, fileInput("array_design", "Experimental design file (tab-delimited)",
                                  accept = c("text/plain", ".txt", ".tsv")))
            ),
            helpText("Each expression file holds one study: the first column has the probe (or gene) IDs, an optional column named 'Symbol' (or 'Gene') gives the gene of each probe, and the other columns are the samples, named as in the 'sample' column of the design file. The design file has the columns 'sample', 'study' and 'condition' (and optionally 'batch'). Values are log2-transformed automatically if they are not log2 already; probes are collapsed to genes by the highest mean expression.")
          ),
          conditionalPanel(
            "input.array_source == 'geo'",
            fluidRow(
              column(8, textInput("geo_ids", "GEO series accessions, separated by commas", placeholder = "e.g. GSE82120, GSE63352")),
              column(4, br(), actionButton("geo_fetch", "Fetch from GEO", icon = icon("download")))
            ),
            uiOutput("geo_config")
          ),
          actionButton("array_analyze", "Start Analysis!", icon = icon("laptop"))
        )),
        uiOutput("array_preparation"),
        uiOutput("array_results")
      ), #close microarray tab item
      ##### DEMO ANALYSIS TAB #####
      #############################
      tabItem(tabName = "demo",
              fluidRow(mod_settings_ui("demo_settings")),
              fluidRow(mod_batch_check_ui("demo_batch_check")),
              fluidRow(box(title = "Reproducible R script", status = "primary", solidHeader = TRUE,
                           width = 12, collapsible = TRUE,
                           downloadButton("demo_script", "Download R script"),
                           helpText("An R script that repeats the demo analysis with the current settings, using the demo data bundled with the package. For your own data, use the button in the Comparisons box of the RNA-Seq tab; its input file paths are placeholders that you replace with the paths of your files."))),
              mod_results_ui("demo")
      ), #close demo tab item
      
      ##### METHODS GUIDE TAB #####
      tabItem(tabName = "methods", methods_guide_ui()),

      ##### USER GUIDE TAB #####
      ##########################
      tabItem(tabName = "manual",
              fluidPage(
                box(title = "Detailed User Guide",
                    solidHeader = T, width = 12, status = "danger",
                    h3("TranscriptMetaAnalyst: Datailed User Manual"),
                    span("A detailed user manual could be downloaded via this"), ##### ADD LINK OF PDF FILE
                    tags$a(href = "https://drive.google.com/file/d/1Ro-DhYj4ctxnVgNov2ZchJbeZYx_cUR_/view?usp=sharing",
                           "file",
                           target = "_blank"),
                    hr(),
                    h3("Contact us"),
                    span("For inquiries and/or feedback, please contact me via"),
                    a("email",
                      href = "mailto:e.abdelsalam@lmu.de"),
                    hr(),
                    h3("How to cite?"),
                    p("Please cite TranscriptMetaAnalyst as follows:")
                    )
              )
      ) #close guide tab item
    ), #close all tab items
    ##################################################################################################################
    
    #add footer and version
    tags$footer(
      tags$p("\U00A9 Developed and supported by Eslam M. Abdel-Salam, Ludwig Maximilians University of Munich (LMU)"),
      align = "right"), # Dasbboardbody close 
  ) #close the body
) #close the whole page
