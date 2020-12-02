# load packages -----------------------------------------
library(shiny)
library(rstudioapi)
library(fgsea)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DT)

# load data ---------------------------------------------
setwd(dirname(getActiveDocumentContext()$path))

## load Enrichr gene list data
load("data/LINCS_L1000_Chem_Pert.Rda", verbose = TRUE)

## load MOA dataframe
df_new_moa <- read.csv("data/new_moa.csv",
                       sep = ",", 
                       header = TRUE, 
                       stringsAsFactors = FALSE)

## allow for up to 20 MB input file to be uploaded
options(shiny.maxRequestSize = 20*1024^2)

# load data ---------------------------------------------
ui <- fluidPage(
  # App title ----
  titlePanel("Connectivity Mapping Application"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("file1", "Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),

      # Horizontal line ----
      tags$hr(),
      
      # Input: Checkbox if file has header ----
      checkboxInput("header", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      # Input: Select quotes ----
      radioButtons("quote", "Quote",
                   choices = c(None = "",
                               "Double " = '"',
                               "Single Quote" = "'"),
                   selected = '"'),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Select number of rows to display ----
      # radioButtons("disp", "Display",
      #              choices = c(Head = "head",
      #                          All = "all"),
      #              selected = "head")
      
      downloadButton("downloadData", "Download"),
      width = 2),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      tabsetPanel(
        tabPanel("Plot", plotOutput("plot")), 
        tabPanel("Summary", DT::dataTableOutput("summary")), 
        tabPanel("Compounds", DT::dataTableOutput("compound")), 
        tabPanel("MOAs", DT::dataTableOutput("moa"))
      )
      
    )
    
  )
)


# global functions --------------------------------------

## run FGSEA function on rank ordered list
get_gsea <- function(input, seed=123, perms=1000){
  ## run GSEA on up-regulated gene list
  set.seed(seed)
  fgsea_up <- fgsea(pathways=up_gene_matrix_rev.list, stats=input, nperm=perms, nproc=1)
  
  ## run GSEA on down-regulated gene list
  set.seed(seed)
  fgsea_down <- fgsea(pathways=down_gene_matrix_rev.list, stats=input, nperm=perms, nproc=1)
  
  ## calculate weighted connectivity scores
  fgsea_es_down <- fgsea_down[, c("pathway", "ES", "pval", "padj")]
  colnames(fgsea_es_down) <- c("pathway", "ES_down", "pval_down", "padj_down")
  
  fgsea_es_up <- fgsea_up[, c("pathway", "ES", "pval", "padj")]
  colnames(fgsea_es_up) <- c("pathway", "ES_up", "pval_up", "padj_up")
  
  fgsea_es <- merge(fgsea_es_up, fgsea_es_down, by = "pathway")
  
  fgsea_es$w <- ifelse(sign(fgsea_es$ES_up) == sign(fgsea_es$ES_down), 0,
                       (fgsea_es$ES_up - fgsea_es$ES_down)/2)
  fgsea_es <- fgsea_es[order(fgsea_es$w, decreasing = TRUE), ]
  
  ## add drug name columns extracted from pathway
  fgsea_es$Drug <- unlist(lapply(fgsea_es$pathway, function(x) strsplit(x, "_")[[1]][1]))
  fgsea_es$Drug <- tolower(fgsea_es$Drug)
  
  ## add cell line extracted from pathway
  fgsea_es$Cell <- unlist(lapply(fgsea_es$pathway, function(x) strsplit(x, "_")[[1]][3]))
  
  return(fgsea_es)
}

create_top_down <- function(input, padj = 0.25, pval = 0.05){
  top_down <- input %>% 
    filter(w < 0, padj_down <= padj, pval_down <= pval) %>% 
    as.data.frame()
  
  top_down <- merge(top_down, df_new_moa[, c("Drugs", "moa_final")], 
                    by.x = "Drug", by.y = "Drugs", all.x = TRUE)
  colnames(top_down)[11] <- "MOA"
  
  top_down$MOA <- unlist(lapply(top_down$MOA, function(x) strsplit(x, ",")[[1]][1]))
  top_down <- top_down[order(top_down$w), ]
  return(top_down)
}

## create scatterplots
create_plot <- function(input, max=300){
  if (dim(input)[1] > 300){
    input_rev <- input[1:300, ]
  } else {
    input_rev <- input
  }
  table1 <- as.data.frame(table(input_rev$MOA))
  ggplot(input_rev, aes(ES_down, -(ES_up))) +
    geom_point(aes(col = MOA, size = -log10(padj_down))) +
    geom_text_repel(size=2, aes(label = toupper(Drug)), point.padding = unit(0.3, "lines")) +
    labs(title = "Relative size indicates -log(adjusted p-value)",
         subtitle = paste(dim(input)[1], "significant perturbagens by weighted connectivity scores, top 300 shown", sep = " ")) +
    ylab("Absolute Enrichment Score from Up-Regulated List") +
    xlab("Enrichment Score from Down-Regulated List") +
    guides(size = FALSE) +
    scale_color_discrete(name="Primary MOA",
                         labels = paste(table1$Var1, table1$Freq, sep = ": ")) +
    theme(text = element_text(size=12),
          #legend.text = element_text(size = 2),
          #legend.title = element_text(size = 2),
          legend.position = "none")
}


# define server functions -------------------------------
server <- function(input, output) {
  
  ## turn dataframe input into rank ordered list
  rank_order <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    
    df$rank <- Map(function(num1, num2) sign(num1)*-log10(num2), 
                   df[, 2],
                   df[, 3]) %>% unlist()
    df <- df[order(df$rank, decreasing = TRUE), ]
    
    rank_order <- as.numeric(df$rank)
    names(rank_order) <- df[, 1]
    
    return(rank_order)
  })
  
  ## run FGSEA algorithm
  top_down <- reactive({
    fgsea_es <- get_gsea(input = rank_order())
    top_down <- create_top_down(input = fgsea_es)
    return(top_down)
  })
  
  output$plot <- renderPlot({
    if (is.null(top_down())) {return()}
    create_plot(top_down())
  })  
  
  output$summary <- DT::renderDataTable({
    if (is.null(top_down)) {return()}
    df <- top_down()
    df$Drug <- toupper(df$Drug)
    drops <- c('pathway', 'pval_up', 'pval_down')
    df <- df[, !(names(df) %in% drops)]
    datatable(df, 
              rownames = FALSE,
              options = list(
                scrollX = TRUE,
                autoWidth = TRUE, pageLength = 50,
                columnDefs = list(list(width = '200px', targets = "_all")))) %>%
      formatRound(c(2:6), 5) %>% 
      formatStyle(columns = c(2:6), 'text-align' = 'center')
  },
)
  
  output$compound <- DT::renderDataTable({
    if (is.null(top_down())) {return()}
    df <- as.data.frame(table(top_down()$Drug))
    colnames(df) <- c("Drug", "Frequency")
    df <- df[order(df$Frequency, decreasing = TRUE), ]
    df$Drug <- toupper(df$Drug)
    datatable(df, rownames = FALSE, options = list(pageLength = 50))
  })
  
  output$moa <- DT::renderDataTable({
    if (is.null(top_down())) {return()}
    df <- as.data.frame(table(top_down()$MOA))
    colnames(df) <- c("MOA", "Frequency")
    df <- df[order(df$Frequency, decreasing = TRUE), ]
    datatable(df, rownames = FALSE, options = list(pageLength = 50))
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = 'cmap_output.csv',
    content = function(file) {
      write.csv(top_down(), file, row.names = FALSE)
    }
  )
}


# create Shiny app object -------------------------------
shinyApp(ui = ui, server = server)