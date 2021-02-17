# --------------------------------------------------------------------------
# Single Cell Visualisation App
# Ben Crabtree, 2021
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Load Libraries
# --------------------------------------------------------------------------

library(shiny)
library(dplyr)
library(heatmaply)
library(pheatmap)
library(RColorBrewer)
library(BBmisc)

# --------------------------------------------------------------------------
# Necessary to get plots to load with pheatmap
# (DON'T NEED WITH HEATMAPLY)
# --------------------------------------------------------------------------

graphics.off()

# --------------------------------------------------------------------------
# Increase max file upload size to 150mb 
# --------------------------------------------------------------------------

options(shiny.maxRequestSize = 150*1024^2)

# --------------------------------------------------------------------------
# Load Single Cell Data
# --------------------------------------------------------------------------

sce_kalluri = readRDS(file = "KalluriDS.rds")
sce_miao = readRDS(file = "MiaoDS.rds")
sce_tms = readRDS(file = "TMSDS.rds")

# --------------------------------------------------------------------------
# User Interface
# --------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel(title = "Single Cell Visualiser"),
  titlePanel(title = span("Ben Crabtree, 2021", 
                                style = "font-size: 14px")),
  #titlePanel("Ben Crabtree, 2021"),
  sidebarLayout(
    sidebarPanel(
      # Input: Select dataset(s) to visualise
      checkboxGroupInput(inputId = "datasets", label = "Select Dataset(s):",
                         choices = c("Kalluri", "Miao", "TMS"),
                         selected = "Kalluri"),
      # Input: Upload a user dataset
      fileInput(inputId = "dsupload", label = "Upload Dataset: ",
                #multiple = TRUE,
                accept = c(".RData")),
      # Input: Select whether to display heatmap with
      # datasets or subjects as rows
      selectizeInput(inputId = "displayby", label = "Display By:",
                     choices = c("Dataset", "Subject")),
      # Input: Select gene to visualise
      selectizeInput(inputId = "gene", label = "Select Gene:",
                  choices = sort(rownames(sce_kalluri)))
    ),
    
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  tabPanel("Static Heatmap", plotOutput(outputId = 'heatmap_static')),
                  tabPanel("Interactive Heatmap", plotlyOutput('heatmap_interactive'))
      )
      
      # Print debug info
      # verbatimTextOutput("debug")
    )
  )
) # end ui

# --------------------------------------------------------------------------
# Server
# --------------------------------------------------------------------------

server <- function(input, output, session)
{
  # -------------------------------------------------------------
  # UPDATE THIS PART TO ADD NEW PRELOADED DATASET (sce object)
  # -------------------------------------------------------------
  
  # Make all_genes dataframe with genes and which dataset they come from
  miao_genes = data.frame("dataset" = "Miao", "genes" = sort(rownames(sce_miao)))
  kalluri_genes = data.frame("dataset" = "Kalluri", "genes" = sort(rownames(sce_kalluri)))
  tms_genes = data.frame("dataset" = "TMS", "genes" = sort(rownames(sce_tms)))
  
  all_genes = bind_rows(kalluri_genes, miao_genes, tms_genes)

  # MAKE A FULL SCE LIST
  all_sces = list("Kalluri" = sce_kalluri, "Miao" = sce_miao, "TMS" = sce_tms)
  
  # -------------------------------------------------------------
  # End update section
  # -------------------------------------------------------------
  
  # User uploaded sce object
  sce_upload = reactive({
    sessionEnvir <- sys.frame()
    if (!is.null(input$dsupload))
    {
      sce_upload = readRDS(input$dsupload$datapath, sessionEnvir)
      if (is.null(sce_upload$Subject))
      {
        sce_upload$Subject = "upload_subj"
      }
      if (!is.null(sce_upload$predRes_scClassify) 
          & is.null(sce_upload$CellType))
      {
        sce_upload$CellType = sce_upload$predRes_scClassify
      }
      return(sce_upload)
    }
  })
  
  # Dataframe of genes from uploaded dataset
  # (to make gene input choices reactive)
  sce_upload_genes = reactive({
    if (!is.null(input$dsupload))
    {
      sce_upload_genes = data.frame("genes" = sort(rownames(sce_upload())))
      return(sce_upload_genes)
    }
  })
  
  # Make gene options reactive in dropdown menu.
  # When multiple datasets selected, show UNION of lists of genes
  observe({
    updateSelectInput(
      session,
      "gene",
      choices = all_genes %>%
        filter(dataset %in% input$datasets) %>%
        select(genes) %>%
        vctrs::vec_c(sce_upload_genes()) %>%
        unique()
    )
  })
  
  # SUBSET FULL SCE LIST all_sces BASED ON input$datasets
  sces = reactive({
    filtered_sces = all_sces[input$datasets]
    if (is.null(input$dsupload))
    {
      return(filtered_sces)
    }
    else
    {
      filtered_sces[["sce_upload"]] = sce_upload()
      return(filtered_sces)
    }
  })
  
  sce_names = reactive({
    return(names(sces()))
  })
  
  subject_list = reactive({
    subject_list = vector()
    for (name in sce_names())
    {
      subject_list = c(subject_list, sces()[name][[1]]$Subject)
    }
    return(sort(unique(subject_list)))
  })
  
  # For use in colour coding subject rows in heatmap
  subject_df = reactive({
    subject_df = data.frame(matrix("NA", ncol = 1,
                                   nrow = length(subject_list()) ))
    colnames(subject_df) = "Dataset"
    rownames(subject_df) = subject_list()
    for (subject in subject_list())
    {
      for (name in sce_names())
      {
        if (subject %in% unique(sces()[name][[1]]$Subject))
        {
          subject_df[subject,] = name
        }
      }
    }
      return(subject_df)
  })
  
  cell_types = reactive({
    cell_types = vector()
    for (name in sce_names())
    {
      cell_types = c(cell_types, sces()[name][[1]]$CellType)
    }
    return(sort(unique(cell_types)))
  })
  
  # HELPER FUNCTIONS
  
  make_empty_matrix <- function(rows, cols, names)
  {
    heatmap_matrix = matrix(0, nrow = rows, 
                            ncol = cols, 
                            dimnames = names)
    return(heatmap_matrix)
  }
  
  fill_matrix_by_ds <- function(matrix, celltypes, sces)
  {
    # sces is reactive list of sce objects (subset of full list)
    ds_names = names(sces)
    ds_index = 1
    for (sce in sces)
    {
      for (cell in celltypes)
      {
        # Names of cells of this cell type from this subject
        cell_names = colnames(sce)[sce$CellType == cell]
        # Check if cell_names is empty list
        if (length(cell_names) == 0)
        {
          matrix[ds_names[ds_index], cell] = matrix[ds_names[ds_index], cell] + 0
        }
        # Check if gene is in sce dataset
        else if (!(input$gene %in% rownames(sce)))
        {
          matrix[ds_names[ds_index], cell] = matrix[ds_names[ds_index], cell] + 0
        }
        # Mean of expression levels of this group of cells
        # goes into appropriate place in the matrix
        # (added to current value for cumulative effect across datasets)
        else
        {
          matrix[ds_names[ds_index], cell] = matrix[ds_names[ds_index], cell] + mean(assays(sce[input$gene, cell_names])$logcounts)
        }
      } # end for cell in celltypes
      ds_index = ds_index + 1
    } # end for sce in sces
    return(matrix)
  }
  
  fill_matrix_by_subj <- function(matrix, subjlist, celltypes, sces)
  {
    # sces is reactive list of sce objects (subset of full list)
    for (sce in sces)
    {
      for (subj in subjlist)
      {
        for (cell in celltypes)
        {
          # Names of cells of this cell type from this subject
          cell_names = colnames(sce)[sce$Subject == subj & sce$CellType == cell]
          # Check if cell_names is empty list
          if (length(cell_names) == 0)
          {
            matrix[subj, cell] = matrix[subj, cell] + 0
          }
          # Check if gene is in sce dataset
          else if (!(input$gene %in% rownames(sce)))
          {
            matrix[subj, cell] = matrix[subj, cell] + 0
          }
          # Mean of expression levels of this group of cells
          # goes into appropriate place in the matrix
          # (added to current value for cumulative effect across datasets)
          else
          {
            matrix[subj, cell] = matrix[subj, cell] + mean(assays(sce[input$gene, cell_names])$logcounts)
          }
        } # end for cell in celltypes
      } # end for subj in subjlist
    } # end for sce in sces
    return(matrix)
  }
  
  # CREATE HEATMAP MATRIX
  
  heatmap_matrix = reactive({
    if (input$displayby == "Dataset")
    {
      empty_matrix = make_empty_matrix(length(sce_names()), # Num rows
                                       length(cell_types()), # Num cols
                                       list(sce_names(), cell_types()))
      heatmap_matrix = fill_matrix_by_ds(empty_matrix, cell_types(), sces())
    }
    else if (input$displayby == "Subject")
    {
      empty_matrix = make_empty_matrix(length(subject_list()),
                                       length(cell_types()),
                                       list(subject_list(), cell_types()))
      heatmap_matrix = fill_matrix_by_subj(empty_matrix, subject_list(), 
                                 cell_types(), sces())
    }
    return(heatmap_matrix)
  })

  # Check if there is only one subject
  # Or only one dataset selected
  clustered = reactive({
    if ( ( input$displayby == "Dataset" & length(sce_names()) == 1 )
       | ( input$displayby == "Subject" & length(subject_list()) == 1 ) )
    {
      clustered = FALSE
    }
    else
    {
      clustered = TRUE
    }
    return(clustered)
  })

  row_colours = reactive({
    if (input$displayby == "Dataset")
    {
      #row_colours = as.factor(sce_names())
      row_colours = data.frame(matrix("NA", ncol = 1,
                                      nrow = length(sce_names()) ))
      colnames(row_colours) = "Dataset"
      rownames(row_colours) = sce_names()
      for (name in sce_names())
      {
        row_colours[name,] = name
      }
    }
    else if (input$displayby == "Subject")
    {
      row_colours = subject_df()
    }
    return(row_colours)
  })
  
  # PLOT HEATMAP
  
  output$heatmap_static <- renderPlot({
    pheatmap(mat = heatmap_matrix(), cluster_rows = clustered(),
             annotation_row = row_colours(),
             color = colorRampPalette(c("white", 
                                        brewer.pal(n = 7, name = "Reds")))(100))
  })
  #outputOptions(output, "heatmap_static", suspendWhenHidden=FALSE)
  
  # Interactive Heatmap
  output$heatmap_interactive <- renderPlotly({
    heatmaply(heatmap_matrix(), Rowv = clustered(),
              row_side_colors = row_colours(),
              xlab = "Cell Types",
              ylab = input$displayby,
              colors = colorRampPalette(c("white", 
                                          brewer.pal(n = 7, name = "Reds")))(100))
    #heat.colors(40))
  })
  #outputOptions(output, "heatmap_interactive", suspendWhenHidden=FALSE)
  
  
  # Debug info
  output$debug <- renderPrint({
    
    #for (name in sce_names())
    #{
    #  print(unique(sces()[name][[1]]$Subject))
    #}

  })
 
} # end server

# --------------------------------------------------------------------------
# Run App
# --------------------------------------------------------------------------

shinyApp(ui = ui, server = server)

