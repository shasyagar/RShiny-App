library(shiny)
library(tidyverse)
library(colourpicker) 
library(gridExtra)
library(RColorBrewer)
library(gplots)
library(SummarizedExperiment)
library(DT)
library(ggbeeswarm)


ui <- fluidPage(
  
  #--Background color of the page-----------------------------------------------
  tags$head(
    tags$style(
      HTML(
        "
        body {
          background-color: #E9F1F9; /* Lighter Pastel Blue */
        }
        "
      )
    )
  ),
  
  #--Theme for all boxed texts--------------------------------------------------
  tags$head(
    tags$style(HTML("
      .boxed-text {
        border: 2px solid #A1C6D8; /* Softer Pastel Blue Border */
        padding: 10px;
        background-color: #F4FAFC; /* Very Light Pastel Blue Background */
        border-radius: 5px;
        line-height: 1.1;
        margin-bottom: 10px;
        text-align: center;
      }
    "))
  ),
  
  #--Theme for title text-------------------------------------------------------
  tags$head(
    tags$style(HTML("
      .title-text {
        border: 2px solid #A1C6D8;
        padding: 5px;
        background-color: #D0E9F7; /* Soft Pastel Blue */
        border-radius: 5px;
        width: 50%;
        margin: 0 auto;
        text-align: center;
        line-height: 1;
        margin-top: 10px;
        margin-bottom: 10px;
      }
    "))
  ),
  
  #--Theme for underline text---------------------------------------------------
  tags$head(
    tags$style(HTML("
      .underlined-text {
        text-decoration: underline;
      }
    "))
  ),
  
  #--Theme for sidebars---------------------------------------------------------
  tags$style(
    HTML(
      "
      .sidebar-border {
        border: 1px solid #A1C6D8; /* Softer Pastel Blue Border */
        padding: 10px;
        background-color: #F4FAFC; /* Very Light Pastel Blue */
        margin-left: 5px;
      }
      "
    )
  ),
  
  #--Theme for tabs-------------------------------------------------------------
  tags$style(
    HTML(
      "
      .tab-border {
        border: 1px solid #A1C6D8;
        padding: 10px 20px;
        margin-left: -10px;
        margin-right: -17px;
      }
      .tab-border strong {
        font-weight: bold;
      }
      "
    )
  ),
  
  #--Theme for credit box-------------------------------------------------------
  tags$head(
    tags$style(HTML("
      .credit-text {
        border: 2px solid #A1C6D8;
        padding: 15px;
        background-color: #F4FAFC;
        border-radius: 5px;
        line-height: 1;
        margin-top: 10px;
        margin-bottom: 10px;
      }
    "))
  ),
  
  #--Title text-----------------------------------------------------------------
  
  div(class = "boxed-text",
      h4(strong("Huntington's Disease RNA-Seq Data: Insights and Analysis"))
  ),
  
  #--Main Panel-----------------------------------------------------------------
  div(
    style = "border: 2px solid #A1C6D8; border-radius: 5px; background-color: #F9FBFF; height: auto; width: 90%; margin: 0 auto; display: flex; flex-direction: column;",
    tabsetPanel(
      tabPanel(title = tags$span(class = "tab-border",tags$strong("Sample Information")),
               sidebarLayout(
                 sidebarPanel(width = 3,
                              tags$head(
                                tags$style(".btn-file {background-color:#D0E9F7;}.progress-bar{color:black;background-color:#E9F1F9;}")),
                              fileInput(inputId = "sample_file",label = HTML(paste0( icon("upload"),"  Upload Metadata File ")),accept = c(".csv",".txt")),
                              class="sidebar-border"
                 ),
                 
                 #Show a plot of the generated distribution
                 mainPanel(
                   
                   tabsetPanel(
                     tabPanel(title = tags$span(class = "tab-border",tags$strong("Summary")),class="tab-border",
                              div(DT::dataTableOutput("sample_summary"),style="width:100%;")
                     ),
                     tabPanel(title = tags$span(class = "tab-border",tags$strong("Metadata")),class="tab-border",
                              div(DT::dataTableOutput("sample_table"), style = "font-size:80%; width: 100%;")
                     ),
                     tabPanel(title = tags$span(class = "tab-border", tags$strong("Scatter Plot")),
                              sidebarPanel(
                                selectInput(inputId = "sample_x", label = "Change X Variable",
                                            choices = c("Age of death", "Post mortem interval", "mRNA-seq reads", "RNA integrity number"),
                                            selected = "mRNA-seq reads"),
                                selectInput(inputId = "sample_y", label = "Change Y Variable",
                                            choices = c("Age of death", "Post mortem interval", "mRNA-seq reads", "RNA integrity number"),
                                            selected = "RNA integrity number"),
                                div(
                                  submitButton(text = "Plot", icon = icon("chart-line")),
                                  style = "width: 100%; text-align: center;"
                                ),
                                class = "sidebar-border"
                              ),
                              mainPanel(
                                plotOutput("sample_scatter_plot")
                              )
                     ),
                     tabPanel(title = tags$span(class = "tab-border",tags$strong("Violin Plot")),
                              sidebarPanel(
                                selectInput(inputId = "sample_x", label = "Change X Variable",
                                            choices = c("Age of death", "Post mortem interval", "mRNA-seq reads","RNA integrity number"),
                                            selected = "mRNA-seq reads"),
                                selectInput(inputId = "sample_y", label = "Change Y Variable",
                                            choices = c("Age of death", "Post mortem interval", "mRNA-seq reads","RNA integrity number"),
                                            selected = "RNA integrity number"),
                                div(
                                  submitButton(text = "Plot", icon = icon("chart-area")),
                                  style = "width: 100%; text-align: center;"  
                                ),
                                class = "sidebar-border",
                              ),  
                              mainPanel(
                                plotOutput("sample_plot")
                              )
                     )
                   )
                 )
               )
      ),
      
      #--Counts tab-------------------------------------------------------------
      
      tabPanel(title = tags$span(class = "tab-border",tags$strong("Counts Data")),
               sidebarLayout(
                 sidebarPanel( width = 3,
                               tags$head(
                                 tags$style(".btn-file {background-color:#D0E9F7;}.progress-bar{color:black;background-color:#E9F1F9;}")),
                               #input count matrix
                               fileInput(inputId = "counts_file", label = HTML(paste0( icon("upload"),"  Upload Counts Data ")), accept = ".csv"),
                               # Add slider inputs
                               tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #D0E9F7; color: black}")),
                               sliderInput(inputId = "slid_var",label = "Choose the gene percentile of variance", min = 0, max = 100, value = 60, step = 1),
                               tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: #D0E9F7; color: black}")),
                               sliderInput(inputId = "slid_zero",label = "Choose genes which have samples with non-zero values", min = 0, max = 69, value = 55, step = 1 ),
                               class="sidebar-border",
                               div(
                                 submitButton(text = "Filter", icon = icon("filter")),
                                 style = "width: 100%; text-align: center;"  
                               )
                 ),
                 
                 # Show a plot of the generated distribution
                 mainPanel(
                   tabsetPanel(
                     tabPanel(title = tags$span(class = "tab-border",tags$strong("Filtered Data")),class="tab-border",
                              tableOutput("filter_count")
                     ),
                     tabPanel(title = tags$span(class = "tab-border",tags$strong("Scatter Plot")),class="tab-border",
                              plotOutput("count_scatter")
                     ),
                     tabPanel(title = tags$span(class = "tab-border",tags$strong("Heatmap")),class="tab-border",
                              plotOutput("clus_heatmap",width = "80%", height = "500px")
                     ),
                     tabPanel(title = tags$span(class = "tab-border",tags$strong("PCA")),
                              sidebarPanel(
                                tags$style(HTML(".js-irs-2 .irs-single, .js-irs-2 .irs-bar-edge, .js-irs-2 .irs-bar {background: #D0E9F7; color: black}")),
                                sliderInput("top_PC",label = "Choose number of top PCs", min = 2, max = 20, value = 5, step = 1),
                                div(
                                  submitButton(text = "Plot", icon = icon("chart-gantt")),
                                  style = "width: 100%; text-align: center;"  
                                ),
                                class="sidebar-border"
                              ),
                              mainPanel(
                                plotOutput("pca_plot",width = "120%", height = "400px")
                              )
                     )
                   )
                 )
               )
      ),
      
      #--DE tab-----------------------------------------------------------------
      
      tabPanel(title = tags$span(class = "tab-border", tags$strong("Differential Expression")),
               sidebarLayout(
                 sidebarPanel(width = 3,
                              # Input count matrix
                              fileInput(inputId = "deseq_file", 
                                        label = HTML(paste0(icon("upload"), "  Upload Expression File ")), 
                                        accept = c(".csv", ".txt")),
                              class = "sidebar-border"
                 ),
                 # Show a plot of the generated distribution
                 mainPanel(
                   tabsetPanel(
                     tabPanel(title = tags$span(class = "tab-border", tags$strong("Data Table")), class = "tab-border",
                              div(DT::dataTableOutput("DE_summary"), style = "font-size:80%; width: 100%;")
                     ),
                     tabPanel(title = tags$span(class = "tab-border", tags$strong("Plot")),
                              sidebarPanel(
                                radioButtons(inputId = "x_axis", 
                                             label = "Choose X-axis variable", 
                                             choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), 
                                             selected = "log2FoldChange"),
                                radioButtons(inputId = "y_axis", 
                                             label = "Choose Y-axis variable", 
                                             choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), 
                                             selected = "padj"),
                                # Add color inputs
                                colourInput(inputId = "base", label = "Choose base colour", value = "#B0E0E6"),
                                colourInput(inputId = "highlight", label = "Choose highlight colour", value = "#003366"),
                                # Add slider inputs
                                tags$style(HTML(".js-irs-3 .irs-single, .js-irs-3 .irs-bar-edge, .js-irs-3 .irs-bar {background: #D0E9F7; color: black}")),
                                sliderInput(inputId = "padj_slider", 
                                            label = "Choose a padj value as threshold", 
                                            min = -35, max = 0, value = -5, step = 1),
                                # Add a submit button
                                div(
                                  submitButton(text = "Generate Volcano Plot", icon = icon("chart-line")),
                                  style = "width: 100%; text-align: center;"
                                ),
                                class = "sidebar-border"
                              ),
                              mainPanel(
                                tabsetPanel(
                                  tabPanel(title = tags$span(class = "tab-border", tags$strong("Volcano Plot")),
                                           plotOutput("volcano", width = "110%", height = "500px")
                                  ),
                                  tabPanel(title = tags$span(class = "tab-border", tags$strong("Filtered Table")), class = "tab-border",
                                           div(DT::dataTableOutput("volcano_table"), style = "font-size:80%; width: 100%;")
                                  )
                                )
                              )
                     )
                   )
                 )
               )
      )
      ,
      
      #--GSEA tab---------------------------------------------------------------
      tabPanel(title = tags$span(class = "tab-border", tags$strong("GSEA")),
               sidebarLayout(
                 sidebarPanel(width = 3,
                              # Input FGSEA results
                              fileInput(inputId = "fgsea_file",
                                        label = HTML(paste0(icon("upload"), "  Upload FGSEA Results File ")),
                                        accept = ".csv"),
                              class = "sidebar-border"
                 ),
                 # Main layout for plots and tables
                 mainPanel(
                   tabsetPanel(
                     tabPanel(title = tags$span(class = "tab-border", tags$strong("Barplot")),
                              sidebarLayout(
                                sidebarPanel(width = 3,
                                             tags$style(HTML(".js-irs-4 .irs-single, .js-irs-4 .irs-bar-edge, .js-irs-4 .irs-bar {background: #D0E9F7; color: black}")),
                                             sliderInput(inputId = "pth_threshold",
                                                         label = "Choose Padj threshold",
                                                         min = -25, max = 0, value = -18, step = 1),
                                             div(
                                               submitButton(text = "Plot Barplot", icon = icon("chart-bar")),
                                               style = "width: 100%; text-align: center;"
                                             ),
                                             class = "sidebar-border"
                                ),
                                mainPanel(
                                  plotOutput("fgsea_bars", width = "110%", height = "500px")
                                )
                              )
                     ),
                     tabPanel(title = tags$span(class = "tab-border", tags$strong("Filtered Table")),
                              sidebarLayout(
                                sidebarPanel(width = 3,
                                             tags$style(HTML(".js-irs-5 .irs-single, .js-irs-5 .irs-bar-edge, .js-irs-5 .irs-bar {background: #D0E9F7; color: black}")),
                                             sliderInput(inputId = "path_slid",
                                                         label = "Choose Padj threshold",
                                                         min = -15, max = 0, value = -9, step = 1),
                                             radioButtons(inputId = "all_path",
                                                          label = "Choose pathways to display",
                                                          choices = c("Positive", "Negative", "All"),
                                                          selected = "All"),
                                             div(
                                               submitButton(text = "Filter Table", icon = icon("filter")),
                                               style = "width: 100%; text-align: center;"
                                             ),
                                             div(
                                               downloadButton(outputId = "download_fgsea_table", label = "Download"),
                                               style = "width: 100%; text-align: center; margin-top: 10px;"
                                             ),
                                             class = "sidebar-border"
                                ),
                                mainPanel(
                                  div(DT::dataTableOutput("fgsea_filt_table"), style = "font-size:80%; width: 130%;overflow-x: auto;")
                                )
                              )
                     ),
                     tabPanel(title = tags$span(class = "tab-border", tags$strong("Scatter Plot")),
                              sidebarLayout(
                                sidebarPanel(width = 3,
                                             tags$style(HTML(".js-irs-6 .irs-single, .js-irs-6 .irs-bar-edge, .js-irs-6 .irs-bar {background: #D0E9F7; color: black}")),
                                             sliderInput(inputId = "scatter_slid",
                                                         label = "Choose Padj threshold",
                                                         min = -15, max = 0, value = -7, step = 1),
                                             div(
                                               submitButton(text = "Generate Scatter Plot", icon = icon("chart-area")),
                                               style = "width: 100%; text-align: center;"
                                             ),
                                             class = "sidebar-border"
                                ),
                                mainPanel(
                                  plotOutput("NES_scatter", width = "110%", height = "500px")
                                )
                              )
                     )
                   )
                 )
               )
      )
      ,
      
      #--Credits tab------------------------------------------------------------ 
      tabPanel(
        title = tags$span(class = "tab-border", tags$strong("Acknowledgement")),
        div(
          class = "credit-text",
          p(
            "Credit: ",
            tags$a(
              href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4670106/",
              "Labadorf A, Hoss AG, Lagomarsino V, Latourelle JC et al. RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression. PLoS One 2015;10(12):e0143563. PMID: 26636579"
            )
          ),
          p(
            "Dataset: ",
            tags$a(
              href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810",
              "GSE64810"
            )
          ),
          p(
            strong(
              "This application was developed by Sanjana Hasyagar as the final project for the BF591 course. Special thanks to Dr.Joseph Orofino and Dr.Adam Labadorf for their guidance and support throughout this project."
            )
          )
        )
      )
    )
  ),
)

# Define server 
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 40*1024^2) # set max file size to 40 MB
  
  #--Tab 1 sample info--------------------------------------------------------    
  
  #Read in sample data
  sample_data <- reactive({
    req(input$sample_file)
    data <-read.csv(input$sample_file$datapath, sep="\t", header = FALSE, stringsAsFactors = FALSE)%>%as_tibble()
    data_info <- data %>% 
      mutate(V1 = apply(data["V1"],1, function(x) gsub("!", "", x))) %>% 
      as.data.frame()
    data_info<-data_info[c(1:3,6,8:14),]%>%
      mutate(Columns = c("GEO accession","Status","Submission date", "Channel count","Organism source","Tissue source","Diagnosis", "Post mortem interval","Age of death","RNA integrity number","mRNA-seq reads"), .before = 1) %>%
      select(-V1)%>%
      t()
    colnames(data_info) <- data_info[1,]
    data_info<- data_info[-1,]%>%
      apply(2, function(x) gsub("tissue: ", "", x)) %>%
      apply(2, function(x) gsub("diagnosis: ", "", x)) %>%
      apply(2, function(x) gsub("pmi: ", "", x)) %>%
      apply(2, function(x) gsub("age of death: ", "", x)) %>%
      apply(2, function(x) gsub("rin: ", "", x)) %>%
      apply(2, function(x) gsub("mrna-seq reads: ", "", x)) %>%
      as.data.frame() 
    rownames(data_info) <- NULL
    data_info <-mutate(data_info,across(c(4,8:11), as.double))
  })
  
  #--Tab 2 metadata table-------------------------------------------------------  
  
  summary_tablef <- function(data) {
    # Count number of rows and columns
    n_rows <- nrow(data[,-1])
    n_cols <- ncol(data)
    
    # Create a data frame with column information
    col_info <- data.frame(
      "Column Name" = names(data),
      "Type" = sapply(data, class),
      stringsAsFactors = FALSE
    )
    # Replace dots in column names with spaces
    names(col_info) <- gsub("\\.", " ", names(col_info))
    
    # Add mean or distinct values for columns
    for (i in 1:n_cols) {
      if (is.numeric(data[[i]])) {
        col_mean <- mean(data[[i]], na.rm = TRUE)
        col_info$Mean[i] <- round(col_mean, 3)
      } 
      else { col_info$Mean[i] <- "N/A"
      }
    }
    return(col_info%>%as_tibble())
  }
  
  #--Tab 1 scatter plot --------------------------------------------------------    
  
  scatter_plot <- function(df, x_var, y_var) {
    ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]], color = .data[[x_var]])) +
      geom_point(alpha = 0.7, size = 2.5, shape = 21) +
      geom_smooth(method = "lm", color = "blue", se = FALSE) +
      labs(x = x_var, y = y_var, color = x_var) +
      theme_minimal()
  }
  
  #--Tab 1 violin plot----------------------------------------------------------    
  
  #Generate violin plot of sample
  violin <- function(df, x_var, y_var) {
    ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) +
      geom_violin(color = "black", alpha = 0.5) +
      geom_point(stat = "summary", fun = "mean", color = "black", size = 2.5, shape = 21) +
      labs(x = x_var, y = y_var, fill = x_var) +
      theme_minimal()
  }
  
  #--Tab 2 counts input data----------------------------------------------------  
  
  counts_data <- reactive({
    req(input$counts_file)  # Ensure file is uploaded
    data <- read.csv(input$counts_file$datapath)  # Read the CSV without setting row.names
    # Rename the first column to "gene"
    colnames(data)[1] <- "gene"
    # Check if data is loaded correctly (you can see this in the console)
    #print(head(data))  # Print the first few rows to inspect
    
    return(data)
  })
  
  # --Tab 2 filtered table -----------------------------------------------------
  
  filter_table <- function(data, pass_filter1, pass_filter2) {
    count_sum <- data[,-1]
    count_sum$non_zeros <- apply(data[,-1], 1, function(x) sum(x != 0))
    geneVar <- apply(data[,-1], 1, var, na.rm = TRUE)  # Ensure var() handles NA values
    
    # Count the percentile variance
    xPercentile <- quantile(geneVar, pass_filter1/100, na.rm = TRUE)  # Remove NA values during quantile calculation
    
    # Find out which rows pass the filter1
    pass_xPercentile <- which(geneVar >= xPercentile)
    
    # Subset the matrix to let the rows that pass remain
    count_sum <- count_sum[pass_xPercentile, ]
    count_sum$pass_zero <- ifelse(count_sum$non_zeros >= pass_filter2, "Pass", "Fail")
    
    num_passing <- sum(count_sum$pass_zero == "Pass")
    perc_passing <- round(num_passing / nrow(data[,-1]) * 100, 2)
    
    # Create a data frame with filter information
    sum_info <- tibble(
      metric = c("Number of samples", "Number of genes", "Number of genes passing", " % of genes passing", "Number of genes not passing", "% of genes not passing"),
      value = c(ncol(data[,-1]), nrow(data[,-1]), num_passing, paste0(perc_passing, "%"), nrow(data[,-1]) - num_passing, paste0(100 - perc_passing, "%"))
    )
    
    return(sum_info)
  }
  
  #--Tab 2 scatterplots---------------------------------------------------------
  
  plot_variance_vs_median <- function(data, pass_filter1) {
    # Calculate median and variance
    new_tib <- tibble(count_median = apply(data[,-1], 1, median)) %>%
      add_column(variance = apply(data[,-1], 1, var), rank = rank(.))
    
    # Compute log of variance
    new_tib$log_variance <- log10(new_tib$variance + 1)  # +1 to avoid log(0)
    
    # Pass filter threshold for variance
    xPercentile <- quantile(new_tib$log_variance, pass_filter1 / 100)
    
    # Determine passing or failing genes based on variance filter
    new_tib <- new_tib %>%
      mutate(determine = ifelse(log_variance >= xPercentile, "Pass", "Fail"))
    
    # Scatter plot: Median count vs Log(Variance)
    ggplot(new_tib, aes(x = count_median, y = log_variance, color = determine)) +
      geom_point(size = 2) +
      scale_color_manual(values = c("Fail" = "#89CFF0", "Pass" = "darkblue")) +  # Light gray for "Fail" and dark blue for "Pass"
      theme_minimal() +
      labs(title = "Median Count vs Log(Variance)", x = "Median Count", y = "Log(Variance)") #+
    #theme(legend.position = "none")  # Hide legend, if desired
  }
  
  plot_nonzero_vs_median <- function(data, pass_filter2) {
    # Calculate median and number of non-zero entries
    new_tib <- tibble(count_median = apply(data[,-1], 1, median)) %>%
      add_column(non_zeros = apply(data[, -1], 1, function(x) sum(x != 0)), rank = rank(.)) %>%
      mutate(determine = ifelse(non_zeros >= pass_filter2, "Pass", "Fail"),
             num_zeros = ncol(data) - non_zeros)
    
    # Scatter plot: Median count vs Number of Zeros
    ggplot(new_tib, aes(x = count_median, y = num_zeros, color = determine)) +
      geom_point(size = 2) +
      scale_color_manual(values = c("Fail" = "#89CFF0", "Pass" = "darkblue")) +  # Light gray for "Fail" and dark blue for "Pass"
      theme_minimal() +
      labs(title = "Median Count vs Number of Zeros", x = "Median Count", y = "Number of Zeros") #+
    #theme(legend.position = "none")  # Hide legend, if desired
  }
  
  
  
  
  #--Tab 2 heatmap--------------------------------------------------------------
  
  #generate filter matrix for heatmap
  filter_res <- function(data, pass_filter1, pass_filter2) {
    filt_res <- data %>% 
      mutate(non_zeros = apply(data[, -1], 1, function(x) sum(x != 0)),
             variance = apply(data[,-1], 1, var))
    xPercentile <- quantile(filt_res$variance, pass_filter1/100)
    filt_res <- filt_res %>% filter(variance > xPercentile & non_zeros > pass_filter2)%>% as.data.frame()
    rownames(filt_res) <- filt_res[,1]
    return (filt_res[,-c(1,71,72)])
  }
  
  #generate count heatmap after filtering
  plot_heatmap <- function(filter_data) {
    # Define custom pastel blue color palette
    custom_palette <- colorRampPalette(c("#1B3B6F", "#4682B4", "#B3CDE0", "#DDEAF6", "#FFFFFF"))(100)
    
    # Transform the data for the heatmap
    num_matrix <- filter_data %>% as.matrix() %>% log2()
    num_matrix[!is.finite(num_matrix)] <- NA  # Replace non-finite values with NA
    
    # Plot the heatmap
    heatmap.2(
      num_matrix, 
      col = custom_palette, 
      trace = "none", 
      xlab = "Samples", 
      ylab = "Genes", 
      margins = c(5, 8),
      key = TRUE, 
      key.title = "Expression level", 
      key.xlab = "Expression"
    )
  }
  
  #--Tab 2 PCA beeswarmplot-----------------------------------------------------
  
  #generate PCA beeswarmplot
  plot_beeswarm <- function(data, N) {
    pca_results <- prcomp(scale(t(data[,-1]%>%as.data.frame())), center=FALSE, scale=FALSE)
    plot_tibble <- as_tibble(pca_results$x) %>%
      add_column(sample = rownames(pca_results$x), .after = 0)
    meta <- tibble(sample = rownames(pca_results$x)) %>%
      mutate(Diagnosis = if_else(row_number() <= 49, "normal", "Huntington's Disease"))
    biplot <- dplyr::left_join(plot_tibble, meta, by = "sample")
    biplot$Diagnosis <- factor(ifelse(biplot$Diagnosis == "normal" & seq_len(nrow(biplot)) <= 49, "normal", "Huntington's Disease"))
    biplot_select <- dplyr::select(biplot, 1:N+1, 71)
    # Define a custom color palette with repeated colors
    all_colors <- rep(c("#023E8A", "#0077B6", "#0096C7", "#00B4D8", "#48CAE4", "#90E0EF", "#ADE8F4", "#CAF0F8", "#005F73", "#0A9396"), 2)
    
    top_var <- head(summary(pca_results)$importance[2, ], N)
    top_var_percent <- 100 * top_var
    # Create a vector of the top N principal components with their percentage contribution
    pcs <- paste0("PC", seq_along(top_var_percent), " (", round(top_var_percent, 2), "%)")
    
    beeswarm_plot <- biplot_select %>%
      pivot_longer(cols = PC1:N, names_to = "PC", values_to = "value") %>%
      ggplot(aes(x = factor(PC, levels = paste0("PC", 1:N)), y = value, color = PC)) +
      geom_beeswarm()+
      geom_beeswarm(shape=21,color="black")+
      scale_color_manual(values = all_colors, labels = pcs) +
      theme_classic() +
      labs(x = "Principal Components", y = "Values")
    beeswarm_plot
  }
  
  #--Tab 3 DE input data--------------------------------------------------------
  
  #Read in DESeq data
  deseq_data <- reactive({
    req(input$deseq_file)
    data <-read.table(input$deseq_file$datapath, sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%
      as_tibble()%>%dplyr::rename(gene= X)
    return(data)
  })
  
  
  #--Tab 3 volcano plot and filter table----------------------------------------
  
  #Generate Volcano plot 
  volcano_plot <- function(dataf, x_axis, y_axis, slider, base, highlight) {
    slider_value <- 10^slider
    
    # Convert columns to numeric and handle NA values
    dataf[[x_axis]] <- as.numeric(dataf[[x_axis]])
    dataf[[y_axis]] <- -log10(as.numeric(dataf[[y_axis]]))
    
    
    # Implement a simple volcano plot using ggplot2
    ggplot(dataf, aes(x = !!rlang::sym(x_axis), y = !!rlang::sym(y_axis), color = color_by_condition)) +
      geom_point(size = 2) +
      labs(
        x = x_axis,
        y = y_axis,
        color = paste("padj <", slider_value)  # Set legend title
      ) +
      scale_color_manual(name=paste("padj <", slider_value),values = c("FALSE" = base, "TRUE" = highlight, "NA" = "grey")) +
      theme_minimal() +
      theme(legend.position = "bottom")
  }
  
  #Generate padj filtered table in DE tab
  draw_table <- function(dataf, slider) {
    slider_value <- 10^slider
    
    # Filter data based on condition
    filtered_data <- dataf
    
    # Create a reactive variable for color by condition
    filtered_data$color_by_condition <- ifelse(filtered_data$padj < slider_value, "TRUE", "FALSE")
    filtered_data$color_by_condition[is.na(filtered_data$padj)] <- "NA"
    
    # Format p-value and p-adjusted columns with more digits
    filtered_data$pvalue <- formatC(filtered_data$pvalue, format = "e", digits = 5)
    filtered_data$padj <- formatC(filtered_data$padj, format = "e", digits = 5)
    
    # Return the filtered data
    return(filtered_data)
  }
  
  #--Tab 4 fgsea input data-----------------------------------------------------
  
  #Read in fgsea data
  fgsea_data <- reactive({
    req(input$fgsea_file)
    data <-read.csv(input$fgsea_file$datapath, sep = '\t', header = TRUE, stringsAsFactors = FALSE)%>%
      as_tibble()
    return(data)
  })
  
  #--Tab 4 filter table---------------------------------------------------------
  
  #Read in fgsea data
  fgsea_data <- reactive({
    req(input$fgsea_file)
    data <-read.csv(input$fgsea_file$datapath, header = TRUE, stringsAsFactors = FALSE)%>%
      as_tibble()
    return(data)
  })
  
  #--Tab 4 filter pathway-------------------------------------------------------
  
  #Generate the filter pathway tables in fgsea
  gsea_table <- reactive({
    filter_df<-fgsea_data()
    filtered_df <- dplyr::filter(filter_df, padj < 10^(input$path_slid))
    if (input$all_path == "All") {
      filtered_df <- filtered_df
    } else if (input$all_path == "Positive") {
      filtered_df <- dplyr::filter(filtered_df, NES > 0)
    } else if (input$all_path == "Negative") {
      filtered_df <- dplyr::filter(filtered_df, NES < 0)
    }
    filtered_df <- filtered_df %>%
      dplyr::mutate(pval = formatC(ifelse(is.na(pval), 0, pval), format = "e"),
                    padj = formatC(ifelse(is.na(padj), 0, padj), format = "e"))
    return(filtered_df)
  })
  
  #--Tab 4 barplot--------------------------------------------------------------
  
  #Generate barplot for top pathways of fgsea
  fgsea_top_pathways <- function(fgsea_results, threshold){
    top_positive_nes <- arrange(fgsea_results, desc(NES)) %>%
      dplyr::filter(padj < 10^(threshold) & NES > 0)
    top_negative_nes <- arrange(fgsea_results, desc(NES)) %>%
      dplyr::filter(padj < 10^(threshold) & NES < 0)
    NES_barplot <- dplyr::bind_rows(top_positive_nes, top_negative_nes)%>%
      ggplot() +
      geom_col(aes(x=reorder(pathway,+NES), y=NES, fill = NES > 0))+
      scale_fill_manual(values =c('TRUE' = '#4682B4', 'FALSE' = '#89CFF0')) +
      theme_minimal() +
      coord_flip()+
      theme(legend.position = "none", axis.text.y = element_text( hjust =1 ,size= 7),axis.title.x = element_text(size = 10))+ 
      labs(title="FGSEA Results for Human MSig C2 Gene Set", x= "",y= "Normalised Enrichment Score (NES)")
    return(NES_barplot)
  }
  
  #--Output Tab 1 summary table-------------------------------------------------
  
  #generate sample summary table
  output$sample_summary <- DT::renderDataTable({
    dataf <- sample_data()
    if(is.null(dataf))
      return(NULL)
    
    samplesum_tab <- summary_tablef(dataf)
    
    datatable(
      samplesum_tab,
      options = list(
        columnDefs = list(
          list(className = "dt-center", targets = "_all"),
          list(className = "dt-right", targets = "_all", render = JS("function(data, type, full, meta) {return '<div>' + data + '</div>';}"))
        ),
        ordering = TRUE,
        scrollX = TRUE  # Enable horizontal scrolling
      )
    )
  })
  
  #--Output Tab 1 metadata------------------------------------------------------
  
  #Generate sample file as table
  output$sample_table <- DT::renderDataTable({
    dataf <- sample_data()
    if(is.null(dataf))
      return(null)
    dataf
  }, options = list(
    columnDefs = list(
      list(className = "dt-center", targets = "_all"),
      list(className = "dt-right", targets = "_all", render = JS("function(data, type, full, meta) {return '<div>' + data + '</div>';}"))
    ),
    ordering = TRUE,
    scrollX = TRUE  # Enable horizontal scrolling
  ))
  
  #--Output Tab 1 scatter plot--------------------------------------------------
  
  output$sample_scatter_plot <- renderPlot({
    dataf <- sample_data()  # Get the dataset
    if (is.null(dataf))      # Check if the dataset is null
      return(NULL)            # If null, stop the plot rendering
    
    scatter_plot_obj <- scatter_plot(dataf, input$sample_x, input$sample_y)
    scatter_plot_obj  # Display the density plot
  })
  
  #--Output Tab 1 violin plot---------------------------------------------------
  
  #output sample violin plot
  output$sample_plot <- renderPlot({
    dataf <- sample_data() 
    if(is.null(dataf))
      return(NULL)
    violin_plot <-violin(dataf, input$sample_x,input$sample_y)
    violin_plot
  })
  
  # --Output Tab 2 filter table ------------------------------------------------
  
  #output count filter table
  output$filter_count <- renderTable({
    dataf <- counts_data()
    if(is.null(dataf))
      return(null)
    result_tab <- filter_table(dataf, input$slid_var ,input$slid_zero)
    trans_result <- pivot_wider(
      data = result_tab,
      names_from = metric,
      values_from = value
    )
    trans_result
  }) 
  
  
  #--Output Tab 2 scatterplot---------------------------------------------------
  
  #output count scatter plot
  output$count_scatter <- renderPlot({
    dataf <- counts_data()
    if(is.null(dataf))
      return(null)
    plot1 <- plot_variance_vs_median(dataf,input$slid_var)
    plot2 <- plot_nonzero_vs_median (dataf,input$slid_zero)
    #ggplotly(plot) # 
    grid.arrange(plot1, plot2, nrow=2)
  }) 
  
  #--Output Tab 2 heatmap-------------------------------------------------------
  
  #Output filtering count heatmap
  output$clus_heatmap <- renderPlot({
    dataf <- counts_data()
    if(is.null(dataf))
      return(null)
    num_matrix <- filter_res(dataf, input$slid_var, input$slid_zero)
    plot_heatmap(num_matrix)
  })
  
  #--Output Tab 2 PCA beeswarmplot----------------------------------------------
  
  #Output count PCA beeswarmplot
  output$pca_plot <- renderPlot({
    dataf <- counts_data()
    if(is.null(dataf))
      return(null)
    beeswarm_plot <- plot_beeswarm(dataf,input$top_PC)
    beeswarm_plot
  }) 
  
  #--Output Tab 3 DE table------------------------------------------------------
  
  #output DESeq summary table
  output$DE_summary <- DT::renderDataTable({
    dataf <- deseq_data()
    if(is.null(dataf))
      return(null)
    dataf
  }, options = list(ordering = TRUE,
                    scrollX = TRUE  # Enable horizontal scrolling)
  )) 
  
  
  filtered_data <- reactive({
    slider_value <- 10^input$padj_slider
    dataf <- deseq_data()
    return(draw_table(dataf, input$padj_slider))
  })
  
  #--Output Tab 3 volcano plot--------------------------------------------------  
  
  #output DE volcano plot
  output$volcano <- renderPlot({
    volcano_plot(
      filtered_data(),
      input$x_axis,      # Reactively updates when x-axis input changes
      input$y_axis,      # Reactively updates when y-axis input changes
      input$padj_slider, # Reactively updates when slider changes
      input$base,        # Reactively updates when base changes
      input$highlight    # Reactively updates when highlight changes
    )
  })
  
  #Output Tab 3 filt table------------------------------------------------------
  
  # output the padj filtered table in DE tab
  output$volcano_table <- DT::renderDataTable({
    input$padj_slider
    dataf <- deseq_data()
    isolate({
      filtered_data <- draw_table(dataf, input$padj_slider)
      # Display only rows where color_by_condition is TRUE
      filtered_data <- filtered_data[filtered_data$color_by_condition == "TRUE", ]
      # Remove the unnecessary column before displaying the table
      filtered_data[, !(colnames(filtered_data) %in% c("color_by_condition"))]
    })
  }, options = list(ordering = TRUE,
                    scrollX = TRUE  # Enable horizontal scrolling)
  )) 
  
  #--Output Tab 4 fgsea bar plot------------------------------------------------
  
  #Output fgsea barplot
  output$fgsea_bars <- renderPlot({
    dataf <- fgsea_data()
    if(is.null(dataf))
      return(null)
    bar_plot <- fgsea_top_pathways(dataf,input$pth_threshold)
    bar_plot
  }) 
  
  #--Output Tab 4 download button-----------------------------------------------
  
  #Output download fgsea table
  output$download_fgsea_table <- downloadHandler(
    filename = function() {
      paste("filtered_fgsea_table", ".csv", sep="")
    },
    content = function(file) {
      write.csv(gsea_table(), file, row.names = FALSE)
    }
  )
  
  #--Output Tab 4 fgsea filt table----------------------------------------------
  
  #Output fgsea table
  observeEvent(input$path_slid, {
    output$fgsea_filt_table <- DT::renderDataTable({
      filter_df <- fgsea_data()
      filtered_df <- dplyr::filter(filter_df, padj < 10^(input$path_slid))
      if (input$all_path == "All") {
        filtered_df <- filtered_df
      } else if (input$all_path == "Positive") {
        filtered_df <- dplyr::filter(filtered_df, NES > 0)
      } else if (input$all_path == "Negative") {
        filtered_df <- dplyr::filter(filtered_df, NES < 0)
      }
      filtered_df <- filtered_df %>%
        dplyr::mutate(pval = formatC(ifelse(is.na(pval), 0, pval), format = "e"),
                      padj = formatC(ifelse(is.na(padj), 0, padj), format = "e"))
      return(DT::datatable(filtered_df, 
                           options = list(
                             ordering = TRUE, 
                             scrollX = TRUE,  # Enable horizontal scroll
                             autoWidth = TRUE,
                             columnDefs = list(
                               list(width = '300px', targets = "_all")  # Forces all columns to fit
                             ))))
    })
  })
  
  
  #--Output Tab 4 fgsea scatterplot---------------------------------------------
  
  #Output fgsea scatter plot
  output$NES_scatter <- renderPlot({
    dataf <- fgsea_data()
    slider<-input$scatter_slid
    df <- dplyr::mutate(dataf, new_padj = -log10(padj)) %>%
      dplyr::mutate(status = dplyr::case_when(padj < 10^(slider) ~ "TRUE",
                                              padj >= 10^(slider) ~ "FALSE"))
    # specify color based on the slider value
    df$colors <- ifelse(df$status == "FALSE", "#B0E0E6", "#003153")
    # plotting scatter plot
    scatter <- ggplot(df, aes(x = NES, y = new_padj, color = colors)) +
      geom_point(size = 1) +
      scale_color_manual(values = c("#003153","#B0E0E6"),
                         labels = c("TRUE", "FALSE")) +
      labs(x = "NES", y = "-log10(padj)",color = paste0( "padj < 10^",slider )) +
      theme_bw()+
      theme(legend.position = "bottom") # move legend to bottom of plot
    return(scatter)
  }) 
}

# Run the application
shinyApp(ui = ui, server = server)