library(shiny)
library(ggplot2)
library(DT)
library(dplyr)
library(tidyr)

# UI Definition
ui <- fluidPage(
    titlePanel("Differential Expression Analysis Viewer"),
    
    # Add custom CSS for better styling
    tags$head(
        tags$style(HTML("
            .well { background-color: #f8f9fa; }
            .chart-title { margin-bottom: 20px; }
        "))
    ),
    
    sidebarLayout(
        sidebarPanel(
            width = 3,
            # Cutoff inputs
            numericInput("pvalCutoff", 
                        "P-value cutoff",
                        min = 0, 
                        max = 1, 
                        value = 0.05, 
                        step = 0.01),
            
            numericInput("fcCutoff", 
                        "log2 Fold Change cutoff",
                        min = 0, 
                        max = 5, 
                        value = 1, 
                        step = 0.1),
            
            # Color scheme selection
            selectInput("colorScheme", 
                       "Color Scheme",
                       choices = list(
                           "Red/Blue" = "rb",
                           "Orange/Purple" = "op",
                           "Green/Red" = "gr"
                       ),
                       selected = "rb"),
            
            # Point size control
            sliderInput("pointSize", 
                       "Point Size",
                       min = 1, 
                       max = 5, 
                       value = 2, 
                       step = 0.5),
            
            # Download buttons
            downloadButton("downloadPlot", "Download Plot"),
            tags$br(),
            tags$br(),
            downloadButton("downloadData", "Download Results")
        ),
        
        mainPanel(
            width = 9,
            tabsetPanel(
                tabPanel("Volcano Plot",
                    plotOutput("volcanoPlot", 
                              height = "600px",
                              # Add click, hover, and brush
                              click = "plot_click",
                              hover = "plot_hover",
                              brush = "plot_brush"),
                    
                    # Hover info
                    verbatimTextOutput("hover_info"),
                    
                    # Summary statistics
                    tags$div(
                        class = "well",
                        h4("Summary Statistics"),
                        verbatimTextOutput("summary_stats")
                    )
                ),
                
                tabPanel("Results Table",
                    # Table filters
                    fluidRow(
                        column(4,
                            selectInput("geneFilter", 
                                      "Show genes:",
                                      choices = c("All", "Upregulated", "Downregulated"),
                                      selected = "All")
                        ),
                        column(4,
                            numericInput("topN", 
                                       "Show top N genes:",
                                       value = 50,
                                       min = 1,
                                       max = 1000)
                        )
                    ),
                    # Interactive table
                    DTOutput("geneTable")
                )
            )
        )
    )
)

# Server logic
server <- function(input, output, session) {
    # Read and process data
    results <- reactive({
        req(input$pvalCutoff, input$fcCutoff)
        
        data <- read.csv("results/deg_results.csv")
        
        # Add -log10(p-value) column
        data$neglog10pval <- -log10(data$adj.P.Val)
        
        # Add regulation status
        data$regulation <- "NS"
        data$regulation[data$adj.P.Val < input$pvalCutoff & 
                       data$logFC > input$fcCutoff] <- "UP"
        data$regulation[data$adj.P.Val < input$pvalCutoff & 
                       data$logFC < -input$fcCutoff] <- "DOWN"
        data$regulation <- factor(data$regulation, levels = c("UP", "DOWN", "NS"))
        
        return(data)
    })
    
    # Color schemes
    colorSchemes <- list(
        rb = c("UP" = "red", "DOWN" = "blue", "NS" = "grey"),
        op = c("UP" = "orange", "DOWN" = "purple", "NS" = "grey"),
        gr = c("UP" = "green", "DOWN" = "red", "NS" = "grey")
    )
    
    # Volcano plot
    output$volcanoPlot <- renderPlot({
        req(results())
        
        ggplot(results(), aes(x = logFC, 
                             y = neglog10pval,
                             color = regulation,
                             label = gene_id)) +
            geom_point(alpha = 0.6, size = input$pointSize) +
            theme_minimal(base_size = 12) +
            scale_color_manual(values = colorSchemes[[input$colorScheme]]) +
            geom_hline(yintercept = -log10(input$pvalCutoff), 
                      linetype = "dashed", 
                      color = "darkgrey") +
            geom_vline(xintercept = c(-input$fcCutoff, input$fcCutoff), 
                      linetype = "dashed", 
                      color = "darkgrey") +
            labs(title = "Differential Expression Analysis",
                 subtitle = paste0("FDR < ", input$pvalCutoff, 
                                 ", |log2FC| > ", input$fcCutoff),
                 x = "log2 Fold Change",
                 y = "-log10 adjusted p-value",
                 color = "Regulation") +
            theme(
                plot.title = element_text(size = 16, face = "bold"),
                plot.subtitle = element_text(size = 10),
                legend.position = "right",
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 10),
                axis.title = element_text(size = 12),
                axis.text = element_text(size = 10)
            )
    })
    
    # Hover information
    output$hover_info <- renderPrint({
        req(input$plot_hover)
        hover <- input$plot_hover
        
        # Get the point under the hover using the correct column names
        point <- nearPoints(results(), 
                          hover, 
                          xvar = "logFC",
                          yvar = "neglog10pval",
                          threshold = 5, 
                          maxpoints = 1)
        
        if(nrow(point) > 0) {
            cat("Gene ID:", point$gene_id, "\n")
            cat("log2 Fold Change:", round(point$logFC, 3), "\n")
            cat("Adjusted p-value:", format(point$adj.P.Val, scientific = TRUE), "\n")
            cat("Regulation:", as.character(point$regulation))
        }
    })
    
    # Summary statistics
    output$summary_stats <- renderPrint({
        res <- results()
        cat("Total genes analyzed:", nrow(res), "\n")
        cat("Significantly DE genes (FDR <", input$pvalCutoff, 
            ", |log2FC| >", input$fcCutoff, "):", 
            sum(res$regulation != "NS"), "\n")
        cat("  - Upregulated:", sum(res$regulation == "UP"), "\n")
        cat("  - Downregulated:", sum(res$regulation == "DOWN"), "\n")
    })
    
    # Results table
    output$geneTable <- renderDT({
        res <- results()
        
        # Apply gene filter
        if(input$geneFilter == "Upregulated") {
            res <- res[res$regulation == "UP",]
        } else if(input$geneFilter == "Downregulated") {
            res <- res[res$regulation == "DOWN",]
        }
        
        # Sort by adjusted p-value and take top N
        res <- res[order(res$adj.P.Val),][1:min(input$topN, nrow(res)),]
        
        # Format table
        datatable(res,
                  options = list(pageLength = 10),
                  rownames = FALSE) %>%
            formatSignif(columns = c("adj.P.Val"), digits = 3) %>%
            formatRound(columns = c("logFC"), digits = 2)
    })
    
    # Download handlers
    output$downloadPlot <- downloadHandler(
        filename = function() {
            paste0("volcano_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
        },
        content = function(file) {
            ggsave(file, plot = last_plot(), width = 10, height = 8)
        }
    )
    
    output$downloadData <- downloadHandler(
        filename = function() {
            paste0("deg_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
        },
        content = function(file) {
            write.csv(results(), file, row.names = FALSE)
        }
    )
}

# Run the app
shinyApp(ui = ui, server = server)