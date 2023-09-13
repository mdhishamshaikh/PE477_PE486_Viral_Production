# Make sure you have the necessary libraries installed:
# install.packages("shiny")
# BiocManager::install(c("flowCore", "flowWorkspace", "ggcyto"))

library(shiny)
library(flowCore)
library(flowWorkspace)
library(ggcyto)
library(tidyverse)

ui <- fluidPage(
  titlePanel("Flow Cytometry Gating with ggcyto"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fcs_file", "Choose FCS File",
                multiple = FALSE,
                accept = c(".fcs")),
      actionButton("plot_btn", "Plot FCS Data")
    ),
    mainPanel(
      plotOutput("cyto_plot")
    )
  )
)

server <- function(input, output) {
  observeEvent(input$plot_btn, {
    req(input$fcs_file)
    
    # Read FCS file
    fr <- read.FCS(input$fcs_file$datapath)
    
    # Using ggcyto for basic visualization
    p <- ggcyto(read_transform_fs_bv(fr), aes(y = `FL1-H`, x = `SSC-H`)) + geom_hex(bins = 64) + theme_minimal()
    
    output$cyto_plot <- renderPlot({
      print(p)
    })
  })
}

shinyApp(ui, server)

############

# Make sure you have the necessary libraries installed:
# install.packages("shiny")
# BiocManager::install(c("flowCore", "flowWorkspace", "ggcyto"))
library(shiny)
library(flowCore)
library(flowWorkspace)
library(ggcyto)
library(magrittr)

# Define transformation
detectors <- c("FSC-H", "SSC-H", "FL1-H", "FL2-H", "FL3-H")

translist_bv <- transformList(detectors, logTransform())

read_transform_fs_bv <- function(x) { 
  flowCore::read.flowSet(x) %>% # Note: Removed ref_fcs from this line as it was not defined earlier
    flowCore::transform(translist_bv)
}

ui <- fluidPage(
  titlePanel("Flow Cytometry Gating with ggcyto"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fcs_file", "Choose FCS File",
                multiple = FALSE,
                accept = c(".fcs")),
      numericInput("bins", "Number of Bins:", value = 100, min = 10, max = 1000),
      actionButton("plot_btn", "Plot FCS Data")
    ),
    mainPanel(
      plotOutput("cyto_plot")
    )
  )
)

server <- function(input, output) {
  observeEvent(input$plot_btn, {
    req(input$fcs_file)
    
    # Read and transform FCS file
    fs <- read_transform_fs_bv(input$fcs_file$datapath)
    
    # Use only the latest file for visualization
    fr <- fs[[length(fs)]]
    
    # Using ggcyto for basic visualization
    p <- ggcyto(fr, aes(x = `SSC-H`, y = `FL1-H`)) + 
      geom_hex(bins = input$bins) +
      geom_density2d(aes(x = `SSC-H`, y = `FL1-H`)) + # Adding density lines
      theme_minimal()
    
    output$cyto_plot <- renderPlot({
      print(p)
    })
  })
}

shinyApp(ui, server)


######

library(shiny)
library(flowCore)
library(flowWorkspace)
library(ggcyto)
library(magrittr)
library(gridExtra)

# ... [rest of your pre-existing code] ...


ui <- fluidPage(
  titlePanel("Flow Cytometry Gating with ggcyto"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("fcs_file", "Choose FCS File",
                multiple = FALSE,
                accept = c(".fcs")),
      numericInput("bins", "Number of Bins:", value = 64, min = 10, max = 1000),
      actionButton("plot_btn", "Plot FCS Data")
    ),
    mainPanel(
      column(8, offset = 2, plotOutput("combined_plot"))
    )
  )
)

server <- function(input, output) {
  observeEvent(input$plot_btn, {
    req(input$fcs_file)
    
    # Read and transform FCS file
    fs <- read_transform_fs_bv(input$fcs_file$datapath)
    fr <- fs[[length(fs)]]
    
    # Main scatter plot
    p <- ggcyto(fr, aes(x = `SSC-H`, y = `FL1-H`)) + 
      geom_hex(bins = input$bins) + 
      theme_minimal()
    
    # X-axis density
    p_x <- ggplot(data = fr@frames[[1]], aes(x = `SSC-H`)) + 
      geom_freqpoly(binwidth = input$bins) +
      theme_minimal() + 
      theme_void()
    
    # Y-axis density
    p_y <- ggplot(data = fr@frames[[1]], aes(x = `FL1-H`)) + 
      geom_freqpoly(binwidth = input$bins) +
      theme_minimal() + 
      theme_void() +
      coord_flip()
    
    # Combine using cowplot
    combined <- cowplot::plot_grid(p_y, NULL, p, p_x, ncol = 2, nrow = 2, rel_widths = c(1, 4), rel_heights = c(1, 4))
    
    output$combined_plot <- renderPlot({
      print(combined)
    })
  })
}

shinyApp(ui, server)

