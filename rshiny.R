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
      fileInput("file_index", "Choose FCS File",
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
    req(input$file_index)
    
    # Read FCS file
    fr <- read.FCS(input$file_index$datapath)
    
    # Using ggcyto for basic visualization
    ref_fcs_create("vi201130.062")
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
      fileInput("file_index", "Choose FCS File",
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
    req(input$file_index)
    
    # Read and transform FCS file
    fs <- read_transform_fs_bv(input$file_index$datapath)
    
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
      fileInput("file_index", "Choose FCS File",
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
    req(input$file_index)
    
    # Read and transform FCS file
    fs <- read_transform_fs_bv(input$file_index$datapath)
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

#ATTEMPT 4####

library(shiny)

ui <- fluidPage(
  $(document).on('shiny:connected', function() {
    var el = document.getElementById('plot');
    var ctx = el.getContext('2d');
    var rect = {};
    var drag = false;
    
    function init() {
      el.addEventListener('mousedown', mouseDown, false);
      el.addEventListener('mouseup', mouseUp, false);
      el.addEventListener('mousemove', mouseMove, false);
    }
    
    function mouseDown(e) {
      rect.startX = e.pageX - this.offsetLeft;
      rect.startY = e.pageY - this.offsetTop;
      drag = true;
    }
    
    function mouseUp() { drag = false; }
    
    function mouseMove(e) {
      if (drag) {
        rect.w = (e.pageX - this.offsetLeft) - rect.startX;
        rect.h = (e.pageY - this.offsetTop) - rect.startY ;
        ctx.clearRect(0,0,el.width,el.height);
        draw();
      }
    }
    
    function draw() {
      ctx.strokeRect(rect.startX, rect.startY, rect.w, rect.h);
    }
    
    init();
  });
  
  tags$script(src = "drawRect.js"), # Include the JavaScript
  plotOutput("plot", click = "plot_click", dblclick = "plot_dblclick"),
  verbatimTextOutput("coords")
)

server <- function(input, output, session) {
  
  startCoord <- reactiveVal(NULL)
  
  observeEvent(input$plot_click, {
    if (is.null(startCoord())) {
      startCoord(list(x = input$plot_click$x, y = input$plot_click$y))
    } else {
      stopCoord <- list(x = input$plot_click$x, y = input$plot_click$y)
      output$coords <- renderText({
        paste0("Start: (", startCoord()$x, ", ", startCoord()$y, ") ",
               "End: (", stopCoord$x, ", ", stopCoord$y, ")")
      })
      startCoord(NULL) # Reset for next rectangle
    }
  })
  
  output$plot <- renderPlot({
    plot(1:10)
  })
  
}

shinyApp(ui = ui, server = server)

#ATTEMPT5####
library(shiny)
library(plotly)

library(shiny)
library(plotly)

ui <- fluidPage(
  title = "Extract Rectangle Coordinates",

  # Inline JavaScript
  tags$script(HTML("
    $(document).ready(function(){
      var startX, startY, endX, endY;
      var isDrawing = false;
      var plotElem = document.getElementById('distPlot');

      plotElem.addEventListener('mousedown', function(e){
        isDrawing = true;
        startX = e.offsetX;
        startY = e.offsetY;
      });

      plotElem.addEventListener('mousemove', function(e){
        if (!isDrawing) return;
        endX = e.offsetX;
        endY = e.offsetY;
      });

      plotElem.addEventListener('mouseup', function(e){
        if(isDrawing) {
          isDrawing = false;

          var rectData = {
            startX: startX,
            startY: startY,
            endX: endX,
            endY: endY
          };

          Shiny.setInputValue('rect_data', rectData);
        }
      });
    });
  ")),
  
  

  # Display the plot
  plotlyOutput("distPlot"),

  # Display rectangle coordinates
  verbatimTextOutput("rect_coords")
)

server <- function(input, output, session) {
  
  # Generate a plot
  output$distPlot <- renderPlotly({
    x <- rnorm(100)
    y <- rnorm(100)
    plot_ly(x = x, y = y, type = "scatter", mode = "markers")
  })

  # Display rectangle coordinates
  output$rect_coords <- renderPrint({
    input$rect_data
  })
}

shinyApp(ui, server)

####

library(shiny)
library(plotly)

ui <- fluidPage(
  titlePanel("Draw a Rectangle on a Normal Distribution"),
  mainPanel(
    plotlyOutput("distPlot"),
    verbatimTextOutput("rectangleCoords")
  )
)

server <- function(input, output, session) {
  # Generate a normal distribution data
  x <- seq(-5, 5, by=0.1)
  y <- dnorm(x)
  
  output$distPlot <- renderPlotly({
    plot_ly(x = ~x, y = ~y, type = 'scatter', mode = 'lines') %>%
      layout(dragmode = "select")
  })
 
  
  observe({
    selected_data <- event_data("plotly_selected")
    if (!is.null(selected_data)) {
      output$rectangleCoords <- renderText({
        paste(
          "X-range: [", selected_data$xrange[1], ", ", selected_data$xrange[2], "]",
          "Y-range: [", selected_data$yrange[1], ", ", selected_data$yrange[2], "]"
        )
      })
    }
  })
}


shinyApp(ui, server)

#ATTEMPT6 ####
library(shiny)

rect_names <- c('Total Bacteria', 'HNA Bacteria', 'LNA Bacteria', 'Total Viruses', 'V1 Viruses', 'V2 Viruses', 'V3 Viruses')
rect_colors <- c('red', 'blue', 'green', 'purple', 'orange', 'cyan', 'pink')

ui <- fluidPage(
  titlePanel("Draw Individual Rectangles on a Plot"),
  
  sidebarLayout(
    sidebarPanel(
      lapply(1:length(rect_names), function(i) {
        tagList(
          h4(rect_names[i]),
          div(style = "font-size: 80%;", # Further reducing the size
              numericInput(paste0("x1_", i), "Top-left X:", 1, width = "45%"),
              numericInput(paste0("y1_", i), "Top-left Y:", 1, width = "45%"),
              numericInput(paste0("x2_", i), "Bottom-right X:", 5, width = "45%"),
              numericInput(paste0("y2_", i), "Bottom-right Y:", 5, width = "45%"),
              actionButton(paste0("drawBtn_", i), "Draw"),
              actionButton(paste0("resetBtn_", i), "Reset")
          ),
          br()
        )
      })
    ),
    mainPanel(
      plotOutput("plotRectangles")
    )
  )
)

server <- function(input, output) {
  
  rect_list <- reactiveVal(vector("list", length(rect_names)))
  
  lapply(1:length(rect_names), function(i) {
    observeEvent(input[[paste0("drawBtn_", i)]], {
      rects <- rect_list()
      rects[[i]] <- list(
        x1 = input[[paste0("x1_", i)]],
        y1 = input[[paste0("y1_", i)]],
        x2 = input[[paste0("x2_", i)]],
        y2 = input[[paste0("y2_", i)]],
        col = rect_colors[i]
      )
      rect_list(rects)
    })
    
    observeEvent(input[[paste0("resetBtn_", i)]], {
      rects <- rect_list()
      rects[[i]] <- NULL
      rect_list(rects)
    })
  })
  
  output$plotRectangles <- renderPlotly({
    p <- plot_ly(x = c(0, 10), y = c(0, 10), type = "scatter", mode = "lines", line = list(color = "#FFFFFF")) %>% 
      layout(showlegend = FALSE)
    
    rects <- rect_list()
    for (rect in rects) {
      if (!is.null(rect$x1) && !is.null(rect$y1) && !is.null(rect$x2) && !is.null(rect$y2)) {
        p <- p %>% add_shape(
          type = "rect", 
          x0 = rect$x1, x1 = rect$x2,
          y0 = rect$y1, y1 = rect$y2,
          line = list(color = rect$col, width = 2),
          fillcolor = "transparent")
      }
    }
    return(p)
  })

}

shinyApp(ui, server)

#ATTEMPT 7 ####
library(shiny)
library(ggplot2)

rect_names <- c('Total Bacteria', 'HNA Bacteria', 'LNA Bacteria', 'Total Viruses', 'V1 Viruses', 'V2 Viruses', 'V3 Viruses')
rect_colors <- c('red', 'blue', 'green', 'purple', 'orange', 'cyan', 'pink')

ui <- fluidPage(
  titlePanel("Draw Individual Rectangles on a Plot"),
  
  sidebarLayout(
    sidebarPanel(
      overflowY = "scroll",
      lapply(1:length(rect_names), function(i) {
        tagList(
          h4(rect_names[i]),
          div(style = "font-size: 80%;", 
              numericInput(paste0("x1_", i), "Top-left X:", 1, width = "45%"),
              numericInput(paste0("y1_", i), "Top-left Y:", 1, width = "45%"),
              numericInput(paste0("x2_", i), "Bottom-right X:", 5, width = "45%"),
              numericInput(paste0("y2_", i), "Bottom-right Y:", 5, width = "45%"),
              actionButton(paste0("drawBtn_", i), "Draw"),
              actionButton(paste0("resetBtn_", i), "Reset")
          ),
          br()
        )
      })
    ),
    mainPanel(
      fixedPanel(
        plotOutput("plotRectangles", width = "100%", height = "auto")
      )
    )
  )
)

server <- function(input, output) {
  
  rect_list <- reactiveVal(vector("list", length(rect_names)))
  
  lapply(1:length(rect_names), function(i) {
    observeEvent(input[[paste0("drawBtn_", i)]], {
      rects <- rect_list()
      rects[[i]] <- list(
        x1 = input[[paste0("x1_", i)]],
        y1 = input[[paste0("y1_", i)]],
        x2 = input[[paste0("x2_", i)]],
        y2 = input[[paste0("y2_", i)]],
        col = rect_colors[i]
      )
      rect_list(rects)
    })
    
    observeEvent(input[[paste0("resetBtn_", i)]], {
      rects <- rect_list()
      rects[[i]] <- NULL
      rect_list(rects)
    })
  })
  
  output$plotRectangles <- renderPlot({
    p <- ggplot(data.frame(), aes(x = c(0, 10), y = c(0, 10))) + 
      xlim(0, 10) + ylim(0, 10) + 
      theme_minimal()
    
    rects <- rect_list()
    for (rect in rects) {
      if (!is.null(rect$x1) && !is.null(rect$y1) && !is.null(rect$x2) && !is.null(rect$y2)) {
        p <- p + geom_rect(xmin = rect$x1, xmax = rect$x2, ymin = rect$y1, ymax = rect$y2, color = rect$col, fill = NA)
      }
    }
    
    return(p)
  })
}

shinyApp(ui, server)


#ATTEMPT 8####

library(shiny)
library(ggplot2)
library(ggcyto)
library(flowCore)
library(tidyverse)

read_transform_fs_bv <- function(x) {
  flowCore::read.flowSet(c(x, x)) %>%
    flowCore::transform(translist_bv)
}

rect_names <- c('Total Bacteria', 'HNA Bacteria', 'LNA Bacteria', 'Total Viruses', 'V1 Viruses', 'V2 Viruses', 'V3 Viruses')
rect_colors <- c('red', 'blue', 'green', 'purple', 'orange', 'cyan', 'pink')

ui <- fluidPage(
  titlePanel("Gating FCM files for bacteria and viruses"),
  
  # Input for selecting sample based on Sample_Name from metadata
  selectInput("selected_sample", "Select Sample", choices = NULL),
  
  sidebarLayout(
    sidebarPanel(
      
      shinyFiles::shinyDirButton("dir", "Choose a directory", "Please select a directory", FALSE),
      
      # Select index of filename based on Sample_Name column
      selectInput("file_index", "Select File:", choices = .GlobalEnv$metadata$Sample_Name),
      
      # Display filename
      textOutput("display_filename"),
      numericInput("bins", "Number of bins:", value = 30, min = 1),
      downloadButton("downloadData", "Download Coordinates"),
      overflowY = "scroll",
      lapply(1:length(rect_names), function(i) {
        tagList(
          h4(rect_names[i]),
          div(style = "font-size: 80%;", 
              numericInput(paste0("x1_", i), "Top-left X:", 1, width = "45%"),
              numericInput(paste0("y1_", i), "Top-left Y:", 1, width = "45%"),
              numericInput(paste0("x2_", i), "Bottom-right X:", 5, width = "45%"),
              numericInput(paste0("y2_", i), "Bottom-right Y:", 5, width = "45%"),
              actionButton(paste0("drawBtn_", i), "Draw"),
              actionButton(paste0("resetBtn_", i), "Reset")
          ),
          br()
        )
      })
    ),
    mainPanel(
      fixedPanel(
        plotOutput("plotRectangles", width = "100%", height = "auto")
      )
    )
  )
)

server <- function(input, output) {
  
  coordinates_df <- reactiveVal(data.frame())
  
  gsbv_fs <- reactive({
    req(input$file_index)
    
    file_path <- file.path(global$datapath, input$file_index)
    
    read_transform_fs_bv(file_path)
  })
  
  rect_list <- reactiveVal(vector("list", length(rect_names)))
  
  lapply(1:length(rect_names), function(i) {
    observeEvent(input[[paste0("drawBtn_", i)]], {
      rects <- rect_list()
      rects[[i]] <- list(
        x1 = input[[paste0("x1_", i)]],
        y1 = input[[paste0("y1_", i)]],
        x2 = input[[paste0("x2_", i)]],
        y2 = input[[paste0("y2_", i)]],
        col = rect_colors[i]
      )
      rect_list(rects)
      
      current_data <- data.frame(
        sample_name = basename(input$file_index$name)
      )
      for (idx in 1:length(rects)) {
        rect = rects[[idx]]
        current_data[paste0(rect_names[idx], "_x1")] = rect$x1
        current_data[paste0(rect_names[idx], "_y1")] = rect$y1
        current_data[paste0(rect_names[idx], "_x2")] = rect$x2
        current_data[paste0(rect_names[idx], "_y2")] = rect$y2
      }
      coordinates_df(current_data)
    })
  })
  
  output$plotRectangles <- renderPlot({
    req(gsbv_fs())
    
    p <- ggcyto::ggcyto(gsbv_fs()[[2]], aes(x = `SSC-H`, y = `FL1-H`), subset = "root") +
      geom_hex(bins = input$bins) +  
      theme_bw() +
      labs(title = basename(input$file_index$name), x = "SSC-H (Side scatter (a.u.))", y = "FL1-H (Green Fluorescence (a.u.))") +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.background = element_rect(colour="white", fill="white"),
            panel.border = element_rect(colour = "white"))
    
    rects <- rect_list()
    for (rect in rects) {
      if (!is.null(rect$x1) && !is.null(rect$y1) && !is.null(rect$x2) && !is.null(rect$y2)) {
        p <- p + geom_rect(xmin = rect$x1, xmax = rect$x2, ymin = rect$y1, ymax = rect$y2, color = rect$col, fill = NA)
      }
    }
    
    return(p)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("coordinates-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(coordinates_df(), file, row.names = FALSE)
    }
  )
  
  
  
  shinyFiles::shinyDirChoose(
    input,
    'dir',
    roots = c(home = '~'),
    filetypes = c('', 'fcs', "csv", "bw")
  )
  
  global <- reactiveValues(datapath = paste0(getwd(), "/", project_title, "/data/raw_data/", metadata$Sample_Name[input$file_index])
                           )
  
  dir <- reactive(input$dir)
  
  output$dir <- renderText({
    global$datapath
  })
  
  
  observe({
    # This should only be done after metadata is available.
    updateSelectInput(session, "file_index", choices = metadata$Sample_Name)
  })
  
  observeEvent(input$dir, {
    req(input$dir)  # ensure 'dir' input is provided
    directoryPath <- parseDirPath(c(home = '~'), input$dir)
    
    # Check if directory exists
    if (dir.exists(directoryPath)) {
      setwd(directoryPath)
      global$datapath <- normalizePath(directoryPath)
    } else {
      showNotification("Directory does not exist!", type = "error")
    }
  })
  
}

shinyApp(ui, server)
