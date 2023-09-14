
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
  
  
  
)