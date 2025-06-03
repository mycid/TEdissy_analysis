#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(shiny)
library(ggplot2)
library(readxl)
library(dplyr)
library(stringr)
library(purrr)

# Source your processing function (make sure it's in the app folder)
source("process_plate_run.R")

ui <- fluidPage(
  titlePanel(
    # Show title, author, and dynamic date
    HTML(paste0(
      "<h1>PE diagnostics</h1>",
      "<p><em>Author: Trevor Eakes</em></p>",
      "<p><em>Date: ", Sys.Date(), "</em></p>"
    ))
  ),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Set working directory, input and output files, and specify blanks (comma separated)."),
      textInput("wd", "Working Directory", value = getwd()),
      textInput("input_file", "Input Excel file", value = "PE_2 13_5_25.xlsx"),
      textInput("output_file", "Output CSV file", value = "tidy_final_PE_2_13_5_25.csv"),
      textInput("blanks", "Blank wells (comma separated)", value = "A01,A02,A03"),
      actionButton("run", "Run Analysis")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Explanation",
                 withMathJax(),
                 tags$h3("PE Diagnostics"),
                 tags$p("This Script analyzes absorption and fluorescence data from all observations in my study of ",
                        tags$a(href="https://algaex.pe/en/raw-materials/alga-chondracanthus-chamissoi/", "C. Chamissoi")),
                 tags$img(src="https://photos.app.goo.gl/MBReKTdw4xzJo2SB7", height = "200px"),
                 tags$br(),
                 tags$p("Scripts are sourced from the ", tags$code("process_plate_run"), " R script and called upon in this app."),
                 tags$p("You must store the script in the working directory as well as all input and output files."),
                 tags$p("Go ahead and pick the directory you want to use in the sidebar on the left."),
                 
                 tags$h4("PART i. Calculating PE"),
                 tags$p("The ", tags$code("process_plate_run"), " function processes spectral absorbance data from a multi-wavelength microplate reader assay to quant
