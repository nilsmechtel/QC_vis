library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(shinycssloaders)
library(plotly)


## header ##
header <- dashboardHeader(title = 'Quality Control',
                          titleWidth = 250)


## sidebar ##
sidebar <- dashboardSidebar(disable = TRUE)

## body ##
author.help <- fluidRow(
  column(9,
         uiOutput(outputId = 'help.text')
  ),
  column(3,
         tags$h5(tags$strong('Designed by Nils Mechtel'), align = 'right'),
         tags$h5(tags$strong('Latest update: 14.03.2021'), align = 'right'),
         tags$img(src='MCTP_Logo.png', width='60%', align = 'right')
  )
)

upload.box <- box(width = 12,
                  fluidRow(
                    # new row
                    column(5,
                           fluidRow(
                             column(9,
                                    fileInput(inputId = 'files', label = 'Upload QC results',
                                              multiple = TRUE, accept = '.txt')
                             ),
                             column(3, style = 'margin-top: 25px;',
                                    actionButton(inputId = 'add', label = 'Add files', width = '120px')
                             )
                           ),
                           fluidRow(
                             column(9,
                                    uiOutput(outputId = 'uploadedFiles')
                             ),
                             column(3,
                                    actionButton(inputId = 'calc', label = 'Start processing', width = '120px')
                             )
                           ),
                           fluidRow(
                             tags$br(),
                             tags$br(),
                             column(9,
                                    checkboxInput(inputId = 'help',
                                                  label = HTML('<b>Show help text</b>'))
                             ),
                             column(3,
                                    actionButton(inputId = 'save', label = 'Save data', width = '120px')
                             )
                           )
                      
                    ),
                    column(3, align = 'right',
                           radioButtons(inputId = 'sort', label = 'Order compounds by',
                                        choices = list('Name', 'Retention time'),
                                        selected = 'Name')
                    ),
                    column(4,
                           uiOutput(outputId = 'compoundPicker'),
                           fileInput(inputId = 'compounds.file', label = 'Upload compounds',
                                     multiple = FALSE, accept = '.xlsx')
                    ),
                    # new row
                    column(12,
                           DT::dataTableOutput(outputId = 'overview') %>% withSpinner(color='#428bca')
                    )
                  )
)

vis.tabBox <- tabBox(width = 12, id = 'vis.tab',
                     tabPanel('PCA', uiOutput(outputId = 'PCA')),
                     tabPanel('Alert', uiOutput(outputId = 'Alert')),
                     tabPanel('Match Frequency', uiOutput(outputId = 'Match.Frequency')),
                     tabPanel('Matches', uiOutput(outputId = 'Matches')),
                     tabPanel('Regression', uiOutput(outputId = 'Regression')),
                     tabPanel('Control Cards', uiOutput(outputId = 'Control.Cards')),
                     tabPanel('Ratios', uiOutput(outputId = 'Ratios')),
                     tabPanel('Delete', uiOutput(outputId = 'Delete'))
)

body <- dashboardBody(
  includeCSS('style.css'),
  author.help,
  
  # Data upload #
  fluidRow(
    upload.box
  ),
  
  # compound picker for Ratios
  fluidRow(
    uiOutput(outputId = 'compoundPickerRatios')
  ),
  
  # Visualisation #
  tags$style(HTML('.man_made_class{text-align: center;}')),
  fluidRow(
    vis.tabBox
  )
)

## Dashboard Page ##
dashboardPage(
  skin = 'black',
  header,
  sidebar,
  body
)