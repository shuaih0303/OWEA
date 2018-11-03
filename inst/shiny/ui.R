library(shiny)
library(OWEA)
# Define UI for application that draws a histogram
shinyUI(
  navbarPage(
    
    id = 'model', title = 'OWEA : Design Generator', selected = 'dropout',
    
    
    # UI for crossover dropout 
    
    tabPanel(
      
      value = 'dropout',title = 'Crossover Dropout',
      
      # App title ----
      h3("Crossover Model with Subject Dropout", align = 'center'),
      withMathJax(),
      h4(
      "$$y_{ij}=\\mu + \\zeta_{i} + \\pi_j + \\tau_{d(i,j)} + \\gamma_{d(i,j-1)} + \\epsilon_{ij}$$"),
      
      helpText("$$i:\\text{subject index }, j:\\text{period index }, \\zeta:\\text{subject effect }, 
                     \\pi:\\text{period effect }, \\tau:\\text{treatment effect }, \\gamma:\\text{carryover effect}$$"),
      tags$hr(),
      
      # Sidebar layout with input and output definitions ----
      sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
          
          # Input: Specify the number of periods ----
          fluidRow(
            column(width = 6, offset = 0, 
              selectInput(inputId = "opt_drop", "opt-criterion", 
                      
                      choices = c('D-optimal' = 0, 'A-optimal' = 1))),
          
          
          # Input: Specify the number of runs ----
            column(width = 6, offset = 0,
              numericInput(inputId = "n_drop", "total runs", 16))),
          
          # Input: Specify the number of periods ----
          fluidRow(
            column(width = 6, offset = 0, 
              numericInput(inputId = "p_drop", "periods", 4)),
          
          # Input: Specify the number of treatments ----
            column(width = 6, offset = 0,
              numericInput(inputId = "t_drop", "treatments", 4))),
          
          # Input: Specify the dropout mechanism ----
          textInput(inputId = "drop_drop", "dropout mechanism (separated by ',')",
                    '0,0,0.5,0.5'),
          
          # Include clarifying text ----
          helpText("Note: length of dropout mechanism and
                   number of periods must match!"),
          
          tags$hr(),
          
          # Input: Specify the random seed and max iters ----
          fluidRow(
            column(width = 6, offset = 0, 
              numericInput(inputId = "max_iter_drop", "max iterations", 40)),
            column(width = 6, offset = 0, 
              numericInput(inputId = "seed_drop", "random seed", 232))),
          
          helpText('Note: Usually leave them as default.'),
          
          wellPanel(
          # Input: Specify whether efficiency is needed ----
          checkboxInput(inputId = 'effNeed','Lower Bound Efficiency'),
          
          
          # Include clarifying text ----
          helpText("Note: computation for lower bound may take a while.")
          ),
          
          
          #submitButton(text = 'RUN')
          actionButton(inputId = "run_drop", "RUN", 
                       
                       style = 'background-color: #00BFFF'),
          
          downloadButton("download_drop", "Export")
          
          ),
        
        # Main panel for displaying outputs ----
        mainPanel(
          
          # output exact design ----
          h4('Exact Design for Direct Treatment'),
          
          textOutput('efficiency_drop'),
          
          textOutput('efflb'),
          
          tags$hr(),
          
          tableOutput('exact_design_drop'), 
          
          # output approximate design ----
          tags$hr(),
          
          h4('Approximate Design'),
          textOutput('get_drop'),
          tags$hr(),
          
          tableOutput('approximate_design_drop')
        )
      )
    ),
    
    # ============== UI for proportional model ============== #
    
    tabPanel(
      
      value = 'proportional',title = 'Crossover Proportional',
      
      # App title ----
      h3("Crossover Model Proportional Carryover Effects", align = 'center'),
      
      withMathJax(),
      h4(
        "$$y_{ij}=\\mu + \\zeta_{i} + \\pi_j + \\tau_{d(i,j)} + \\lambda\\tau_{d(i,j-1)} + \\epsilon_{ij}$$"),
      
      helpText("$$i:\\text{subject index }, j:\\text{period index }, \\zeta:\\text{subject effect }, 
               \\pi:\\text{period effect }, \\tau:\\text{treatment effect }, \\lambda:\\text{proportional coefficient}$$"),
      

      tags$hr(),
      
      # Sidebar layout with input and output definitions ----
      sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
          fluidRow(
            column(width = 6, offset = 0,
          # Input: Specify the number of periods ----
          selectInput(inputId = "opt_prop", "opt-criterion", 
                      
                      choices = c('D-optimal' = 0, 'A-optimal' = 1))),
          
          
          # Input: Specify the number of runs ----
            column(width = 6, offset = 0,
          numericInput(inputId = "n_prop", "total runs", 100))
          ),
          
          fluidRow(
            column(width = 6, offset = 0,
              # Input: Specify the number of periods ----
              numericInput(inputId = "p_prop", "periods", 3)),
            column(width = 6, offset = 0,
              # Input: Specify the number of treatments ----
              numericInput(inputId = "t_prop", "treatments", 3))
          ),
          # input :Specify Assumed covariance matrix
          
          tags$i('Specify Initial Values for Parameters:'),
          
                     # Input: Specify initial lambda
                     numericInput(inputId = 'lambda_prop', 'proportion coefficient', value = 0.2),
                     
                     # Input: Specify the initial parameter ----
                     textInput(inputId = "tau_prop", "treatment effects(separated by ',')",
                               '2,2,2'),
          
          tags$hr(),
          
          fluidRow(
            column(width = 6, offset = 0, 
                   numericInput(inputId = "max_iter_prop", "max iterations", 40)),
            column(width = 6, offset = 0, 
                   numericInput(inputId = "seed_prop", "random seed", 123))),
          #helpText('Note: Usually leave it as default.'),

          
          
          #submitButton(text = 'RUN')
          actionButton(inputId = "run_prop", "RUN", 
                       
                       style = 'background-color: #00BFFF'),
          
          downloadButton("download_prop", "Export")
          
          ),
        
        # Main panel for displaying outputs ----
        mainPanel(
          
          # output exact design ----
          h4('Exact Design for Direct Treatment'),
          
          textOutput('efficiency_prop'),
          
          tags$hr(),
          
          tableOutput('exact_design_prop'), 
          
          # output approximate design ----
          tags$hr(),
          
          h4('Approximate Design'),
          textOutput('get_prop'),
          tags$hr(),
          
          tableOutput('approximate_design_prop')
        )
      )
    ),
    
  # UI for interference #
  tabPanel(value = 'interference',title = 'Interference Model',
           
           # App title ----
           h3("Interfence Model",align = 'center'),
           
           withMathJax(),
           h4(
             "$$y_{ij} = \\mu + \\gamma_i + \\tau_{d(i,j)} + \\lambda_{d(i,j-1)} + \\rho_{d(i,j+1)} + \\epsilon_{ij}$$"),
           
           helpText("$$i:\\text{block index }, j:\\text{plot index }, \\gamma:\\text{block effect }, 
                     \\tau:\\text{treatment effect }, \\lambda:\\text{left effect }, \\rho: \\text{right effect}$$"),
           
           
           tags$hr(),
           
           # Sidebar layout with input and output definitions ----
           sidebarLayout(
             
             # Sidebar panel for inputs ----
             sidebarPanel(
               fluidRow(
                 column(width = 6, offset = 0, 
                   # Input: Specify the number of periods ----
                   selectInput(inputId = "opt_inter", "opt-criterion",
                           choices = c('D-optimal' = 0,
                                       'A-optimal' = 1)
               )),
               
                 column(width = 6, offset = 0, 
                   # Input: Specify the number of runs ----
                   numericInput(inputId = "n_inter", "total runs", 100))
               ),
               
              fluidRow(
                column(width = 6, offset = 0, 
               # Input: Specify the number of periods ----
               numericInput(inputId = "p_inter", "block size", 4)
               ),
                column(width = 6, offset = 0,
               # Input: Specify the number of treatments ----
               numericInput(inputId = "t_inter", "treatments", 4))
              ),
              
              tags$hr(),
              
               fluidRow(
                 column(width = 6, offset = 0, 
                        numericInput(inputId = "max_iter_inter", "max iterations", 40)),
                 column(width = 6, offset = 0, 
                        numericInput(inputId = "seed_inter", "random seed", 456))),
               
               
               
               #submitButton(text = 'RUN')
               actionButton(inputId = "run_inter", "RUN",
                            style = 'background-color: #00BFFF'),
               
               downloadButton("download_inter", "Export")
               
               
             ),
             
             # Main panel for displaying outputs ----
             mainPanel(
               
               # output exact design ----
               h4('Exact Design for Direct Treatment'),
               
               textOutput('efficiency_inter'),
               
               tags$hr(),
               
               tableOutput('exact_design_inter'), 
               
               # output approximate design ----
               tags$hr(),
               
               h4('Approximate Design'),
               textOutput('get_inter'),
               tags$hr(),
               
               tableOutput('approximate_design_inter')
             )
           )
  ),
  
  
  tabPanel('About',
           verticalLayout(
           h5("Bug Report: shuaih0303@gmail.com", align = 'center'),
           h5('Author: Shuai Hao, Wei Zheng, Min Yang', align = 'center'),
           tags$hr()
           )
    )
  )
)
  