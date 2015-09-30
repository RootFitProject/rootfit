# Guillaume Lobet - University of Liege


# ROOT-FIT
# Julkowska et al. (2014). 
# Capturing Arabidopsis root architecture dynamics with ROOT-FIT reveals 
# diversity in responses to salinity. 
# Plant Physiology. doi:10.1104/pp.114.248963


library(shiny)

shinyUI(fluidPage(
  
  # Application title
  titlePanel(h1("- ROOT-FIT -")),
  
  fluidRow(
    column(3, wellPanel(
      helpText("ROOT-FIT is descriptive model of fitted growth functions to root system architecture dynamics. ROOT-FIT was developed by Magdalena Julkowska and Christa Testerink. The R version was made by Guillaume Lobet and Magdalena Julkowska."),
      tags$hr(),      
      
      fileInput('data_file', 'Choose CSV File', accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
      
      selectInput("fitting", label = "Fitting function:",
                  choices = c("Find best", "Linear", "Quadratic", "Exponential"), selected = "Quadratic"),
      
#       selectizeInput('to_analyse', 
#                      label = "Variable to analyse",
#                      choices = NULL,
#                      multiple = T, 
#                      options = list(placeholder = 'select a variable')
#       ),      
      
      tags$hr(),      
      
      numericInput("time_of_treatment", label = "Time of treatment", value = 6), # updated with the datafile
      
      selectInput("ref_gen", label = "Reference genotype", choices = c("Please load datafile")), # updated with the datafile
      
      selectInput("ref_cond", label = "Reference treatment", choices = c("Please load datafile")), # updated with the datafile
      
      numericInput("sample", label = "Find best fit on # samples",
                   value = 10, min = 10, max = 100, step = 10)
      
      
    )),
    
    
    column(3, wellPanel(
      h4("Plotting  options"),
      
      selectInput("plotting", label = "Plot by:",
                  choices = c("Treatment", "Genotypes"), selected = "Treatment"),


      #checkboxGroupInput("to_plot",
      #                   "Genotypes to plot",
      #                   c("Select data file" = 1)),
      
      selectizeInput('to_plot', 
                     label = "Genotypes to plot",
                     choices = c("Select data file" = 1),
                     multiple = T, 
                     options = list(placeholder = 'select a genotype')
      ),
      
      numericInput("plot_width", label = "Plot width", value = 600), 
      numericInput("plot_height", label = "Plot height", value = 800), 
      
      actionButton(inputId = "runROOTFIT", label="Unleash ROOT-FIT"),
      
      
      tags$hr(),      
      h5("Reference:"),
      helpText("Capturing Arabidopsis root architecture dynamics with ROOT-FIT reveals diversity in responses to salinity. 
               \n Julkowska M. M., Hoefsloot H. C. J., Mol S., Feron R., de Boer G. J., Haring M. A. & Testerink C. 
               \n (2014), Plant Physiology. 166, 1387-1402")      


      
      
#       tags$hr(),      
#       textInput("name_sep", label = "Name separator", value = "_"),
#       textInput("gen_tag", label = "Genotype tag", value = "GEN"),
#       textInput("tr_tag", label = "Treatment tag", value = "TR"),
#       textInput("dag_tag", label = "Day tag", value = "DAS"),
      

    )),  
  
  

    # Show a plot of the generated distribution
    column(6,
      tabsetPanel(     
        tabPanel("Model visualization",
                 helpText("Comparison between the data and the ROOT-FIT model. Dotted lines represent the modelled values. Plain lines represent the data. Error bars represent the standart devation."),
                 downloadButton('downloadPlot', 'Download Plot'),                 
                 tags$hr(),                       
                 plotOutput("modelPlot", width = "100%"),
                 value=1
        ),
        
        tabPanel("Factor Comparison",
                 helpText("Comparison between the factor estimated for each parameters"),        
                 downloadButton('downloadPlot2', 'Download Plot'),                 
                 plotOutput("factorPlot", width = "100%"),
                 value=2
        ), 
        
        tabPanel("Relative Factor Comparison",
                 helpText("Comparison between the relative factor estimated for each parameters"),        
                 downloadButton('downloadPlot3', 'Download Plot'),                 
                 plotOutput("relativeFactorPlot", width = "100%"),
                 value=3
        ),         

        tabPanel("Fitting results",
                 helpText("Table with the sum of r-squared value obtained with each fitting methods. The method yielding the highest value was used in the fitting"),        
                 tags$hr(),
                 tableOutput('fitting_results'),
                 value = 4),
        
        tabPanel("Simulation results",
                 helpText("Results from the simulation with the estimated factors"),                 
                 tags$hr(),
                 downloadButton('downloadData', 'Download'),
                 tags$hr(),                 
                 tableOutput('results'),
                 value = 5),
        
        tabPanel("Growth factors",
                 helpText("Value of the different estimated factors."),
                 tags$hr(),
                 downloadButton('downloadFactors', 'Download'),
                 tags$hr(),                 
                 tableOutput('factors'),
                 value = 6),  
        
        tabPanel("Relative growth factors",
                 helpText("Value of the different estimated factors, relative to the refence genotype / treatment."),
                 tags$hr(),
                 h4("Relative to genotype"),
                 downloadButton('downloadFactorsRelGen', 'Download'),
                 tableOutput('relfactorsgen'),
                 tags$hr(),                 
                 h4("Relative to treatment"),
                 downloadButton('downloadFactorsRelTr', 'Download'),
                 tableOutput('relfactorstr'),
                 value = 7),          
        id="tabs1"
      )
    )
  )
))