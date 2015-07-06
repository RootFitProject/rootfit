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
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(

      helpText("ROOT-FIT is descriptive model of fitted quadratic growth functions to root system architecture dynamics. ROOT-FIT was developped by Magdalena Julkowska and Christa Testerink. The R version was made by Guillaume Lobet."),
      tags$hr(),      
      
      selectInput("software", label = "Data coming from:",
                  choices = c("RSML_reader", "EZ-Rhizo"), selected = "EZ-Rhizo"),
      
      fileInput('data_file', 'Choose CSV File', accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),

      selectInput("fitting", label = "Fitting function:",
                  choices = c("Find best", "Linear", "Quadratic", "Exponential"), selected = "Find best"),

      actionButton(inputId = "runROOTFIT", label="Unleash ROOT-FIT"),
      
      tags$hr(),      
      h4("More options"),
      
      selectInput("plotting", label = "Plot by:",
                  choices = c("Treatment", "Genotypes"), selected = "Treatment"),
      selectInput("n_plot", label = "Number of genotype to plot",
                  choices = c(1:6), selected = 3),
      selectInput("sample", label = "Find best fit on # samples",
                  choices = seq(10, 100, by=10), selected = 10),
      
      tags$hr(),      
      textInput("name_sep", label = "Name separator", value = "_"),
      textInput("gen_tag", label = "Genotype tag", value = "GEN"),
      textInput("tr_tag", label = "Treatment tag", value = "TR"),
      textInput("dag_tag", label = "Day tag", value = "DAS"),
      
      tags$hr(),      
      h5("Reference:"),
      helpText("Capturing Arabidopsis root architecture dynamics with ROOT-FIT reveals diversity in responses to salinity. 
               \n Julkowska M. M., Hoefsloot H. C. J., Mol S., Feron R., de Boer G. J., Haring M. A. & Testerink C. 
               \n (2014), Plant Physiology. 166, 1387-1402")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(     
        tabPanel("Model visualization",
                 helpText("Comparison between the data and the ROOT-FIT model. Dotted lines represent the modelled values. Plain lines represent the data. Error bars represent the standart devation."),
                 downloadButton('downloadPlot', 'Download Plot'),                 
                 tags$hr(),                       
                 plotOutput("modelPlot"),
                 value=1
        ),
        
        tabPanel("Factor Comparison",
                 helpText("Comparison between the factor estimated for each parameters"),        
                 downloadButton('downloadPlot2', 'Download Plot'),                 
                 plotOutput("factorPlot"),
                 value=2
        ),        

        tabPanel("Fitting results",
                 helpText("Table with the sum of r-squared value obtained with each fitting methods. The method yielding the highest value was used in the fitting"),        
                 tags$hr(),
                 tableOutput('fitting_results'),
                 value = 3),
        
        tabPanel("Download results",
                 helpText("Results from the simulation with the estimated factors"),                 
                 tags$hr(),
                 downloadButton('downloadData', 'Download'),
                 tags$hr(),                 
                 tableOutput('results'),
                 value = 4),
        
        tabPanel("Download factors",
                 helpText("Value of the different estimated factors."),
                 tags$hr(),
                 downloadButton('downloadFactors', 'Download'),
                 tags$hr(),                 
                 tableOutput('factors'),
                 value = 5),        
        id="tabs1"
      )
    )
  )
))