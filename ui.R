library(shiny)

shinyUI(fluidPage(
  # Application title
  titlePanel("angsd-wrapper graph"),
  tabsetPanel(
    tabPanel(
      "Thetas",
      sidebarLayout(
        sidebarPanel(
          fileInput('userThetas',
                    label= 'Choose Thetas File'
          ),
          fileInput('userThetasH',
                    label= 'Choose H Thetas File'
          ),
          fileInput('userThetasL',
                    label= 'Choose L Thetas File'
          ),
          selectInput("thetaChoice",
                      label = "Choose estimator of theta to graph", 
                      choices = c("Watterson's Theta", "Pairwise Theta", "Fu and Li's Theta", "Fay's Theta", "Maximum likelihood (L) Theta"),
                      selected = "Watterson's Theta"
          ),
          selectInput("selectionChoice",
                      label = "Choose a neutrality test statistic to graph", 
                      choices = c("Tajima's D", "Fi and Li's D", "Fu and Li's F", "Fay and Wu's H", "Zeng's E"),
                      selected = "Tajima's D"
          ),
          uiOutput('thetaChroms'),
          uiOutput('thetaChromsH'),
          uiOutput('thetaChromsL'),
          hr(),
          checkboxInput("thetaLowess","Theta Lowess", value=FALSE),
          checkboxInput("selectionLowess","Neutrality Test Statistic Lowess", value=FALSE),
          hr(),
          numericInput("WinCenterLow", "Base Start Position", value=0),
          numericInput("WinCenterHigh", "Base End Position", value=10000),
          checkboxInput("subset","Toggle subset data", value=FALSE),
          hr(),
          checkboxInput("rm.nsites", "Remove data where nSites < x", value=FALSE),
          numericInput("nsites", "x:",value=0),
          hr(),
          fileInput('userAnnotations',
                    label= 'Choose GFF File'
          ),
          checkboxInput("annotations","Toggle GFF annotations", value=FALSE)
        ),
        # Show a plot of the thetas
        mainPanel(
          plotOutput("thetaPlot"), 
          hr(),
          plotOutput("thetaPlotH"),
          hr(),
          plotOutput("thetaPlotL"),
          hr(),
          plotOutput("selectionPlot"),
          hr(),
          plotOutput("selectionPlotH"),
          hr(),
          plotOutput("selectionPlotL")
        )
      )
    )
  )
)
)
