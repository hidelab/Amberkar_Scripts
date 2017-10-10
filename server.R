library(shiny)
library(datasets)
mpgData <- mtcars
mpgData$am <- factor(mpgData$am, labels = c("Automatic", "Manual"))

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  formulaText <- reactive({
    paste("mpg ~", input$variable)
  })
  output$caption <- renderText({
    formulaText()
  })
  output$mpgPlot <- renderPlot({
    boxplot(as.formula(formulaText()), 
            data = mpgData,
            outline = input$outliers)
  })
})