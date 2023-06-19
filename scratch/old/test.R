library(shiny)

ui <- fluidPage(
  radioButtons(inputId = "plot_type" , label = "Select the plot", choices = c("scatter", "bar", "hist" )),
  plotOutput("myplot")
  
)

server <- function(input, output, session) {
  
  output$myplot <- renderPlot({
    if (input$plot_type == "scatter") {
      ggplot(mtcars, aes(wt, mpg)) + geom_point()
    } else if (input$plot_type == "bar") {
      ggplot(mtcars, aes(cyl)) + geom_bar()
    }
      else if (input$plot_type == "hist") {
      ggplot(mtcars, aes(mpg)) + geom_histogram()
    }
  })
  
}

shinyApp(ui, server)
