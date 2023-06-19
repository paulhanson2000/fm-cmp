library(data.table)
library(shiny)

selectInput(inputId = "locus",
            label = strong("Locus"),
            choices = locus_configs$locus)

radioButtons(inputId = "method1",
             choices = c("SuSiEx", "PAINTOR"),
             label = strong("Methods to compare"), inline=T)
radioButtons(inputId = "method2",
             choices = c("SuSiEx", "PAINTOR"),
             label = NULL, inline=T)
