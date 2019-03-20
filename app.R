#' This shiny app script redirects visitors to the new site.
#' If you're looking for the code the runs the actualy app, look inside the archive directory

library(shiny)

server <- function(input, output, session) {}
ui <- fluidPage(
    # tags$head(includeScript("google-analytics.js")),
    h4("The app has moved. This URL will stop working soon..."),
    h4("Visit and bookmark/save the new site:", a("http://escal.site", href="http://escal.site"))
    )
shinyApp(ui = ui, server = server)