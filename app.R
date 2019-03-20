library(shiny)

server <- function(input, output, session) {}
ui <- fluidPage(
    # tags$head(includeScript("google-analytics.js")),
    h4("The app has moved. This URL will stop working soon..."),
    h4("Bookmark/save this new site:", a("http://escal.site", href="http://escal.site"))
    )
shinyApp(ui = ui, server = server)