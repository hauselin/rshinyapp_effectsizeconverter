library(shiny)

es <- function(d = NULL, r = NULL, R2 = NULL, f = NULL, oddsratio = NULL, logoddsratio = NULL, auc = NULL, fishersz = NULL, decimal = 2, msg = TRUE) {

  effectsizes <- data.frame(matrix(NA, nrow = length(c(d, r, R2, f, oddsratio, logoddsratio, auc, fishersz)), ncol = 8)) # dataframe version
  names(effectsizes) <- c("d", "r", "R2", "f", "oddsratio", "logoddsratio", "auc", "fishersz")

  if (length(c(d, r, R2, f, oddsratio, logoddsratio, auc, fishersz)) < 1) {
    stop("Please specify one effect size!")
  }

  if (is.numeric(d)) {
    if (msg) {message(paste0("d: ", d, " ")) }
    effectsizes$d <- d
    effectsizes$f <- effectsizes$d / 2
    effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(r)) {
    if (msg) {message(paste0("r: ", r, " ")) }
    effectsizes$d <- (2 * r) / (sqrt(1 - r^2))
    effectsizes$r <- r
    effectsizes$f <- effectsizes$d / 2
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(f)) {
    if (msg) {message(paste0("f: ", f, " ")) }
    effectsizes$d <- f * 2
    effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
    effectsizes$f <- f
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(R2)) {
    if (msg) {message(paste0("R2: ", R2, " ")) }
    effectsizes$r <- sqrt(R2)
    effectsizes$d <- (2 * effectsizes$r) / (sqrt(1 - effectsizes$r^2))
    effectsizes$f <- effectsizes$d / 2
    effectsizes$R2 <- R2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(oddsratio)) {
    if (msg) {message(paste0("odds ratio: ", oddsratio, " "))}
    effectsizes$d <- log(oddsratio) * (sqrt(3) / pi)
    effectsizes$f <- effectsizes$d / 2
    effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- oddsratio
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(logoddsratio)) {
    if (msg) {message(paste0("log odds ratio: ", logoddsratio, " ")) }
    effectsizes$logoddsratio <- logoddsratio
    effectsizes$d <- effectsizes$logoddsratio * (sqrt(3) / pi)
    effectsizes$f <- effectsizes$d / 2
    effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(auc)) { # also known as common language (CL) effect size statistic
    if (msg) {message(paste0("auc: ", auc, " ")) }
    effectsizes$auc <- auc
    effectsizes$d <- stats::qnorm(auc, 0, 1) * sqrt(2) # assumes equal sample size (Ruscio 2008)
    effectsizes$f <- effectsizes$d / 2
    effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
    effectsizes$R2 <- effectsizes$r^ 2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$fishersz <- .5 * log((1+effectsizes$r)/ (1-effectsizes$r))
  } else if (is.numeric(fishersz)) {
    if (msg) {message(paste0("fishersz: ", auc, " ")) }
    effectsizes$fishersz <- fishersz
    effectsizes$r <- (exp(effectsizes$fishersz / 0.5) - 1) / (exp(effectsizes$fishersz / 0.5) + 1)
    effectsizes$d <- (2 * effectsizes$r) / (sqrt(1 - effectsizes$r^2))
    effectsizes$f <- effectsizes$d / 2
    effectsizes$R2 <- R2
    effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
    effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    effectsizes$auc <- stats::pnorm(effectsizes$d/sqrt(2), 0, 1)
  }

  round(effectsizes, decimal)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    tags$head(includeScript("google-analytics.js")),

   # Application title
   titlePanel("Effect size convertor (to be deprecated soon!)"),

   h3(a("NEW: Bookmark and use the faster app here at http://escal.site", href="http://escal.site")),

   h5("The es() function in the",
      a("hausekeep R package", href="https://hauselin.github.io/hausekeep/"),
      "also performs these conversions."),
   h6("Source code is available at ",
      a("GitHub", href="https://github.com/hauselin/rshinyapp_effectsizeconverter"),
      " and R tutorials are available ",
      a("here", href="https://hausetutorials.netlify.com/"),
      "."),

   # Sidebar with a slider input for number of bins
   sidebarLayout(

      sidebarPanel(
        numericInput(inputId = "effectsize", label = "Input effect size", value = 0.3, min = 0, max = 1e6, step = 0.001),
        # selectInput(inputId = "effect_size_type", label = "Input effect size measure",
        #             choices = c("Cohen's d", "correlation r", "r-squared", "Cohen's f",
        #                         "odds ratio", "log odds ratio", "area-under-curve"),
        #             selected = "Cohen's d"),
        radioButtons(inputId = "effect_size_type", label = "Input effect size measure",
                    choices = c("Cohen's d", "correlation r", "r-squared", "Cohen's f",
                                "odds ratio", "log odds ratio", "area-under-curve (auc)",
                                "Fisher's z (z')"),
                    selected = "Cohen's d", inline = FALSE)
      ),

      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Converted effect sizes", tableOutput(outputId = "outputeffectsizes"),
                             tags$p("All conversions assume equal sample size groups."),
                             uiOutput(outputId = "formulae")),
                    tabPanel("References", tags$p(""), htmlOutput("references")),
                    tabPanel("About", tags$p(""), tags$p("Developed and maintained by Hause Lin"),
                             tags$p("hauselin@gmail.com"),
                             tags$p("Lin, H. (2019). hauselin/rshinyapp_effectsizeconverter: shiny effect size converter v0.0.1 (Version v0.0.1). Zenodo. https://doi.org/10.5281/zenodo.2563830"))
        )
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$outputeffectsizes <- renderTable({
     if (input$effect_size_type == "Cohen's d") {
       es(d = input$effectsize, decimal = 4)
     } else if (input$effect_size_type == "correlation r") {
       es(r = input$effectsize, decimal = 4)
     } else if (input$effect_size_type == "r-squared") {
       es(R2 = input$effectsize, decimal = 4)
     } else if (input$effect_size_type == "Cohen's f") {
       es(f = input$effectsize, decimal = 4)
     } else if (input$effect_size_type == "odds ratio") {
       es(oddsratio = input$effectsize, decimal = 4)
     } else if (input$effect_size_type == "log odds ratio") {
       es(logoddsratio = input$effectsize, decimal = 4)
     } else if (input$effect_size_type == "area-under-curve (auc)") {
       es(auc = input$effectsize, decimal = 4)
     } else if (input$effect_size_type == "Fisher's z (z')") {
       es(fishersz = input$effectsize, decimal = 4)
     }
   }, digits = 3)

  output$formulae <- renderUI({
    withMathJax(
      helpText("By convention, Cohen's d of 0.2, 0.5, 0.8 are considered small, medium and large effect sizes respectively (correlation r = 0.1, 0.3, 0.5)."),
      helpText("Cohen's d to correlation r"),
      helpText("$$r = \\frac{d}{\\sqrt{d^2 + 4}}$$"),
      helpText("Cohen's d to Cohen's f"),
      helpText("$$f = \\frac{d}{2}$$"),
      helpText("Cohen's d to log odds ratio"),
      helpText("$$\\lg{odds ratio} = \\frac{d}{\\frac{\\sqrt{3}}{\\pi}}$$"),
      helpText("Cohen's d to area-under-curve (auc)"),
      helpText("$$auc = \\phi\\frac{d}{\\sqrt{2}}, [R syntax: pnorm(d/sqrt(2), 0, 1)]$$"),
      helpText("correlation r to Fisher's z (z')"),
      helpText("$$z' = 0.5 * (\\lg{(1 + r)} - \\lg{(1 - r)})$$")
    )
    })

  output$references <- renderUI({
    HTML(paste(
      "Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). Introduction to meta-analysis. Chichester, West Sussex, UK: Wiley.",
      "",
      "Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.), Hillsdale, NJ: Erlbaum.",
      "",
      "Lin, H. (2019). hauselin/hausekeep: third release (Version v0.0.0.9002-alpha). Zenodo. https://doi.org/10.5281/zenodo.2555874",
      "",
      "Rosenthal, R. (1994). Parametric measures of effect size. In H. Cooper & L. V. Hedges (Eds.), The Handbook of Research Synthesis. New York, NY: Sage.",
      "",
      "Ruscio, J. (2008). A probability-based measure of effect size: Robustness to base rates and other factors. Psychological Methods, 13(1), 19-30. doi:10.1037/1082-989x.13.1.19",
      sep = "<br/>"))
  })

}

# Run the application
shinyApp(ui = ui, server = server)