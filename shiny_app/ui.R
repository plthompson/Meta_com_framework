library(shiny)
library(shinyWidgets)

load("data/Metacom_shiny.RData")

shinyUI(fluidPage(
  h4("Metacommunity model"),
  h5("Patrick L. Thompson and sTURN group"),
  fluidRow(
    column(2,
           wellPanel(
             sliderTextInput(inputId = "disp_select", 
                             label = "Dispersal rate:", 
                             choices = unique(results_long$dispersal), selected = 0.01),
             sliderTextInput(inputId = "niche", 
                             label = "Niche breadth:", 
                             choices = unique(results_long$sig_niche)[order(unique(results_long$sig_niche))], selected = unique(results_long$sig_niche)[order(unique(results_long$sig_niche))][5]),
             selectInput(inputId = "type", 
                         label = "Competition structure:", 
                         list("equal" = "equal", "stabilizing" = "0.5", "mixed" = "1.5",
                              "CC tradeoff" = "patch_dynamics"), selected = "equal"),
             selectInput(inputId = "measure", 
                         label = "Heat plot response variable:", 
                         list("alpha richness" = "alpha_div", 
                              "gamma richness (space)" = "gamma_div",
                              "gamma richness (time)" = "gamma_div_time",
                              "beta richness (space)" = "beta_space",
                              "beta richness (time)" = "beta_time",
                              "community abundance" = "alpha_N",
                              "local abundance invariability" = "stability",
                              "regional abundance invariability" = "stability_gamma",
                              "environmental matching" = "env_match_raw"), selected = "alpha_div")
             )
    ),
  column(6, offset = 0.5,
         "Time series of example dynamics in 6 different habitat patches.",
         "Thin lines represent individual species, coloured by their environmental optima.",
         "The thick straight line at the bottom of each panel shows the local environmental conditions through time.",
         plotOutput("plot1",width = "100%")
  )
  ,
  column(4, offset = 0.5,
         plotOutput("plot2",width = "100%")
  )
)
)
)
