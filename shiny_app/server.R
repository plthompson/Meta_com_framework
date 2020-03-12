library(tidyverse)
library(ggExtra)
library(viridis)
library(cowplot)

#setwd(dir = "~/Documents/Scientific Work/Projects/sTURN/meta_com_shiny/")
load("data/Metacom_shiny.RData")

shinyServer(function(input, output, session) {
  data_select<-reactive({
    sub_dynamics<-model.sub %>% filter(as.character(dispersal) == as.character(input$disp_select), as.character(sig_niche) == as.character(input$niche), alpha == input$type)
    
    sub_dynamics <- left_join(data.frame(Time = 1:max(env.df$Time), Patch = rep(1:6,each = max(env.df$Time))),sub_dynamics)
    
    sub_dynamics <- sub_dynamics %>% 
      complete(Species, Patch, Time, fill = list(N = 0)) 
    
    sub_dynamics <- sub_dynamics %>% 
      group_by(Species) %>%
      mutate(z = max(z, na.rm = TRUE))
    
    sub_dynamics
  })
  
  data_select_heat <-reactive({
    hold_results<-results_long
    hold_results <- hold_results[hold_results$measure == input$measure,]
    hold_results$select <- FALSE
    hold_results$select[as.character(hold_results$dispersal) == as.character(input$disp_select) & as.character(hold_results$sig_niche) == as.character(input$niche) & hold_results$alpha == input$type] <- TRUE
    
    hold_results
    })
  
  
  output$plot1 <- renderPlot({
    ggplot(data_select(), aes(x = Time, y = N, color = z, group = Species))+
      geom_line()+
      geom_path(data = env.df, aes(y = -5, color = env, group = NULL), size = 3)+
      facet_wrap(~Patch)+
      scale_color_viridis(option = "D", name = "env. optima", limits = c(0,1), end = 0.95)+
      theme_classic()+
      theme(strip.background = element_blank())
  })
  
  output$plot2 <- renderPlot({
      ggplot(data_select_heat(), aes(x=dispersal, y=as.numeric(sig_niche), fill=value))+
      geom_tile()+
      facet_wrap(~competition)+
      scale_fill_viridis(option = "B", discrete = FALSE, name = "", trans = "log10")+  
      ylab(expression("abiotic niche breadth ("~sigma[italic(i)]~")"))+
      xlab(expression(paste("dispersal (",italic(a[i]),")",sep = "")))+
      scale_y_log10()+
      scale_x_log10()+
      geom_point(aes(color = select), size = 2)+
      scale_color_manual(values = c(NA, "white"), guide = F)+
      theme_classic()+
      theme(strip.background = element_blank())
  })
})
