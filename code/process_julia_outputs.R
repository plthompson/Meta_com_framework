library(tidyverse)
library(data.table)
library(cowplot)
library(scales)

all.results.df <- data.frame()
for(i in 1:15){
  print(i)
  load(paste("./outputs/outputfile_", i, ".RData", sep = ""))
  
  model.df <- model.df %>% 
    filter(Time > 20)
  
  
  model_summary<- model.df %>% 
    mutate(env_match_raw = (1-abs(z-env))*N) %>% 
    group_by(dispersal, sig_niche, alpha, Time, Patch) %>% 
    dplyr::summarise(alpha_div = sum(N>0),
                     alpha_N = sum(N),
                     env_match_raw = sum(env_match_raw)/sum(N)) %>% 
    ungroup() %>% 
    group_by(dispersal, sig_niche, alpha, Patch) %>% 
    dplyr::summarise(alpha_div = mean(alpha_div),
                     stability = mean(alpha_N, na.rm = TRUE)/sd(alpha_N,na.rm = TRUE),
                     alpha_N = mean(alpha_N),
                     env_match_raw = mean(env_match_raw)) %>% 
    ungroup() %>% 
    mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
    group_by(dispersal, sig_niche, alpha) %>% 
    dplyr::summarise(alpha_div = mean(alpha_div),
                     stability = mean(stability, na.rm = TRUE),
                     alpha_N = mean(alpha_N),
                     env_match_raw = mean(env_match_raw))
  
  model_summary_gamma<- model.df %>% 
    group_by(dispersal, sig_niche, alpha, Time, Species) %>% 
    dplyr::summarise(N = sum(N)) %>% 
    ungroup() %>% 
    group_by(dispersal, sig_niche, alpha, Time) %>% 
    dplyr::summarise(gamma_div = sum(N>0),
                     N = sum(N)) %>% 
    ungroup() %>% 
    group_by(dispersal, sig_niche, alpha) %>% 
    dplyr::summarise(stability_gamma = mean(N)/sd(N),
                     gamma_div = mean(gamma_div))
  
  model_summary_t_gamma<- model.df %>% 
    group_by(dispersal, sig_niche, alpha, Patch, Species) %>% 
    dplyr::summarise(N = sum(N)) %>% 
    ungroup() %>% 
    group_by(dispersal, sig_niche, alpha, Patch) %>% 
    dplyr::summarise(gamma_div_time = sum(N>0)) %>% 
    ungroup() %>% 
    group_by(dispersal, sig_niche, alpha) %>% 
    dplyr::summarise(gamma_div_time = mean(gamma_div_time))
  
  coexist.df <- model.df %>% 
    group_by(dispersal, sig_niche, alpha, Time) %>% 
    dplyr::summarise(dens_fit_cov = cov(N,lambda),
                     env_cov = cov(N,env_match),
                     abund_densfit_cov = cov(N,1/density)) %>% 
    ungroup() %>% 
    group_by(dispersal, sig_niche, alpha) %>% 
    dplyr::summarise(dens_fit_cov = mean(dens_fit_cov,na.rm = TRUE),
                     abund_densfit_cov = mean(abund_densfit_cov,na.rm = TRUE),
                     env_cov = mean(env_cov,na.rm = TRUE))
  
  all_model_summary <- left_join(left_join(left_join(model_summary,model_summary_gamma),
                                           model_summary_t_gamma),
                                 coexist.df)
  
  all_model_summary <- all_model_summary %>% 
    mutate(beta_time = gamma_div_time/alpha_div, 
           beta_space = gamma_div/alpha_div)
  
  all_model_summary$rep <- i
  
  all.results.df <- bind_rows(all.results.df, all_model_summary)
}

#save(all.results.df, file = "./data/julia_outputs.RData")
#load data ####