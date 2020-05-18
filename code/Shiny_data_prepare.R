load("./data/julia_outputs.RData")
library(tidyverse)

all.results.df <- do.call(data.frame,lapply(all.results.df, function(x) replace(x, is.infinite(x),NA)))
all.results.df <- do.call(data.frame,lapply(all.results.df, function(x) replace(x, is.nan(x),NA)))

data.median<- all.results.df %>%
  #filter(dispersal < max(dispersal)) %>% 
  group_by(dispersal, sig_niche, alpha) %>% 
  mutate(reps = length(unique(rep))) %>% 
  select(-rep) %>% 
  summarise_all(funs(median(., na.rm = TRUE))) %>% 
  ungroup() %>% 
  filter(reps > 10)

data.median <- data.median %>% 
  mutate(competition = recode(alpha, `0.5` = "stabilizing competition",
                              `1.5` = "mixed competition",
                              equal = "equal competition",
                              patch_dynamics = "CC trade-off"))

data.median$competition <- factor(data.median$competition, levels = c("equal competition", "stabilizing competition", "mixed competition", "CC trade-off"), ordered = TRUE)


results_long <- data.median %>%
  select(alpha_div, gamma_div, gamma_div_time, beta_space,env_match_raw, stability_gamma, stability, alpha_N, beta_time, dispersal, sig_niche, competition, alpha) %>% 
  gather(key = measure, value = value, alpha_div:beta_time)

load("./data/outputs/outputfile_1.RData")
env.df <- model.df %>% 
  filter(Patch < 7) %>% 
  ungroup() %>% 
  select(Time, Patch, env) %>% 
  unique()

head(model.sub)

model.sub <- model.df %>% 
  filter(Patch < 7) %>% 
  select(N, env, Species, Time, Patch, dispersal, sig_niche, alpha, z) %>% 
  mutate(competition = recode(alpha, `0.5` = "stabilizing competition",
                              `1.5` = "mixed competition",
                              equal = "equal competition",
                              patch_dynamics = "CC trade-off"))

save(model.sub, env.df, results_long, file = "./data/Metacom_shiny_updated.RData")

