library(tidyverse)
library(data.table)
library(cowplot)
library(scales)

#load data ####
load("./data/julia_outputs.RData")

#set plot theme####
plt_theme <- theme_bw() + theme(
  plot.background = element_blank(),#element_rect(),    # Background of the entire plot
  
  panel.background = element_rect(),   # Background of plotting area
  panel.border = element_rect(),       # Border around plotting area.
  # fill argument should be NA
  
  panel.grid = element_blank(),         # All grid lines
  panel.grid.major = element_blank(),   # Major grid lines
  panel.grid.minor = element_blank(),   # Minor grid lines
  
  panel.grid.major.x = element_blank(), # Vertical major grid lines
  panel.grid.major.y = element_blank(), # Horizontal major grid lines
  panel.grid.minor.x = element_blank(), # Vertical minor grid lines
  panel.grid.minor.y = element_blank(),  # Vertical major grid lines
  
  strip.background = element_rect(colour=NA, fill=NA),
  strip.text.x = element_text(size = 12),
  strip.text.y = element_text(size = 12),
  axis.title=element_text(size = 12)
)

theme_set(plt_theme)
options(scipen = 99999)

#clean up data####

all.results.df <- do.call(data.frame,lapply(all.results.df, function(x) replace(x, is.infinite(x),NA)))
all.results.df <- do.call(data.frame,lapply(all.results.df, function(x) replace(x, is.nan(x),NA)))

all.results.df[all.results.df$sig_niche == unique(all.results.df$sig_niche)[9] & all.results.df$dispersal == unique(all.results.df$dispersal)[1],]

data.median<- all.results.df %>%
  filter(dispersal < max(dispersal)) %>% 
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

#Figure 4####
data.median %>%
  select(alpha_div, gamma_div, beta_space, beta_time, dispersal, sig_niche, competition) %>% 
  gather(key = measure, value = value, alpha_div:beta_time) %>%
  mutate(scale_greek = recode(measure, gamma_div = "\u03B3",alpha_div = "\u03B1", beta_space = "\u03B2 space",  beta_time = "\u03B2 time")) %>% 
  ggplot(aes(x = dispersal, y = sig_niche, fill = value)) +
  geom_tile()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = c(0.01,0.1,1,10))+
  ylab(expression("abiotic niche breadth ("~sigma[italic(i)]~")"))+
  xlab(expression(paste("dispersal (",italic(a[i]),")",sep = "")))+
  scale_fill_viridis_c(option = "B", trans = "log10", name = "species\nrichness", breaks = c(1,2,5,10,20,45))+
  facet_grid(competition~scale_greek)
ggsave("./figures/Figure 4.pdf", height = 12*0.6, width = 14*0.6)

#Figure 3####
hold <- all.results.df %>%
  filter(sig_niche > 9 | sig_niche > 0.4 & sig_niche < 0.5) %>%
  filter(dispersal < max(dispersal)) %>% 
  mutate(competition = recode(alpha, `0.5` = "stabilizing competition",
                              `1.5` = "mixed equilibria",
                              equal = "equal competition",
                              patch_dynamics = "CC trade-off")) %>%
  mutate(competition = factor(competition, levels = c("equal competition", "stabilizing competition", "mixed equilibria", "CC trade-off"), ordered = TRUE)) %>% 
  select(dispersal, alpha_div, gamma_div, beta_space, beta_time, competition, sig_niche, rep) %>% 
  group_by(dispersal, sig_niche, competition) %>% 
  mutate(reps = length(unique(rep))) %>% 
  filter(reps > 10) %>% 
  gather(key = scale, value = diversity, -dispersal, -competition, -rep, -reps, -sig_niche) %>% 
  group_by(dispersal, competition, sig_niche, scale) %>% 
  summarise(diversity_low = quantile(diversity, probs = 0.25), diversity_high = quantile(diversity, probs = 0.75), diversity = quantile(diversity, probs = 0.5)) %>% 
  mutate(scale_greek = recode(scale, gamma_div = "\u03B3",alpha_div = "\u03B1", beta_space = "\u03B2 space",  beta_time = "\u03B2 time"))

hold$sig_niche_names <- "narrow abiotic niche (\u03C3 = 0.5)"
hold$sig_niche_names[hold$sig_niche > 9] <- "wide abiotic niche (\u03C3 = 10)"

ggplot(hold, aes(x = dispersal, y = diversity, color = scale_greek, fill = scale_greek, group = interaction(sig_niche, scale_greek))) +
  geom_ribbon(aes(ymin = diversity_low, ymax = diversity_high), color = NA, alpha = 0.25)+
  geom_line(size = 1)+
  scale_y_log10()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(values = c("black","dodgerblue4", "orange1", "red1"), name = "scale")+
  scale_fill_manual(values = c("black","dodgerblue4", "orange1", "red1"), guide = FALSE)+
  ylab("species richness")+
  xlab(expression(paste("dispersal (",italic(a[i]),")",sep = "")))+
  facet_grid(competition~sig_niche_names)
ggsave("./figures/Figure 3.png", width = 8*0.8, height = 13*0.55)
ggsave("./figures/Figure 3.pdf", width = 8*0.8, height = 13*0.55)

#Figure 5####
data.median %>% 
  mutate(ab_niche = factor(round(sig_niche,2), levels = unique(round(sig_niche,2))[order(unique(round(sig_niche,2)), decreasing = F)], ordered = TRUE)) %>%
  mutate(ab_niche = recode(ab_niche,
                           "0.02" = "0.02",
                           "0.05" = "0.05",
                           "0.1" = "0.1",
                           "0.22" = "0.2",
                           "0.46" = "0.5",
                           "1" = "1",
                           "2.15" = "2",
                           "4.64" = "5",
                           "10" = "10")) %>% 
  ggplot(aes(x = dispersal, y = stability_gamma, fill = ab_niche, group = ab_niche))+
  geom_path(aes(color = ab_niche),size = 1)+
  geom_point(pch = 21, size = 2)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  #scale_fill_gradient(low = "gray90", high = "black",trans = "log10", name = expression("abiotic niche breadth ("~sigma[italic(i)]~")"))+
  scale_fill_brewer(palette = "Greys", breaks = c("0.1", "0.5", "1","5", "10"), name = expression("abiotic niche breadth ("~sigma[italic(i)]~")"))+
  scale_color_brewer(palette = "Greys", breaks = c("0.1", "0.5", "1","5", "10"), name = expression("abiotic niche breadth ("~sigma[italic(i)]~")"))+
  #scale_colour_gradient(low = "gray90", high = "black",trans = "log10", name = expression("abiotic niche breadth ("~sigma[italic(i)]~")"))+
  #scale_fill_viridis_c(end = 0.85,trans = "log10", guide = F)+
  #scale_color_viridis_c(end = 0.85, trans = "log10", name = expression("abiotic niche breadth ("~sigma[italic(i)]~")"))+
  facet_wrap(~competition)+
  ylab("metacommunity stability (\u03B3 invariability)")+
  xlab(expression(paste("dispersal (",italic(a[i]),")",sep = "")))
ggsave("./figures/Figure 5.png", height = 6.5*0.8, width = 8*0.9)
ggsave("./figures/Figure 5.pdf", height = 6.5*0.8, width = 8*0.9)

#Figure 6####
data.median %>% 
  mutate(ab_niche = factor(round(sig_niche,2), levels = unique(round(sig_niche,2))[order(unique(round(sig_niche,2)), decreasing = F)], ordered = TRUE)) %>%
  mutate(ab_niche = recode(ab_niche,
                           "0.02" = "0.02",
                           "0.05" = "0.05",
                           "0.1" = "0.1",
                           "0.22" = "0.2",
                           "0.46" = "0.5",
                           "1" = "1",
                           "2.15" = "2",
                           "4.64" = "5",
                           "10" = "10")) %>% 
ggplot(aes(x = alpha_div, y = alpha_N, fill = factor(round(dispersal,5)), group = ab_niche))+
  geom_path(aes(color = ab_niche), size = 1)+
  geom_point(pch = 21, size = 2)+
  ylab("community abundance (N)")+
  xlab("\u03B1 richness")+
  scale_fill_viridis_d(breaks = c(1e-5,1e-4,1e-3,1e-2,1e-1), name = expression(paste("dispersal (",italic(a[i]),")")))+
  scale_color_brewer(palette = "Greys", breaks = c("0.1", "0.5", "1","5", "10"), name = expression("abiotic niche breadth ("~sigma[italic(i)]~")"))+
  facet_wrap(~competition, scales = "free")
ggsave("./figures/Figure 6.png", height = 6.5*0.8, width = 8*0.9)
ggsave("./figures/Figure 6.pdf", height = 6.5*0.8, width = 8*0.9)


#Figure S1####
env.df <- read.csv("./data/landscape_data/env_2.csv", row.names = 1)

env.df$patch <- env.df %>% 
  group_by(x,y) %>% 
  group_indices()

env.sub <- env.df %>% 
  filter(patch %in% c(1,5,20,40,80,100))

ggplot(filter(env.df, time == 1), aes(x = x, y = y, fill = env1))+
  geom_point(pch = 21, size = 3)+
  scale_fill_viridis_c(name = "environment", breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1))+
  coord_equal()
ggsave("./figures/Figure S1.pdf", height = 5, width = 5)

#Figure S2####
ggplot(env.df, aes(x = time, y = env1, group = patch))+
  geom_line(color = "grey")+
  geom_line(data = env.sub, aes(col = factor(patch)), size = 1)+
  scale_color_brewer(palette = "Set1", name = "patch")+
  ylab("environment")
ggsave("./figures/Figure S2.pdf", height = 4, width = 7)



