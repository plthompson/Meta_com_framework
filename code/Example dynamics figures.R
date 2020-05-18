library(tidyverse)

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
  
  strip.background = element_blank(),
  strip.text.x = element_text(size = 14),
  strip.text.y = element_text(size = 14),
  axis.title=element_text(size = 14)
)

theme_set(plt_theme)

load("./data/outputs/outputfile_1.RData")

model.sub <- model.df %>% filter(Patch == 7)

env.df <- model.sub %>% 
  group_by(Time, Patch) %>% 
  summarise(env = mean(env))

model.sub <- left_join(data.frame(Time = 1:max(env.df$Time)), model.sub)

model.sub <- model.sub %>% 
  complete(Species, Patch, Time, dispersal, alpha, sig_niche, fill = list(N = 0)) 

model.sub <- model.sub %>% 
  group_by(Species) %>%
  mutate(z = max(z, na.rm = TRUE))

model.sub <- model.sub %>% 
  mutate(competition = recode(alpha, `0.5` = "stabilizing competition",
                              `1.5` = "mixed competition",
                              equal = "equal competition",
                              patch_dynamics = "CC trade-off"))

model.sub$competition <- factor(model.sub$competition, levels = c("equal competition", "stabilizing competition", "mixed competition", "CC trade-off"), ordered = TRUE)

model.sub %>% 
  filter(dispersal == unique(model.sub$dispersal)[9], sig_niche %in% unique(model.sub$sig_niche)[c(4,6,10)]) %>% 
  mutate(sig_n2 = round(sig_niche,1)) %>% 
  ggplot(aes(x = Time, y = N, color = z, group = Species))+
  geom_line()+
  geom_path(data = env.df, aes(y = -5, color = env, group = NULL), size = 3)+
  scale_color_viridis_c(name = "env. optima", limits = c(0,1), end = 0.95)+
  theme_classic()+
  theme(strip.background = element_blank())+
  facet_grid(competition ~ sig_n2)
ggsave("./figures/Figure S3.pdf", width = 8, height = 7)

model.sub <- model.df %>% filter(Patch == 4)

env.df <- model.sub %>% 
  group_by(Time, Patch) %>% 
  summarise(env = mean(env))

model.sub <- left_join(data.frame(Time = 1:max(env.df$Time)), model.sub)

model.sub <- model.sub %>% 
  complete(Species, Patch, Time, dispersal, alpha, sig_niche, fill = list(N = 0)) 

model.sub <- model.sub %>% 
  group_by(Species) %>%
  mutate(z = max(z, na.rm = TRUE))

model.sub <- model.sub %>% 
  mutate(competition = recode(alpha, `0.5` = "stabilizing competition",
                              `1.5` = "mixed competition",
                              equal = "equal competition",
                              patch_dynamics = "CC trade-off"))

model.sub$competition <- factor(model.sub$competition, levels = c("equal competition", "stabilizing competition", "mixed competition", "CC trade-off"), ordered = TRUE)


dispersal.contrast <- model.sub %>% 
  filter(dispersal %in% unique(model.sub$dispersal)[c(1,9,15)], sig_niche %in% unique(model.sub$sig_niche)[6]) %>% 
  mutate(dispersal2 = formatC(dispersal, format = "e", digits = 0))

dispersal.contrast$dispersal2 <- factor(dispersal.contrast$dispersal2, levels = c("1e-05","5e-03","5e-01"), ordered = TRUE)

ggplot(dispersal.contrast, aes(x = Time, y = N, color = z, group = Species))+
  geom_line()+
  geom_path(data = env.df, aes(y = -5, color = env, group = NULL), size = 3)+
  scale_color_viridis_c(name = "env. optima", limits = c(0,1), end = 0.95)+
  theme_classic()+
  theme(strip.background = element_blank())+
  facet_grid(competition ~ dispersal2)
ggsave("./figures/Fig S4.pdf", width = 8, height = 7)
