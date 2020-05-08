source("./code/model_functions.R")

reps <- 30
for(r in 1:reps){
time_steps <- 2000
patches <- 100
species <- 50
env1Scale <- 500

landscape_rep<-landscape_generate(patches = patches, time_steps = time_steps, env1Scale = env1Scale, spatial_env = TRUE, temporal_env = TRUE, burn_in = 200)
env.df<-landscape_rep$env.df
disp_mat<-landscape_rep$disp_mat
landscape<-landscape_rep$landscape

write.csv(x = env.df, file = paste("./data/landscape_data/env_",r,".csv",sep = ""))
write.csv(x = disp_mat, file = paste("./data/landscape_data/disp_mat_",r,".csv",sep = ""))
write.csv(x = landscape, file = paste("./data/landscape_data/landscape_",r,".csv",sep = ""))
write.csv(x = data.frame(Tmax = landscape_rep$Tmax, burn_in = landscape_rep$burn_in), file = paste("./data/landscape_data/time_",r,".csv",sep = ""))
}

