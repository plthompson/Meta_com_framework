landscape_generate<-function(patches = 40, time_steps = 2000, env1Scale = 500, spatial_env = TRUE, temporal_env = TRUE, burn_in = 200) {
  #require(igraph)
  library(ggplot2)
  library(som.nn)
  #require(viridis)
  Tmax<-time_steps +burn_in
  
  #spatially explicit metacomm####
  # if(ex.space == TRUE) {
  #   success<-F
  #   while(!success){
  #     landscape<-round(data.frame(x = runif(patches, min = 1, max = 1000), y = runif(patches, min = 1, max = 1000)))
  #     distance_mat1<-as.matrix(dist(landscape,method = "euclidean",diag = T,upper=T))
  #     
  #     distance_mat<-1*(distance_mat1<300)
  #     diag(distance_mat)<-0
  #     connections<-distance_mat
  #     distance_mat[upper.tri(distance_mat)]<-0
  #     
  #     graph<-as.undirected(graph.adjacency(distance_mat))
  #     graph<-set.vertex.attribute(graph,"x coordinate",value=landscape$x)
  #     graph<-set.vertex.attribute(graph,"y coordinate",value=landscape$y)
  #     graph<-set.edge.attribute(graph,"weight",value=distance_mat1[cbind(as.numeric(get.edgelist(graph)[,1]),  as.numeric(get.edgelist(graph)[,2]))])
  #     
  #     
  #     if(components(graph)$no == 1 & sum(duplicated(landscape))==0){success<-T}}
  #   
  #   d<-shortest.paths(graph)
  #   diag(d)<-0
  #   
  #   d_exp<-exp(-0.005*d) - diag(nrow(d))
  #   d_exp[d_exp==1]<-0
  #   disp_mat <- apply(d_exp, 1, function(x) x/sum(x))
  # } else {
  #   n <- patches # number of points you want on the unit circle
  #   landscape <- 1+(data.frame(t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n)))))+1)/2*999
  #   names(landscape)<-c("x","y")
  #   
  #   distance_mat1<-as.matrix(dist(landscape,method = "euclidean",diag = T,upper=T))
  #   
  #   distance_mat<-1*(distance_mat1==unique(c(distance_mat1))[order(unique(c(distance_mat1)))[2]])
  #   diag(distance_mat)<-0
  #   connections<-distance_mat
  #   
  #   graph<-as.undirected(graph.adjacency(distance_mat))
  #   graph<-set.vertex.attribute(graph,"x coordinate",value=landscape$x)
  #   graph<-set.vertex.attribute(graph,"y coordinate",value=landscape$y)
  #   graph<-set.edge.attribute(graph,"weight",value=distance_mat1[cbind(as.numeric(get.edgelist(graph)[,1]),  as.numeric(get.edgelist(graph)[,2]))])
  #   
  #   disp_mat <- connections*0.5
  # }
  
  repeat {
  landscape<-round(data.frame(x = runif(patches, min = 1, max = 100), y = runif(patches, min = 1, max = 100)))
  if(dim(unique(landscape))[1] == dim(landscape)[1]) {break}
  }
  
  dist_mat <- as.matrix(dist.torus(coors = landscape))
  
  disp_mat<-exp(-0.1*dist_mat)
  diag(disp_mat) <- 0
  disp_mat <- apply(disp_mat, 1, function(x) x/sum(x))
  
  
  #generate spatially autocorrelated environment####
  require(RandomFields)
  require(vegan)
  require(dplyr)
  repeat {
  model <- RMexp(var=0.5, scale=env1Scale) + # with variance 4 and scale 10
    RMnugget(var=0) + # nugget
    RMtrend(mean=0.05) # and mean
  
  RF<-RFsimulate(model = model,x = landscape$x*10, y = landscape$y*10, T = (1+burn_in):Tmax)
  env.df<-data.frame(env1 = decostand(RF$variable1,method = "range"), x= landscape$x, y = landscape$y, time = rep((1+burn_in):Tmax,each = patches))
  burn.env1<-env.df[env.df$time == burn_in+1,]
  if((max(burn.env1$env1)-min(burn.env1$env1)) > 0.6) {break}
  }
  burn.env1$time<-NULL
  burn.env<-left_join(data.frame(time = rep(1:burn_in,each = patches),x = landscape$x, y = landscape$y), burn.env1)
  env.df<-bind_rows(burn.env,env.df)
  
  model <- RMexp(var=0.5, scale=100) + # with variance 4 and scale 10
    RMnugget(var=0) + # nugget
    RMtrend(mean=0.05) # and mean
  RF<-RFsimulate(model = model,x = landscape$x, y = landscape$y, T = 1)
  env.df$env2<-c(decostand(RF$variable1,method = "range"))
  
  ecum <- ecdf(env.df$env1)
  env.cum <- ecum(env.df$env1)
  env.df$env1 <- env.cum
  
  ecum <- ecdf(env.df$env2)
  env.cum <- ecum(env.df$env2)
  env.df$env2 <- env.cum
  
  ggplot(env.df, aes(x = time, y = env1, group = interaction(x,y), color=interaction(x,y)))+
    geom_line()+
    scale_color_viridis_d(guide=F)
  
  #par(mfrow=c(1,2))
  #plot.igraph(graph,layout = data.matrix(landscape),vertex.color=viridis(100)[1+env.df$env1[env.df$time==1]*99], vertex.label=NA, main = "env1")
  #plot.igraph(graph,layout = data.matrix(landscape),vertex.color=viridis(100,option = "B")[1+env.df$env2[env.df$time==1]*99], vertex.label=NA, main = "env2")
  #par(mfrow = c(1,1))
  
  return(list(env.df=env.df,disp_mat =disp_mat,landscape = landscape, Tmax  = Tmax, burn_in = burn_in))
}

trait_generate<-function(species = 25, type, sig_p2 = 0.5){
  #environmental trait
  env.traits<-data.frame(z1 = seq(from = 0,to = 1,length = species), z2 = runif(n = species,min = 0,max = 1))
  env.traits <- env.traits[sample(1:species, size = species, replace = F),]
  
  bdiag1<- 1
  if(type == "no_inter"){
    B <- matrix(0, nrow = species, ncol = species)
    
  }
  
  if(type == "neutral"){
    B <- matrix(1, nrow = species, ncol = species)
  }
  
  if(type == "competitive"){
    B <- matrix(runif(n = species*species)*0.5, nrow = species, ncol = species)
  }
  
  if(type == "priority"){
    B <- matrix(runif(n = species*species)*1.5, nrow = species, ncol = species)
  }
  
  if(type == "patch_dynamics"){
    B <- matrix(runif(n = species*species)*0.5, nrow = species, ncol = species)
    B[1:round(species*0.3),]<-1.25+runif(n=round(species*0.3),min = -0.25, max = 0.25)
    B[lower.tri(B)]<- matrix(runif(n = species*species)*0.5, nrow = species, ncol = species)[lower.tri(matrix(runif(n = species*species)*0.5, nrow = species, ncol = species))]
  }
  diag(B)<-bdiag1
  B <- B*0.05
  
  return(list(B=B, env.traits=env.traits))
}