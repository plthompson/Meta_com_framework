library(data.table)
reps <- 30
for(i in 1:reps){
  print(i)
  #model.df <- read_csv("/Volumes/LaCie Mac/outputs/outputfile_1.csv")

  model.df <- fread(paste("/Volumes/LaCie Mac/outputs/outputfile_", i,".csv", sep = ""), header = T, sep = ',') 
  save(model.df, file = paste("./outputs/outputfile_", i, ".RData", sep =""))
}
                     