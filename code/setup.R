source("sim_dat.R")
source("rank_FDR.R")
source("i_FDR_pair.R")
source("single-experiment.R")


suppressPackageStartupMessages({
  library(magrittr)
  library(splines)
  library(robustbase)
  library(ggplot2)
  library(caTools)
  library(randomForest)
  library(doParallel)
  library(foreach)
  library(tidyr)
  library(dplyr)
  library(rpart)
  library(quantregForest)
  library(partykit)
  #library(causalTree)
})

cl <- makeCluster(detectCores())
registerDoParallel(cl)

