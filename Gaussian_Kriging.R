rm(list = ls())
setwd("/Users/nagp/Desktop/env stats/project/Final_scripts/")
source("LMC_functions.R")

## Empirical distribution calculation for loc 1 in test set
num_sim = 50
index_loc = 57
result = empirical.distbn_2d(num_sim,index_loc)
plot(result$x_axis,result$proba)

pred_list_Gaussian = rep(NA,num_sim)
for(sim in 1:num_sim ){
  file_path = paste0("synthetic_data_simulations/training_data/2D_nonGaussian_1200_projection_",sim,"train.csv")
  data_train = read.csv(file_path, header = T)
  file_path = paste0("synthetic_data_simulations/testing_data/2D_nonGaussian_1200_projection_",sim,"test.csv")
  data_test = read.csv(file_path, header = T)
  init.ind<-c(1,1,1,1,1,1,0,0)
  indmat.estim<-optim_indmat_loglik(init.ind)
  init.lmc<-c(indmat.estim$par[c(1:6)],0,0,0,0,indmat.estim$par[7:8])
  fit.Model.lmc <- optim_lmc.loglikelihood(par=init.lmc)
  pred = lmc.pred.summary(fit.Model.lmc$par, data_train, data_test)
  pred_list_Gaussian[sim] = list(pred)
}

saveRDS(pred_list_Gaussian, "pred_list_Gaussian.rds")


