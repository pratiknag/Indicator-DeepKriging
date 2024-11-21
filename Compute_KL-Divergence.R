rm(list = ls())
setwd("/Users/nagp/Desktop/env stats/project/Final_scripts/")
source("LMC_functions.R")
source("density_computation_DNN.R")
library(VGAM)

## function for kl-divergence 
kl_divergence <- function(p, q){
  summation = 0
  for(i in 1:length(p)){
    e = 0.0000001
    summation = summation + (p[i] * log(p[i]+e/q[i]+e))
  }
  return(summation)
}

## calculate distance between functions 
func_distance <- function(p,q){
  return(mean(abs(p-q)))
}


# ## Empirical cdf
# prob_empirical = empirical.distbn_2d(23,27)
# plot(prob_empirical$x_axis,prob_empirical$proba)
# 
# 
#   
# ## DNN 2d-cdf
# prob_dnn = cond.dist.2d_comparison(23,27)
# plot(prob_dnn$x_axis,prob_dnn$proba)
# 
# # kl_divergence(prob_empirical$proba,prob_dnn$proba)
# func_distance(prob_empirical$proba,prob_dnn$proba)
# 
# 
# ## Gaussian 2d-cdf 
# 
# prob_Gaussian = Gaussian.cond.cdf(23,27)
# func_distance(prob_Gaussian$proba,prob_empirical$proba)


## distance between all points for simulation 21 
sim_num = 50 
dist_Gaussian = rep(NA,120)
dist_dnn = rep(NA,120)
for(i in 1:120){
  prob_empirical = empirical.distbn_2d(sim_num,i)
  prob_Gaussian = Gaussian.cond.cdf(sim_num,i)
  prob_dnn = cond.dist.2d_comparison(sim_num,i)
  dist_Gaussian[i] = func_distance(prob_Gaussian$proba,prob_empirical$proba)
  dist_dnn[i] = func_distance(prob_empirical$proba,prob_dnn$proba)
  print(i)
}

boxplot(dist_Gaussian,dist_dnn,
        cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2, xaxt = "n",outline=FALSE)

