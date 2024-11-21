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
# index_loc = 21
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


## compute picp and mpiw
sim_num = 50
file_path = paste0("synthetic_data_simulations/testing_data/2D_nonGaussian_1200_projection_",sim_num,"test.csv")
data_test = read.csv(file_path, header = T)
pred_list_Gaussian = readRDS("pred_list_Gaussian.rds")
file_path = paste0("Results_DNN/2D_nonGaussian_1200_predictions_",sim_num,".csv")
data_test_dnn = read.csv(file_path, header = T)
file_path = paste0("synthetic_data_simulations/training_data/2D_nonGaussian_1200_projection_",sim_num,"train.csv")
data_train = read.csv(file_path, header = T)
total_classes = dim(data_test_dnn)[2] - 11
Z1_node = matrix(NA, nrow = 1, ncol = total_classes )
Z2_node = matrix(NA, nrow = 1, ncol = total_classes )
for(i in 1:total_classes){
  df_sample = data_train[data_train$class == i,]
  Z1_node[,i] = mean(df_sample$var1)
  Z2_node[,i] = mean(df_sample$var2)
}

median_Gaussian = rep(NA,120)
lb_Gaussian = rep(NA,120)
ub_Gaussian = rep(NA,120)
median_dnn = rep(NA,120)
lb_dnn = rep(NA,120)
ub_dnn = rep(NA,120)
for(index_loc in 1:120){
  preds = pred_list_Gaussian[[sim_num]]$pred[c(index_loc,(120+index_loc)),1]
  cond_var = pred_list_Gaussian[[sim_num]]$conditional_var[c(index_loc,(120+index_loc)),
                                                           c(index_loc,(120+index_loc))]
  ro = cond_var[1,2]/sqrt(cond_var[1,1]*cond_var[2,2])
  cond.Gaussian_1d_pred = preds[1] + ro*(cond_var[1,1]/cond_var[2,2])*
    (data_test$var2[index_loc] - preds[2])
  cond.Gaussian_1d_var = (1-ro^2)*cond_var[1,1]
  median_Gaussian[index_loc] = qnorm(0.5,cond.Gaussian_1d_pred,sqrt(cond.Gaussian_1d_var))
  lb_Gaussian[index_loc] = qnorm(0.025,cond.Gaussian_1d_pred,sqrt(cond.Gaussian_1d_var))
  ub_Gaussian[index_loc] = qnorm(0.975,cond.Gaussian_1d_pred,sqrt(cond.Gaussian_1d_var))
  median_dnn[index_loc] = marginal.cond.quantile_z1(0.5,data_test$var2[index_loc],index_loc,
                                                    data_test_dnn,Z1_node, Z2_node)
  lb_dnn[index_loc] = marginal.cond.quantile_z1(0.025,data_test$var2[index_loc],index_loc,
                                                data_test_dnn,Z1_node, Z2_node)
  ub_dnn[index_loc] = marginal.cond.quantile_z1(0.975,data_test$var2[index_loc],index_loc,
                                                data_test_dnn,Z1_node, Z2_node)
  print(index_loc)
}

## mad 

mean(abs(data_test$var1 - median_Gaussian))
mean(abs(data_test$var1 - median_dnn))

## picp
count_Gaussian = 0
count_dnn = 0
for(i in 1:120){
  if(data_test$var1[i]>lb_Gaussian[i] & data_test$var1[i]<ub_Gaussian[i]) {
    count_Gaussian = count_Gaussian + 1}
  if(data_test$var1[i]>lb_dnn[i] & data_test$var1[i]<ub_dnn[i]) {
    count_dnn = count_dnn + 1}
}
count_Gaussian/120
count_dnn/120



set.seed(12345)
p.sample <- rbeta(10^3,252,450)
x.tilde <- rbinom(10^3,7,p.sample)





