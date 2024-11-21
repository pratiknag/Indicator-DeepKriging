### bivariate conditional cdfs
library(VGAM)
cond.cdf_2d <- function(z,index,data, Z1_node, Z2_node){
  sum = 0
  h = 0.2
  z1 = z[1]
  z2 = z[2]
  
  for(i in 12:dim(data)[2]){
    #   if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #     z1_res = 0
    #     z2_res = 0
    #   }else{
    z1_res = z1 - Z1_node[,(i-11)]
    z2_res = z2 - Z2_node[,(i-11)]
    # }
    
    p = data[index,i]
    sum = sum + (p*pbinorm(z1_res,z2_res,var1 = h,var2 = h))
    
  }
  return(sum)
}

cond.pdf_2d <- function(z,index,data, Z1_node, Z2_node){
  sum = 0
  h = 0.05
  z1 = z[1]
  z2 = z[2]
  # Z1_node[,62] = estim.model$pred[1]
  # Z2_node[,62] = estim.model$pred[2]
  for(i in 12:dim(data)[2]){
    #   if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #     z1_res = 0
    #     z2_res = 0
    #   }else{
    z1_res = z1 - Z1_node[,(i-11)]
    z2_res = z2 - Z2_node[,(i-11)]
    # }
    
    p = data[index,i]
    sum = sum + (p*dbinorm(z1_res,z2_res,var1 = h,var2 = h))
    
  }
  return(sum)
}


marginal.z2 <- function(z2,index,data, Z2_node){
  sum = 0
  h = 0.2
  
  for(i in 12:dim(data)[2]){
    
    z2_res = (z2 - Z2_node[,(i-11)])/h
    p = data[index,i]
    sum = sum + (p*dnorm(z2_res))/h
  }
  return(sum)
}

marginal.cond.cdf_z1 <- function(z,index, data, Z1_node, Z2_node){
  sum = 0
  h = 0.2
  z1 = z[1]
  z2 = z[2]
  for(i in 12:dim(data)[2]){
    #   if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #     z1_res = 0
    #     z2_res = 0
    #   }else{
    z1_res = (z1 - Z1_node[,(i-11)])/h
    z2_res = (z2 - Z2_node[,(i-11)])/h
    # }
    
    p = data[index,i]
    sum = sum + (p*pnorm(z1_res)*dnorm(z2_res))/h
    
  }
  
  sum = sum/marginal.z2(z2,index,data,Z2_node)
  return(sum)
}

marginal.cond.quantile_z1 <- function(p,z2_given,index,data,Z1_node,Z2_node){
  # f <- function(z1){
  #   return(marginal.cond.cdf_z1(z1,z2,index))
  # }
  # f2 <- inverse(f, 0.01, 0.999)
  # return(f2(q))
  z1_val = seq(-10,10,length.out = 100)
  z1_cond.cdf = rep(NA,100)
  for(i in 1:100){
    z = c(z1_val[i],z2_given)
    z1_cond.cdf[i] = marginal.cond.cdf_z1(z,index,data, Z1_node, Z2_node)
  }
  return(approx(z1_cond.cdf,z1_val,p)$y)
}

cond.dist.2d_comparison <- function(sim_num, index_loc){
  value_range = seq(-20,40,0.2)
  proba = rep(NA,length(value_range))
  file_path = paste0("Results_DNN/2D_nonGaussian_1200_predictions_",sim_num,".csv")
  data_test = read.csv(file_path, header = T)
  file_path = paste0("synthetic_data_simulations/training_data/2D_nonGaussian_1200_projection_",sim_num,"train.csv")
  data_train = read.csv(file_path, header = T)
  total_classes = dim(data_test)[2] - 11
  Z1_node = matrix(NA, nrow = 1, ncol = total_classes )
  Z2_node = matrix(NA, nrow = 1, ncol = total_classes )
  for(i in 1:total_classes){
    df_sample = data_train[data_train$class == i,]
    Z1_node[,i] = mean(df_sample$var1)
    Z2_node[,i] = mean(df_sample$var2)
  }
  for(value in (1:length(value_range))){
  proba[value] = cond.cdf_2d(rep(value_range[value],2),index_loc,data_test,Z1_node,Z2_node)
  }
  return(list(proba = proba,x_axis = value_range))
}