rm(list = ls())
setwd("/Users/nagp/Desktop/env stats/project/Final_scripts/")
library(ddalpha)
library(ggplot2)
library(RColorBrewer)
library(quantreg)
library(calculus)

# 2d-nonGaussian data generated from tukey g-h transformation over parsimonious matern
num_sim = 50
min_samples_per_class = 20
for(sim in 1:num_sim){
  file_path <- paste0("synthetic_data_simulations/2d_nongaussian_1200_",sim,".csv")
  data <- read.csv(file_path, sep = ",", header = T)
  
  
  # fit <- rq(data$var2 ~ data$var1, tau = 0.5)
  regress = lm(var2  ~ var1, data)
  
  n_model = 5
  # n_splits_per_model = 5
  
  tau = seq(0.05,0.95,length.out = n_model)
  a = rep(NA,n_model)
  b = rep(NA,n_model)
  for(i in 1:n_model){
    fit <- rq(data$var2 ~ data$var1, tau = tau[i])
    a[i] = unname(fit$coefficients[1])
    b[i] = unname(fit$coefficients[2])
  }
  # a <- seq(-2,1.2,length.out = n_model)
  # b <- rep(unname(regress$coefficients[2]),n_model)
  
  df_threshold_list = list()
  threshold = rep(NA,dim(data)[1])
  proj_var1 = rep(NA,dim(data)[1])
  proj_var2 = rep(NA,dim(data)[1])
  for(i in 1:dim(data)[1]){
    
    vec_distance = abs(data$var2[i] - b*data$var1[i] - a)/sqrt(1+b^2)
    index = order(vec_distance)[1]
    
    proj.var1 = (b[index]*data$var2[i] + data$var1[i] - a[index]*b[index])/(1+b[index]^2)
    proj.var2 = b[index]*proj.var1 + a[index]
    
    # r[i] = sqrt(proj.var1^2 + (proj.var2-a[index])^2)
    # theta[i] = atan((proj.var2-a[index])/proj.var1)
    proj_var1[i] = proj.var1
    proj_var2[i] = proj.var2
    threshold[i] = index
  }
  
  data$proj_var1 = proj_var1
  data$proj_var2 = proj_var2
  data$threshold = threshold
  rm(threshold,proj_var1,proj_var2)
  dist_frm_org = rep(NA,dim(data)[1])
  
  for(i in 1:dim(data)[1]){
    dist_frm_org[i] = (data$proj_var2[i] - a[data$threshold[i]])^2 + (data$proj_var1[i]^2)
  }
  
  data$dist_frm_orig = dist_frm_org
  
  data_list = list()
  tot_class_number = 0
  for(i in 1:n_model){
    df = data[data$threshold == i,]
    # threshold2 = rep(NA,dim(df)[1])
    min_dist = min(df$dist_frm_orig)
    max_dist = max(df$dist_frm_orig)
    df = df[order(df$dist_frm_orig),]
    # row.names(df) <- NULL
    if(dim(df)[1] >= min_samples_per_class){
      if(dim(df)[1]%%min_samples_per_class >= 10){
        n_class = dim(df)[1]%/%min_samples_per_class + 1
      }else{
        n_class = dim(df)[1]%/%min_samples_per_class
      }
      }
    else{n_class = 1}
    # for(j in 1:dim(df)[1]){
    #   for(k in 1:n_splits_per_model){
    #     if(df$dist_frm_orig[j] >= (min_dist + (max_dist-min_dist)*(k-1)/n_splits_per_model) & 
    #        df$dist_frm_orig[j] <= (min_dist + (max_dist-min_dist)*(k)/n_splits_per_model)){
    #       threshold2[j] = k
    #     }
    #     # else{
    #     #   # print(df$dist_frm_orig[j])
    #     # }
    #   }
    # }
    batch = c()
    for(class_number in 1:(n_class-1)){
      tot_class_number = tot_class_number + 1
      batch = c(batch,rep((tot_class_number),min_samples_per_class))
    }
    li = length(batch)
    tot_class_number = tot_class_number + 1
    threshold2 = c(batch,rep((tot_class_number),(dim(df)[1] - li)))
    
    
    
    df$threshold2 = threshold2
    # df$constant = rep(((i-1)*n_class),dim(df)[1])
    data_list[[i]] = df
    print(i)
  }
  
  df_threshold2 = do.call(rbind, data_list)
  df_threshold2 = na.omit(df_threshold2)
  
  df_threshold2$class = df_threshold2$threshold2
  df_threshold2$index <- as.numeric(row.names(df_threshold2))
  df_threshold2 = df_threshold2[order(df_threshold2$index), ]
  total_classes = length(unique(df_threshold2$class))
  df = df_threshold2
  Z1_node = matrix(NA, nrow = 1, ncol = total_classes )
  Z2_node = matrix(NA, nrow = 1, ncol = total_classes )
  
  for(i in 1:total_classes){
    df_sample = df[df$class == i,]
    Z1_node[,i] = mean(df_sample$var1)
    Z2_node[,i] = mean(df_sample$var2)
  }
  
  
  ####### dataset projection explanation through plots ##########
  n_model = 5
  # n_splits_per_model = 10
  
  tau = seq(0.05,0.95,length.out = n_model)
  a = rep(NA,n_model)
  b = rep(NA,n_model)
  for(i in 1:n_model){
    fit <- rq(df$var2 ~ df$var1, tau = tau[i])
    a[i] = unname(fit$coefficients[1])
    b[i] = unname(fit$coefficients[2])
  }
  plot(df$var1,df$var2, pch = 20, col = "green", xlab = "Z1", ylab = "Z2")
  abline(a = a[1], b = b[1])
  abline(a = a[2], b = b[2])
  abline(a = a[3], b = b[3])
  abline(a = a[4], b = b[4])
  abline(a = a[5], b = b[5])
  # abline(a = a[6], b = b[6])
  # abline(a = a[7], b = b[7])
  points(Z1_node,Z2_node, pch = 17, col = "blue", cex = 1.7)
  write.csv(x = df_threshold2, file = "2D_biv_matern_8100_projection.csv",row.names = F)
  write.csv(x = df_threshold2,
            file = paste0("synthetic_data_simulations/2D_nonGaussian_1200_projection_",sim,".csv"),
            row.names = F)
}
