setwd("/Users/nagp/Desktop/env stats/project/Final_scripts/")
library(ddalpha)
library(ggplot2)
library(RColorBrewer)
library(quantreg)
library(calculus)

# 2d-gaussian data generated from parsimonious matern
data <- read.csv("2D_biv_matern_8100.csv", sep = ",", header = T)


# fit <- rq(data$var2 ~ data$var1, tau = 0.5)
regress = lm(var2  ~ var1, data)

n_model = 7
n_splits_per_model = 20

tau = seq(0.1,0.9,length.out = n_model)
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

for(i in 1:n_model){
  df = data[data$threshold == i,]
  threshold2 = rep(NA,dim(df)[1])
  min_dist = min(df$dist_frm_orig)
  max_dist = max(df$dist_frm_orig)
  for(j in 1:dim(df)[1]){
    for(k in 1:n_splits_per_model){
      if(df$dist_frm_orig[j] >= (min_dist + (max_dist-min_dist)*(k-1)/20) & 
         df$dist_frm_orig[j] <= (min_dist + (max_dist-min_dist)*(k)/20)){
        threshold2[j] = k
      }
      # else{
      #   # print(df$dist_frm_orig[j])
      # }
    }
  }
  df$threshold2 = threshold2
  df$constant = rep(((i-1)*n_splits_per_model),dim(df)[1])
  data_list[[i]] = df
  print(i)
}

df_threshold2 = do.call(rbind, data_list)
df_threshold2 = na.omit(df_threshold2)

df_threshold2$class = df_threshold2$constant + df_threshold2$threshold2


write.csv(x = df_threshold2, file = "2D_biv_matern_8100_projection.csv",row.names = F)

data_test <- read.csv("median_results.csv", sep = ",", header = T)
rownames(data_test) <- 1:nrow(data_test)

# data_train <- read.csv("train_data.csv", sep = ",", header = T)
data_train <- df_threshold2
# Z05 = matrix(NA, nrow = 2, ncol = dim(data_test)[1])
# Z25 = matrix(NA, nrow = 2, ncol = dim(data_test)[1])
# Z50 = matrix(NA, nrow = 2, ncol = dim(data_test)[1])
# Z75 = matrix(NA, nrow = 2, ncol = dim(data_test)[1])
# Z95 = matrix(NA, nrow = 2, ncol = dim(data_test)[1])

Z1_node = matrix(NA, nrow = n_model, ncol = n_splits_per_model)
Z2_node = matrix(NA, nrow = n_model, ncol = n_splits_per_model)

for(i in 1:7){
  for(j in 1:20){
    df = data_train[data_train$threshold == i,]
    Z1_node[i,j] = mean(df[df$threshold2 == j,]$proj_var1)
    Z2_node[i,j] = mean(df[df$threshold2 == j,]$proj_var2)
  }
}

marginal.cdf_var1 <- function(z1,index){
  sum = 0
  h = 0.07
  for(i in 1:20){
    #   if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #     z1_res = 0
    #     z2_res = 0
    #   }else{
    a_intercept = data_test$var2_mean[index] - b[data_test$threshold[index]]*data_test$var1_mean[index]
    node1 = (b[data_test$threshold[index]]*Z2_node[data_test$threshold[index],i] + 
               Z1_node[data_test$threshold[index],i] - 
               a_intercept*b[data_test$threshold[index]])/(1+b[data_test$threshold[index]]^2)
    node2 = b[data_test$threshold[index]]*node1 + a_intercept
    # }
    
    p = data_test[index,(14+i)]
    sum = sum + (p*pnorm(z1,mean = node1, sd = h))
    
  }
  return(sum)
}

marginal.cdf_var2 <- function(z1,index){
  sum = 0
  h = 0.07
  for(i in 1:20){
    #   if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #     z1_res = 0
    #     z2_res = 0
    #   }else{
    a_intercept = data_test$var2_mean[index] - b[data_test$threshold[index]]*data_test$var1_mean[index]
    node1 = (b[data_test$threshold[index]]*Z2_node[data_test$threshold[index],i] + 
               Z1_node[data_test$threshold[index],i] - 
               a_intercept*b[data_test$threshold[index]])/(1+b[data_test$threshold[index]]^2)
    node2 = b[data_test$threshold[index]]*node1 + a_intercept
    # }
    
    p = data_test[index,(14+i)]
    sum = sum + (p*pnorm(z1,mean = node2, sd = h))
    
  }
  return(sum)
}


######################################################
##### Conditional cdf and pdfs #######################
######################################################

# cdf of z1,z2 | data
cond.cdf_2d <- function(z,index){
  sum = 0
  h = 0.2
  z1 = z[1]
  z2 = z[2]
  for(i in 1:20){
    #   if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #     z1_res = 0
    #     z2_res = 0
    #   }else{
    
    a_intercept = data_test$var2_mean[index] - b[data_test$threshold[index]]*data_test$var1_mean[index]
    node1 = (b[data_test$threshold[index]]*Z2_node[data_test$threshold[index],i] + 
               Z1_node[data_test$threshold[index],i] - 
               a_intercept*b[data_test$threshold[index]])/(1+b[data_test$threshold[index]]^2)
    node2 = b[data_test$threshold[index]]*node1 + a_intercept
    z1_res = z1 - node1
    z2_res = z2 - node2
    # }
    
    p = data_test[index,(14+i)]
    sum = sum + (p*pbinorm(z1_res,z2_res,var1 = h,var2 = h))
    
  }
  return(sum)
}

# pdf of z1,z2 | data
cond.pdf_2d <- function(z1,z2,index){
  sum = 0
  h = 0.05
  # Z1_node[,62] = estim.model$pred[1]
  # Z2_node[,62] = estim.model$pred[2]
  for(i in 1:20){
    #   if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #     z1_res = 0
    #     z2_res = 0
    #   }else{
    a_intercept = data_test$var2_mean[index] - b[data_test$threshold[index]]*data_test$var1_mean[index]
    node1 = (b[data_test$threshold[index]]*Z2_node[data_test$threshold[index],i] + 
               Z1_node[data_test$threshold[index],i] - 
               a_intercept*b[data_test$threshold[index]])/(1+b[data_test$threshold[index]]^2)
    node2 = b[data_test$threshold[index]]*node1 + a_intercept
    z1_res = z1 - node1
    z2_res = z2 - node2
    # }
    
    p = data_test[index,(14+i)]
    sum = sum + (p*dbinorm(z1_res,z2_res,var1 = h,var2 = h))
    
  }
  return(sum)
}

z2_given = -4.270493
index1 = 1

# pdf of z2 | data
marginal.z2 <- function(z2,index){
  sum = 0
  h = 0.2
  for(i in 1:20){
    
    a_intercept = data_test$var2_mean[index] - b[data_test$threshold[index]]*data_test$var1_mean[index]
    node1 = (b[data_test$threshold[index]]*Z2_node[data_test$threshold[index],i] + 
               Z1_node[data_test$threshold[index],i] - 
               a_intercept*b[data_test$threshold[index]])/(1+b[data_test$threshold[index]]^2)
    node2 = b[data_test$threshold[index]]*node1 + a_intercept
    z2_res = (z2 - node2)/h
    p = data_test[index,(14+i)]
    sum = sum + (p*dnorm(z2_res))/h
  }
  return(sum)
}

# cdf of z1|z2,data
marginal.cond.cdf_z1 <- function(z1,z2,index){
  sum = 0
  h = 0.2
  for(i in 1:20){
    #   if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #     z1_res = 0
    #     z2_res = 0
    #   }else{
    a_intercept = data_test$var2_mean[index] - b[data_test$threshold[index]]*data_test$var1_mean[index]
    node1 = (b[data_test$threshold[index]]*Z2_node[data_test$threshold[index],i] + 
               Z1_node[data_test$threshold[index],i] - 
               a_intercept*b[data_test$threshold[index]])/(1+b[data_test$threshold[index]]^2)
    node2 = b[data_test$threshold[index]]*node1 + a_intercept
    z1_res = (z1 - node1)/h
    z2_res = (z2 - node2)/h
    # z1_res = (z1 - Z1_node[,i])/h
    # z2_res = (z2 - Z2_node[,i])/h
    # }
    
    p = data_test[index,(14+i)]
    sum = sum + (p*pnorm(z1_res)*dnorm(z2_res))/h
    
  }
  
  sum = sum/marginal.z2(z2,index)
  return(sum)
}

# quantile of z1|z2,data
marginal.cond.quantile_z1 <- function(q,z2,index){
  # f <- function(z1){
  #   return(marginal.cond.cdf_z1(z1,z2,index))
  # }
  # f2 <- inverse(f, 0.01, 0.999)
  # return(f2(q))
  z1_val = seq(-10,10,length.out = 100)
  z1_cond.cdf = rep(NA,100)
  for(i in 1:100){
    z1_cond.cdf[i] = marginal.cond.cdf_z1(z1_val[i],z2_given,index1)
  }
  return(approx(z1_cond.cdf,z1_val,q)$y)
}


z1_val = seq(-10,10,length.out = 100)
z1_cond.cdf = rep(NA,100)
for(i in 1:100){
  z1_cond.cdf[i] = marginal.cond.cdf_z1(z1_val[i],z2_given,index1)
}

marginal.cond.quantile_z1(0.5,z2_given,index1)





























