library(VGAM)
setwd("/Users/nagp/Desktop/env stats/project/Final_scripts/")
df = read.csv("2D_biv_matern_8100_projection.csv", sep = ",", header = T)
df = read.csv("synthetic_data_simulations/2D_nonGaussian_1200_projection_1.csv", 
              sep = ",", header = T)
# thresholds = unique(df$threshold2)
# thresholds = sort(thresholds)
# emp.var1 = var(df$data1)
# emp.var2 = var(df$data2)
# 
# emp.cov = cov(df$data1,df$data2)


# thres.prob = read.csv("prediction_2D.csv", sep = ",", header = T)

# thres.prob = thres.prob[,2:62]

thres.prob = read.csv("projection_probs_deep.csv", sep = ",", header = T)
thres.prob = thres.prob[,-c(1)]

Z1_node = matrix(NA, nrow = 1, ncol = 140 )
Z2_node = matrix(NA, nrow = 1, ncol = 140 )
# 
# Z1_range = matrix(NA, nrow = 2, ncol = length(thresholds))
# Z2_range = matrix(NA, nrow = 2, ncol = length(thresholds))
# 
for(i in 1:140){
  df_sample = df[df$class == i,]
  Z1_node[,i] = mean(df_sample$proj_var1)
  Z2_node[,i] = mean(df_sample$proj_var2)

  # Z1_range[1,i] = min(df[df$threshold2 == thresholds[i],]$data1)
  # Z1_range[2,i] = max(df[df$threshold2 == thresholds[i],]$data1)
  # 
  # Z2_range[1,i] = min(df[df$threshold2 == thresholds[i],]$data2)
  # Z2_range[2,i] = max(df[df$threshold2 == thresholds[i],]$data2)
}


####### dataset projection explanation through plots ##########
n_model = 7
n_splits_per_model = 20

tau = seq(0.1,0.9,length.out = n_model)
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
abline(a = a[6], b = b[6])
abline(a = a[7], b = b[7])
points(Z1_node,Z2_node, pch = 17, col = "blue", cex = 1.7)



### bivariate conditional cdf
cond.cdf_2d <- function(z,index){
  sum = 0
  h = 0.2
  z1 = z[1]
  z2 = z[2]
  for(i in 1:140){
  #   if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
  #     z1_res = 0
  #     z2_res = 0
  #   }else{
      z1_res = z1 - Z1_node[,i]
      z2_res = z2 - Z2_node[,i]
  # }
    
    p = thres.prob[index,i]
    sum = sum + (p*pbinorm(z1_res,z2_res,var1 = h,var2 = h))
    
  }
  return(sum)
}

cond.pdf_2d <- function(z1,z2,index){
  sum = 0
  h = 0.05
  # Z1_node[,62] = estim.model$pred[1]
  # Z2_node[,62] = estim.model$pred[2]
  for(i in 1:140){
    #   if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #     z1_res = 0
    #     z2_res = 0
    #   }else{
    z1_res = z1 - Z1_node[,i]
    z2_res = z2 - Z2_node[,i]
    # }
    
    p = thres.prob[index,i]
    sum = sum + (p*dbinorm(z1_res,z2_res,var1 = h,var2 = h))
    
  }
  return(sum)
}

z2_given = -4.270493
index1 = 1

marginal.z2 <- function(z2,index){
  sum = 0
  h = 0.2
  for(i in 1:140){
    
      z2_res = (z2 - Z2_node[,i])/h
    p = thres.prob[index,i]
    sum = sum + (p*dnorm(z2_res))/h
  }
  return(sum)
}

marginal.cond.cdf_z1 <- function(z1,z2,index){
  sum = 0
  h = 0.2
  for(i in 1:140){
    #   if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #     z1_res = 0
    #     z2_res = 0
    #   }else{
    z1_res = (z1 - Z1_node[,i])/h
    z2_res = (z2 - Z2_node[,i])/h
    # }
    
    p = thres.prob[index,i]
    sum = sum + (p*pnorm(z1_res)*dnorm(z2_res))/h
    
  }
  
  sum = sum/marginal.z2(z2,index)
  return(sum)
}

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



#### some calculations for pit histogram 

pit_val = rep(NA,8100)
for(i in 1:8100){
  pit_val[i] = marginal.cond.cdf_z1(df$var1[i],df$var2[i],i)
}

##### median absolute deviation #######
median.deep = rep(NA,100)
for(i in 1:100){
  median.deep[i] = marginal.cond.quantile_z1(0.5,df$var2[i],i)
}
mad_deep = mean(abs(median.deep-df$var1[1:100]))

#### run gaussian_estimation.R before running this to get the d1 coordinates
probability_deep = rep(NA,100*100)

for(i in 1:dim(d1)[1]){
  probability_deep[i] = cond.pdf_2d(d1$x[i],d1$y[i], 11)
}


#### Some additional claculations #################

#### run gaussian_estimation.R before running this
### This is for cdf calculation over the region
z1_95_deep = c()
z2_95_deep = c()

for(i in 1:dim(d1)[1]){
  a = cond.cdf_2d(d1$x[i],d1$y[i], 50)
  print(a)
  if( a <= 0.95 & a >= 0.05){
    z1_95_deep = cbind(z1_95_deep,d1$x[i])
    z2_95_deep = cbind(z2_95_deep,d1$y[i])
  }
}

points(z1_95_deep,z2_95_deep,col = "green")

### univariate conditional cdf

### cdf of z1|z2 

marginal.z2 <- function(z2,index){
  sum = 0
  h = 0.12
  for(i in 1:length(thresholds)){
    if(z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
      z2_res = 0
    }else{
      z2_res = z2 - Z2_node[,i]
    }
    p = thres.prob[index,i]
    sum = sum + (p*dnorm(z2, Z2_node[,i] , h))
  }
  return(sum)
}

marginal.z2_cdf <- function(z2,index){
  sum = 0
  h = 0.12
  for(i in 1:length(thresholds)){
    # if(z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #   z2_res = 0
    # }else{
    # sum1 = pnorm((z2 - Z2_range[1,i])/h, 0 , sqrt(emp.var2)) + 
    #   pnorm((z2 - Z2_range[2,i])/h, 0 , sqrt(emp.var2))
    z2_res = z2 - Z2_node[,i]
    # }
    # print(z2_res)
    p = thres.prob[index,i]
    sum = sum + (p*pnorm(z2, Z2_node[,i] , h))
  }
  return(sum)
}

# z1_values = seq(-10,10,length.out = 100)

############
probs = rep(NA,6400)
for(i in 1:6400){
  probs[i] = marginal.z2_cdf(df$data2[i],1)
  # print(probs[i])
  # print(i)
}
# plot(z1_values,probs)
hist(probs)
###############

cond.pdf_z1_z2 <- function(z1,z2,index){
  sum = 0
  h = 1.5
  for(i in 1:length(thresholds)){
    # if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #   z1_res = 0
    #   z2_res = 0
    # }else{
      z1_res = z1 - Z1_node[,i]
      z2_res = z2 - Z2_node[,i]
    # }
    p = thres.prob[index,i]
    sum = sum + (p*dbinorm(z1,z2,Z1_node[,i],Z2_node[,i],h,h,0))
  }
  return(sum/(marginal.z2(z2,index)))
}

cond.cdf_z1_z2 <- function(z1,z21,z22,index){
  sum = 0
  h = 1
  for(i in 1:length(thresholds)){
    # if(z1>=Z1_range[1,i] && z1<=Z1_range[2,i] && z2>=Z2_range[1,i] && z2<=Z2_range[2,i]){
    #   z1_res = 0
    #   z2_res = 0
    # }else{
    z1_res = z1 - Z1_node[,i]
    z2_res = z2 - Z2_node[,i]
    # }
    p = thres.prob[index,i]
    sum = sum + (p*pnorm(z1,Z1_node[,i],h)*(pnorm(z22,Z2_node[,i],h) - pnorm(z21,Z2_node[,i],h)))
  }
  return(sum/(marginal.z2_cdf(z22,index) - marginal.z2_cdf(z21,index)))
}

### integrate the following to get cdf
# integrand <- function(x){
#   return(cond.pdf_z1_z2(x,-1.31,1))
# }


z1_values = seq(-10,10,length.out = 100)
probs = rep(NA,100)
for(i in 1:100){
  probs[i] = cond.pdf_z1_z2(z1_values[i],-1.31,1)
  print(probs[i])
  print(i)
}

integrand <- function(x){
  return(cond.pdf_z1_z2(x,-1.31,1))
}
integrate(integrand,-Inf,100)$value

# df_test = read.csv("test_data.csv",sep = ",", header = T)

probs = rep(NA,3000)#dim(df)[1])



for(i in 1:3000){
  integrand <- function(x){
    return(cond.pdf_z1_z2(x,-1.31,1))
  }
  probs[i] = integrate(integrand,-Inf,df$data1[i])$value
  print(probs[i])
  print(i)
}
hist(probs)

library(pracma)
fun <- function(x, y) cos(x) * cos(y)
integral2(fun, 0.499, 0.511, 0, 1, reltol = 1e-10)



a = rnorm(1000,0,1)
probs = pnorm(a,0,1)








x <- rnorm(50)
y <- runif(30)
# Do x and y come from the same distribution?
a = ks.test(x, y)


