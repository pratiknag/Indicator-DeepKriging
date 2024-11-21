library(fields)
library(geoR)
library(ggplot2)
library(numDeriv)
library(calculus)

setwd("/Users/nagp/Desktop/env stats/project/Final_scripts/")
# exa_data = read.csv("train_data.csv", sep = ",", header = T)
test_location = c(1:100)

###################################################################################
# exa_data = df[c(1:(test_location-1),(test_location+1):8100),]
exa_data = df[c(101:8100),]

df1 = do.call(rbind, Map(data.frame,x = exa_data$x, y = exa_data$y,
                        var1=exa_data$var1,
                        var2 = exa_data$var2))

# mean1 = mean(exa_data$data1)
# mean2 = mean(exa_data$data2)
rownames(df1) <- 1:nrow(df1)
un.grd.train = df1[sample(1:8000,1080),]

# exa_data = read.csv("test_data.csv", sep = ",", header = T)
test.data = do.call(rbind, Map(data.frame,x = df$x, y = df$y, 
                               var1=df$var1,var2 = df$var2))

test.data = test.data[test_location,]

# test.data = test.data[34,]
# test.data = df_threshold2[1,]

True.model.estimation<-function()
{
  ########## Test set #############
  full.data.coords = rbind(as.matrix(test.data[,c(1,2)]),as.matrix(un.grd.train[,c(1,2)]))
  m = as.matrix(dist(full.data.coords))
  
  a1<-0.1
  nu1<-0.8
  sigma1<-1
  a2<-0.2
  nu2<-0.9
  sigma2<-1.3
  b11<-1.28
  b12<-0.36
  b21<-0.63
  b22<-0.9
  nug1<-0
  nug2<-0
  ######## Putting hard constraints on the parameters #############
  #n<-nrow(dist.mat)
  C11<-sigma1*matern(m,a1,nu1)
  C22<-sigma2*matern(m,a2,nu2)
  
  COV11<-(b11^2)*C11+(b12^2)*C22
  COV22<-(b21^2)*C11+(b22^2)*C22
  COV12<-(b11*b21)*C11+(b12*b22)*C22
  
  NUG1<-diag(nug1,nrow = nrow(COV11),ncol=ncol(COV11))
  NUG2<-diag(nug2,nrow = nrow(COV22),ncol=ncol(COV22))
  ##################################################################################
  ############## Creating full covariance matrix ###################################
  ##################################################################################
  
  ############## Inverting C11 ##########
  C<-rbind(cbind(COV11+NUG1,COV12),cbind(t(COV12),COV22+NUG2))
  print(dim(C))
  test.index = c(1:dim(test.data)[1],(dim(un.grd.train)[1]+1):(dim(un.grd.train)[1]+dim(test.data)[1]))
  # test.index = c(1:120,1201:1320)
  # test.index = c(1,1082)
  # test.index = c(1:1190,1201:2390)
  ##################################################################################
  ############## Creating full covariance matrix ###################################
  ##################################################################################
  C.test<-C[test.index,test.index]
  C.train<-C[-test.index,-test.index]
  C.test.train<-C[test.index,-test.index]
  #coh12<-lmc.coh(w=u,a1=a1,nu1=nu1,a2=a2,nu2=nu2,b11=b11,b12=b12,b21=b21,b22=b22)
  ############## Inverting C11 ##########
  #C.train<-C22
  #C.test<-C11
  #C.test.train<-C12 print(dim(C.test))
  print(dim(C.train))
  print(dim(C.test.train))
  
  C.inv = solve(C.train)
  
  
  
  prediction<-C.test.train%*%C.inv%*%c(un.grd.train$var1,un.grd.train$var2) #Conditional mean of a Multivariate Gaussian
  cond_var_MAT<- C.test - C.test.train%*%C.inv%*%t(C.test.train)
  
  return(list(pred = prediction, conditional_var = cond_var_MAT ))
  
  
}

estim.model <- True.model.estimation()

##### conditional distribution of Z1|Z2

cond.meanZ1 = estim.model$pred[1:100] + 
  estim.model$conditional_var[101:200,1:100]%*%
  solve(estim.model$conditional_var[101:200,101:200])%*%
  (df$var2[1:100] - estim.model$pred[101:200])

cond.varZ1 = estim.model$conditional_var[1:100,1:100] - 
  estim.model$conditional_var[101:200,1:100]%*%
  solve(estim.model$conditional_var[101:200,101:200])%*%
  estim.model$conditional_var[1:100,101:200]


median.gaussian = rep(NA,100)

for(i in 1:100){
  median.gaussian[i] = qnorm(0.5,cond.meanZ1[i],sqrt(cond.varZ1[i,i]))
}
median.gaussian[42] = 10

mad_g = mean(abs(median.gaussian-df$var1[1:100]))





###############################################################################
val1 = seq(test.data$var1-3,test.data$var1+3,length.out = 100)
val2 = seq(test.data$var2-3,test.data$var2+3,length.out = 100)

# val1 = seq(test.data$data1-3,test.data$data1+3,length.out = 100)
# val2 = seq(test.data$data2-3,test.data$data2+3,length.out = 100)

d1 = expand.grid(x=val1,y=val2)

probability_gaussian = rep(NA,100*100)

for(i in 1:dim(d1)[1]){
  probability_gaussian[i] = mvtnorm::dmvnorm(c(d1$x[i],d1$y[i]), 
                              mean = c(estim.model$pred) , 
                                             sigma = estim.model$conditional_var)
}

### bivariate conditional cdf
index = test_location

#### run gaussian_estimation.R before running this to get the d1 coordinates
probability_deep = rep(NA,100*100)

for(i in 1:dim(d1)[1]){
  probability_deep[i] = cond.pdf_2d(d1$x[i],d1$y[i], index)
}


##########################################################################################

### # run bivarite_cdf.R pdf probability generation first to get the probability_deep variable
contour_df = data.frame(x=d1$x,y=d1$y,prob_G=probability_gaussian,prob_D = probability_deep)

ggplot(contour_df,aes(x = d1$x, y = d1$y, z = prob_G)) +
  stat_contour(color = 'red') +
  stat_contour(data = contour_df, aes(x = d1$x, y = d1$y, z = prob_D), color = 'green')





for(i in 1:dim(d1)[1]){
  a = mvtnorm::pmvnorm(c(-Inf,-Inf),c(d1$x[i],d1$y[i]), mean = c(estim.model$pred),
                   sigma = estim.model$conditional_var)
  print(a[1])
  if( a[1] <= 0.95 & a[1] >= 0.05){
    z1_95 = cbind(z1_95,d1$x[i])
    z2_95 = cbind(z2_95,d1$y[i])
  }
}

  z1_95 = c(z1_95)
  z2_95 = c(z2_95)
  
  plot(d1$x,d1$y)
points(z1_95,z2_95,col = "red")
points(test.data$var1,test.data$var2, col = "blue", type = "p")

data = read.csv("synthetic_data_simulations/training_data/2D_nonGaussian_1200_projection_1train.csv", header = T)



var1 = seq(-2,2,0.1)
g = 0.5
h = 0.5
Tw = ((exp(g*var1)-1)/g) * exp(h*(var1^2)/2)
plot(var1,Tw,type = "l")





