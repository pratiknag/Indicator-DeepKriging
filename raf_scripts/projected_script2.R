library(fields)
library(geoR)
library(ggplot2)

# var1_probs = rep(NA,10)
# var2_probs = rep(NA,10)
# 
# true_05 = matrix(NA,nrow = 2, ncol = 10)
# true_25 = matrix(NA,nrow = 2, ncol = 10)
# true_50 = matrix(NA,nrow = 2, ncol = 10)
# true_75 = matrix(NA,nrow = 2, ncol = 10)
# true_95 = matrix(NA,nrow = 2, ncol = 10)
# 
# estim_05 = matrix(NA,nrow = 2, ncol = 10)
# estim_25 = matrix(NA,nrow = 2, ncol = 10)
# estim_50 = matrix(NA,nrow = 2, ncol = 10)
# estim_75 = matrix(NA,nrow = 2, ncol = 10)
# estim_95 = matrix(NA,nrow = 2, ncol = 10)

count = 1
# location = 10
for(location in 150:160){
setwd("/Users/nagp/Desktop/env stats/project/")
exa_data = read.csv("train_data.csv", sep = ",", header = T)


test_location = location

df = do.call(rbind, Map(data.frame,x = exa_data$x, y = exa_data$y,
                        var1=exa_data$var1 - mean(exa_data$var1),
                        var2 = exa_data$var2 - mean(exa_data$var2)))

mean1 = mean(exa_data$var1)
mean2 = mean(exa_data$var2)
un.grd.train = df[sample(1:dim(df)[1],1080),]

# exa_data = read.csv("test_data.csv", sep = ",", header = T)
exa_data = read.csv("median_results.csv", sep = ",", header = T)
test.data = do.call(rbind, Map(data.frame,x = exa_data$x, y = exa_data$y, 
                               var1=exa_data$var1,var2 = exa_data$var2))

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
  # test.index = c(1:dim(test.data)[1],(dim(un.grd.train)[1]+1):(dim(un.grd.train)[1]+dim(test.data)[1]))
  test.index = c(1,1082)
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

estim.model$pred <- estim.model$pred + c(mean1,mean2)

# var1_probs[count] = pnorm(Z90[1,test_location],mean = estim.model$pred[1],estim.model$conditional_var[1,1]) - 
#   pnorm(Z10[1,test_location],mean = estim.model$pred[1],estim.model$conditional_var[1,1])
# 
# var2_probs[count] = pnorm(Z90[2,test_location],mean = estim.model$pred[2],estim.model$conditional_var[2,2]) - 
#   pnorm(Z10[2,test_location],mean = estim.model$pred[2],estim.model$conditional_var[2,2])
# 
# 
# true_10[1,count] = qnorm(0.1,estim.model$pred[1],estim.model$conditional_var[1,1])
# 
# true_90[1,count] = qnorm(0.9,estim.model$pred[1],estim.model$conditional_var[1,1])
# 
# true_10[2,count] = qnorm(0.1,estim.model$pred[2],estim.model$conditional_var[2,2])
# 
# true_90[2,count] = qnorm(0.9,estim.model$pred[2],estim.model$conditional_var[2,2])
# 
# estim_10[,count] = Z10[,test_location]
# estim_90[,count] = Z90[,test_location]
# print("###############")
# print(estim.model$pred)
# print("###############")
# print(Z90[,test_location])
# print(Z10[,test_location])
# print("###############")

val1 = seq(test.data$var2-1.5,test.data$var2+1.5,length.out = 100)

probability_gaussian = rep(NA,100)
probability_deep = rep(NA,100)
probability_indep = rep(NA,100)
for(i in 1:100){
  probability_gaussian[i] = pnorm(val1[i],estim.model$pred[2],estim.model$conditional_var[2,2])
  probability_deep[i] = marginal.cdf_var2(val1[i],test_location)
  probability_indep[i] = indep.cdf_var2(val1[i],test_location)
}




# point_true = rep(NA,5)
# 
# point_estim = rep(NA,5)
# 
# point_true[1] = qnorm(0.05,estim.model$pred[1],estim.model$conditional_var[1,1])
# point_true[2] = qnorm(0.25,estim.model$pred[1],estim.model$conditional_var[1,1])
# point_true[3] = qnorm(0.5,estim.model$pred[1],estim.model$conditional_var[1,1])
# point_true[4] = qnorm(0.75,estim.model$pred[1],estim.model$conditional_var[1,1])
# point_true[5] = qnorm(0.95,estim.model$pred[1],estim.model$conditional_var[1,1])
# 
# point_estim[1] = Z05[1,test_location]
# point_estim[2] = Z25[1,test_location]
# point_estim[3] = Z50[1,test_location]
# point_estim[4] = Z75[1,test_location]
# point_estim[5] = Z95[1,test_location]
# 
# probs = c(0.05,0.25,0.50,0.75,0.95)

df2 <- data.frame(model=rep(c("true", 
                              "estimated_biv","estimated_indep"), each=100),
                  N=rep(val1,3),
                  time=c(probability_gaussian,probability_deep,probability_indep))

p<-ggplot(df2, aes(x=N, y=time, group=model)) +
  geom_line(aes(linetype=model,color=model),size=3)+
  geom_point(aes(shape=model,color=model),size=5)+
  labs(y= "probability", x = "variable 2", size = 7)+ 
  theme(legend.text=element_text(size=rel(1.5)), axis.title=element_text(size=22,face="bold"),
        axis.text=element_text(size=17))
print(p)
# plot(point_true,probs, type = "l", col = "red")
# points(point_estim,probs, type = "l", col = "green")

# z1_95 = c()
# z2_95 = c()

# val1 = seq(test.data$var1-3,test.data$var1+3,length.out = 100)
# val2 = seq(test.data$var2-3,test.data$var2+3,length.out = 100)
# 
# # val1 = seq(test.data$data1-3,test.data$data1+3,length.out = 100)
# # val2 = seq(test.data$data2-3,test.data$data2+3,length.out = 100)
# 
# d1 = expand.grid(x=val1,y=val2)
# 
# probability_gaussian = rep(NA,100*100)
# 
# for(i in 1:dim(d1)[1]){
#   probability_gaussian[i] = mvtnorm::dmvnorm(c(d1$x[i],d1$y[i]), 
#                                              mean = c(estim.model$pred) , 
#                                              sigma = estim.model$conditional_var)
# }
# 
# #### run gaussian_estimation.R before running this to get the d1 coordinates
# probability_deep = rep(NA,100*100)
# 
# for(i in 1:dim(d1)[1]){
#   probability_deep[i] = cond.pdf_2d(d1$x[i],d1$y[i], test_location)
# }
# 
# ### # run bivarite_cdf.R pdf probability generation first to get the probability_deep variable
# contour_df = data.frame(x=d1$x,y=d1$y,prob_G=probability_gaussian,prob_D = probability_deep)
# 
# p = ggplot(contour_df,aes(x = d1$x, y = d1$y, z = prob_G)) +
#   stat_contour(color = 'red') +
#   stat_contour(data = contour_df, aes(x = d1$x, y = d1$y, z = prob_D), color = 'green')
# 
# print(p)
count = count + 1
}

df = data.frame(x = 1:10,t_10=true_10[1,],t_90=true_90[1,],e_10=estim_10[1,],e_90=estim_90[1,])

ggplot(df, aes(x)) +                    # basic graphical object
  geom_line(aes(y=t_10), colour="red") +  # first layer
  geom_line(aes(y=e_10), colour="green")+
  geom_line(aes(y=t_90), colour="red")+# second layer
  geom_line(aes(y=e_90), colour="green")

w = mvtnorm::rmvnorm(100, mean = c(estim.model$pred) , sigma = estim.model$conditional_var)



