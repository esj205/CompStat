rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(rstudioapi)
library(microbenchmark)
library(ggplot2)
sourceCpp("q2-1.cpp")


#######################
##Problem 2.1
#######################
###################
##Dataset 설정
data <- c(1.77, -0.23, 2.76, 3.8, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29,
            3.71, -2.4, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)
theta = seq(-50, 50, by=0.01)

###################
##(a)newton-raphson method

#log-likelihood graph 그리기
llfunction <- function(theta, data){
  return(sum(dcauchy(data, theta, scale=1, log = TRUE)))
}

ll_result <- sapply(theta, llfunction, data=data)
plot(theta, ll_result, type='l', ylab='log likelihood', xlab=expression(theta))


newton_list <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
for(i in newton_list){
  print(newton(i, 0.0001, data, 10000))
}

newton(mean(data), 0.0001, data, 10000) #data의 평균을 starting point로

theta[which.max(ll_result)] ##graph에서 max값의 theta값


##gprime graph
gprime <- NULL
for(i in 1:length(theta)){
  gprime[i] <- sum(-2*(theta[i]-data)/((data-theta[i])^2+1))
}

gprime <- as.data.frame(gprime)

ggplot(data=gprime,aes(x=theta,y=gprime)) +
  geom_line() + 
  geom_hline(yintercept=0,col="red",linetype="dashed") +
  theme_bw()



###################
##(b)bisection method

bisec(-1,1,0.0001,data, 10000)
bisec(-10,-5,0.0001, data, 10000)



###################
##(c)fixed-point method
alpha_list <- c(3,1,0.64,0.25,0.1)

for(i in alpha_list){
  print(fixed(1, i, 0.0001, data, 10000))
}


for(i in alpha_list){
  print(fixed(-1, i, 0.0001, data, 10000))
}

for(i in alpha_list){
  print(fixed(-5, i, 0.0001, data, 10000))
}


###################
##(d)secant method
secant_0_list <- c(-1,-2,-3)

for(i in secant_0_list){
  print(secant(i,-1,0.0001, data, 10000))
}

for(i in secant_0_list){
  print(secant(i,1,0.0001, data, 10000))
}

for(i in secant_0_list){
  print(secant(i,3,0.0001, data, 10000))
}


###################
##(e)comparison

microbenchmark(newton(-1, 10e-8, data, max_iter = 1000))
microbenchmark(bisec(-1,1,10e-8,data, max_iter = 1000))
microbenchmark(fixed(-1, 0.64, 10e-8, data, max_iter = 1000))
microbenchmark(secant(-1,1,10e-8, data, max_iter = 1000))





#######################
##Problem 2.2
#######################
###################
##Dataset 설정
data2 = c(3.91,4.85,2.28,4.06,3.70,4.04,5.46,3.53,2.28,1.96,2.53,3.88,2.22,3.47,4.82,2.46,2.99,
         2.54,0.52,2.50)
theta2 <- seq(from=-pi,to=pi,length.out = 1000)

sourceCpp('q2-2.cpp')


###################
##(a)Graph the plot
llfunction2 <- function(theta, data){
  return(sum(log((1-cos(data-theta))/(2*pi))))
}

ll_result2 <- sapply(theta2, llfunction2, data=data2)
plot(theta2, ll_result2, type='l', ylab='log likelihood', xlab=expression(theta))


###################
##(b)MoM Estimator
asin(mean(data2)-pi)
mom <- asin(mean(data2)-pi)


###################
##(c)MLE
theta2[which.max(ll_result2)]
newton2(mom, 10e-5, data2, max_iter = 1000)
newton2(-2.7, 10e-5, data2, max_iter= 1000)
newton2(2.7, 10e-5, data2, max_iter= 1000)

###################
##(d)local mode
theta3 <- seq(from=-pi, to=pi, length.out=200)
local_mode <- c()
num=0
for (i in theta3){
  num=num+1
  local_mode[num]<-round(newton2(i, 10e-5, data2, max_iter = 1000),5)
}
local_mode.df <- cbind(theta3, local_mode)
length(unique(local_mode))



###################
##(e)Comparison between two nearly equal starting values
theta2[which.min(ll_result2)]

newton2(x=2.279935-0.0001, epsilon = 0.0001, data = data2, max_iter = 1000)
newton2(x=2.279935+0.0001, epsilon = 0.0001, data = data2, max_iter = 1000)


#######################
##추가숙제 1
#######################
###################
sourceCpp("q2-3.cpp")
newton3(2, 10e-8,10000)
secant3(2,3,10e-8,10000)
bisec3(2,4,10e-8,10000)
fixed3(2,0.7,10e-8,10000)
