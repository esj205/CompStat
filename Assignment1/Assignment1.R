#dataset
mydata <- c(1.77, -0.23, 2.76, 3.8, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29,
            3.71, -2.4, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)
theta = seq(-60, 60, by=0.01)

#log-likelihood function of cauchy
llfunction <- function(theta, data){
  return(sum(dcauchy(data, theta, scale=1, log = TRUE)))
}


ll_result <- sapply(theta, llfunction, data=mydata)
plot(theta, ll_result, type='l', ylab='log likelihood', xlab=expression(theta))


#cauchy의 mle 구하는 function
mlecauchy=function(x,startvalue, toler=.001){      
  n=length(x);
  thetahatcurr=startvalue;
  g.prime=2*sum((x-thetahatcurr)/(1+(x-thetahatcurr)^2))
  while(abs(g.prime)>toler){
    g.2prime=2*sum(((x-thetahatcurr)^2-1)/(1+(x-thetahatcurr)^2)^2);
    thetahatnew=thetahatcurr-g.prime/g.2prime;
    thetahatcurr=thetahatnew;
    g.prime=2*sum((x-thetahatcurr)/(1+(x-thetahatcurr)^2))
  }
  return(thetahatcurr);
}

mlecauchy(mydata, 1.5, 0.0001)


##startingpoints 비교
startingpoints <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
num=0
newton.list <- as.list(NA)
for (i in startingpoints){
  num = num+1
  newton.list[[num]] <- mlecauchy(mydata, i, toler = 0.1)
  names(newton.list)[num] <- as.character(i)
}
newton.list

##mean data 사용
mean(mydata)
mlecauchy(mydata, mean(mydata), toler=0.001)
## 굳이?

##R에 내장되어 있는 함수 사용
optimize(function(theta) -sum(dcauchy(mydata, location=theta, log=TRUE)),  c(-50,50))


###########################################################
###Bisection
###########################################################
#setup
a = -1
b = 1
x = a+(b-a)/2
itr = 40
g.prime <- function(theta){
  2*sum((mydata-theta)/(1+(mydata-theta)^2))
}

#bisection model
for (i in 1:itr){
  if (g.prime(a)*g.prime(x) < 0) {b = x}
  else {a = x}
  x = a+(b-a)/2
}
x

##different starting points
a.2 = -10
b.2 = -5
x = a.2 + (b.2-a.2)/2
for (i in 1:itr){
  if (g.prime(a.2)*g.prime(x) < 0) {b.2 = x}
  else {a.2 = x}
  x = a.2+(b.2-a.2)/2
}
x

#############################################################
###Fixed Points
#############################################################
#setup
x = 2.29
itr = 100
num = 0
fixed.list <- as.list(NA)

# MAIN
for(alpha in c(1, 0.64, 0.25)){
  for(i in 1:itr){x = alpha*g.prime(x) + x}
  num = num + 1
  fixed.list[[num]] <- x
  names(fixed.list)[num]<-as.character(alpha)
}
fixed.list

x = 1
num = 0
fixed.list <- as.list(NA)
for(alpha in c(1, 0.64, 0.25)){
  for(i in 1:itr){x = alpha*g.prime(x) + x}
  num = num + 1
  fixed.list[[num]] <- x
  names(fixed.list)[num]<-as.character(alpha)
}
fixed.list

fixed.method <- function(f, x, alpha, tol=1e-9, n=500) {
  for (i in 1:n){
    x_new <- alpha*f(x) + x
    if(abs(x_new-x)<tol){
      return(x_new)
    }
    x <- x_new
  }
}
##얘는 값이 잘 안 나오는데 아마 tol값이 문제인거 같음.

#############################################################
######Secant Method
############################################################
secant.method <- function(f, x0, x1, tol=1e-9, n=500) {
  for (i in 1:n){
    x2 <- x1 - f(x1) / ((f(x1) - f(x0)) / (x1 - x0))
    if(abs(x2-x1)<tol){
      return(x2)
    }
    x0 <- x1
    x1 <- x2
  }
}

secant.method(g.prime, -2, -1)
secant.method(g.prime, -3, 3)

