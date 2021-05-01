library(Rcpp)
library(RcppArmadillo)
library(rstudioapi)
sourceCpp('HW3.cpp')

###Multinomial EM Algorithm
x = c(38,34,125)
n = rep(0,4)
theta = 1/2
p = c(1/2-theta/2, theta/4, theta/4, 1/2)
itr = 50

#E-step
e.step = function(x, p){
  n1 = 38
  n2 = 34
  n3 = 125*(p[3]/(p[4]+p[3]))
  n= c(n1,n2,n3)
  return(n)
}


#M-step
m.step =function(n){
  theta = (34+n[3])/(72+n[3])
  p = c(1/2-theta/2, theta/4, theta/4, 1/2)
  return(p)
}


#Iteration
for(i in 1:itr){
  n = e.step(x,p)
  p = m.step(n)
}

#output
p
sum(p)

#RCPP
x = c(38,34,125)
n = rep(0,4)
theta = 1/2
p = c(1/2-theta/2, theta/4, theta/4, 1/2)
itr = 50


multi_out(x,p, itr)


###PEPPERED MOTHS EM Algorithm
x = c(85, 196, 341) #observed data (carbonaria, insularia, typica)
n = rep(0,6) #expected number of each phenotype (CC, CI, CT, II, IT, TT)
p = rep(1/3,3) #probabilities of allele(carbonaria, insularia, typica)
itr = 50 

#E-Step
e.step = function(x,p){
  n.cc = (x[1]*(p[1]^2))/((p[1]^2)+2*p[1]*p[2]+2*p[1]*p[3])
  n.ci = (2*x[1]*p[1]*p[2])/((p[1]^2)+2*p[1]*p[2]+2*p[1]*p[3])
  n.ct = (2*x[1]*p[1]*p[3])/((p[1]^2)+2*p[1]*p[2]+2*p[1]*p[3])
  n.ii = (x[2]*(p[2]^2))/((p[2]^2)+2*p[2]*p[3])
  n.it = (2*x[2]*p[2]*p[3])/((p[2]^2)+2*p[2]*p[3])
  n = c(n.cc,n.ci,n.ct,n.ii,n.it,x[3])
  return(n)
}

#M-step
m.step = function(x,n){
  p.c = (2*n[1]+n[2]+n[3])/(2*sum(x))
  p.i = (2*n[4]+n[5]+n[2])/(2*sum(x))
  p.t = (2*n[6]+n[3]+n[5])/(2*sum(x))
  p = c(p.c,p.i,p.t)
  return(p)
}

## Iteration
for(i in 1:itr){
  n = e.step(x,p)
  p = m.step(x,n)
}

## OUTPUT
p 
sum(p)


##RCPP
x = c(85, 196, 341)
n = rep(0,6)
p = rep(1/3,3)
itr = 50 

pep_out1(x, p, itr)

###4.1(a)
#observed data (carbonaria, insularia, typica, in+ty)
x = c(85, 196, 341, 578) 
#expected number of each phenotype (CC, CI, CT, II, IT, TT)
n = rep(0,6) 
p = rep(1/3,3) #probabilities of allele(carbonaria, insularia, typica)
itr = 30

#E-Step
e.step = function(x,p){
  n.cc = (x[1]*(p[1]^2))/((p[1]^2)+2*p[1]*p[2]+2*p[1]*p[3])
  n.ci = (2*x[1]*p[1]*p[2])/((p[1]^2)+2*p[1]*p[2]+2*p[1]*p[3])
  n.ct = (2*x[1]*p[1]*p[3])/((p[1]^2)+2*p[1]*p[2]+2*p[1]*p[3])
  n.ii = (x[2]*(p[2]^2))/((p[2]^2)+2*p[2]*p[3]) + ((x[4])*(p[2]^2))/((p[2]^2)+2*p[2]*p[3]+p[3]^2)
  n.it = (2*x[2]*p[2]*p[3])/((p[2]^2)+2*p[2]*p[3]) + (2*x[4]*p[2]*p[3])/((p[2]^2)+2*p[2]*p[3]+p[3]^2)
  n.tt = x[3] + ((x[4])*(p[3]^2))/((p[2]^2)+2*p[2]*p[3]+p[3]^2)
  n = c(n.cc,n.ci,n.ct,n.ii,n.it,n.tt)
  return(n)
}

#M-step
m.step = function(x,n){
  p.c = (2*n[1]+n[2]+n[3])/(2*sum(x))
  p.i = (2*n[4]+n[5]+n[2])/(2*sum(x))
  p.t = (2*n[6]+n[3]+n[5])/(2*sum(x))
  p = c(p.c,p.i,p.t)
  return(p)
}

## Iteration
for(i in 1:itr){
  n = e.step(x,p)
  p = m.step(x,n)
}

## OUTPUT
p 
sum(p)


##RCPP
x = c(85, 196, 341, 578) 
n = rep(0,6) 
p = rep(1/3,3)
itr = 30

pep_out2(x, p, itr)


##4.2
#data
obs = c(379.0,299.0,222.0,145.0,109.0,95.0,73.0,59.0,45.0,30.0,24.0,12.0,4.0,2.0,0.0,1.0,1.0)
itr = 100
alpha = 0.5
beta = 0.2
mu = 1.0
lambda = 5.0
param = c(alpha, beta, mu, lambda)

#RCPP
hiv_out(param, obs, itr)




