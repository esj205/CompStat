########################
##1
########################
set.seed(123)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('final.cpp')
sourceCpp('final2.cpp')
sourceCpp('final3.cpp')
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
##(a)
##functions
h_function <- function(x){
  res <- (cos(50*x) + sin(20*x))^2
  return(res)
}

h_prime <- function(x){
  res <- 2*(cos(50*x)+sin(20*x))*(20*cos(20*x)-50*sin(50*x))
  return(res)
}

h_2prime <- function(x){
  res <- 2*(20*cos(20*x)-50*sin(50*x))^2 +
    2*(-2500*cos(50*x)-400*sin(20*x))*(cos(50*x)+sin(20*x))
  return(res)
}

#Newton-raphson function
Newton <- function(max_iteration, epsilon, x){
  diff = 1
  i = 0
  while(abs(diff) > epsilon){
    diff = -h_prime(x)/h_2prime(x)
    x = x + diff
    i = i + 1
    if(i == max_iteration){
      print("Iteration is over. Need More Max_iter")
      break
    }
    if(x <0 | x>1){
      print('Out of Range. Please Check the Starting value')
      break
    }
  }
  return(x)
}

#h(x) plot
x_seq <- seq(0, 1, length.out = 1000)
h_x <- h_function(x_seq)
plot(x_seq, h_x, type='l')


#Optimization
Newton(10000, 0.001, seq(0.36, 0.38, length.out = 10))
abline(v=0.3791384)
h_function(0.3791384)

Newton(10000, 0.001, seq(0.55, 0.58, length.out = 10))
abline(v=0.5633394)
h_function(0.5633394)

Newton(10000, 0.001, seq(0.72, 0.75, length.out = 10))
abline(v=0.7480271)
h_function(0.7480271)


optimize(h_function, c(0,1), tol=0.0001, maximum = TRUE)


##(b)
#setting
N <- 2500
u <- runif(N)
xval1 <- rep(0, 2500)
r <- 0.5

#annealing function
annealing <- function(start_value, r,iteration){
  run.current <- start_value
  run.best <- run.current
  best <- h_function(start_value)
  runs <- c()
  for(i in 1:iteration){
    run_h <- h_function(run.current)
    u <- runif(1, min = max((run.current-r),0), max=min((run.current+r),1))
    stat <- exp(log(i) * (h_function(u) - h_function(run.current)))
    p.current <- min(stat, 1)
    if(runif(1)<p.current){
      run.current <- u
      run_h <- h_function(u)
    }
    if(run_h > best){
      run.best <- run.current
      best <- run_h
    }
    runs[i] <- run.current
  }
  res <- list(x= runs, max = c(run.best, best))
  return(res)
}

#four starting points
anneal1 <- annealing(0.2, 0.5, 10000)
anneal2 <- annealing(0.4, 0.5, 10000)
anneal3 <- annealing(0.6, 0.5, 10000)
anneal4 <- annealing(0.8, 0.5, 10000)

#maximum values
anneal1$max
anneal2$max
anneal3$max
anneal4$max

#trajectory
plot(anneal1$x, h_function(anneal1$x), type='l', lwd=2, xlab='x', ylab='h(x)', main='Trajectory')
lines(anneal2$x, h_function(anneal2$x), col='2', lwd=2)
lines(anneal3$x, h_function(anneal3$x), col='3', lwd=2)
lines(anneal4$x, h_function(anneal4$x), col='4', lwd=2)


########################
##2
########################
##(b)
#setting
rikz <- read.table('rikz.txt', header=TRUE)
x <- rikz$NAP
y <- rikz$Richness
z <- rikz$Beach
uid <- c(1:45)

#EM Algorithm
EM <- function(y, x, z, uid, iteration){
  y <- as.matrix(y)
  x <- as.matrix(x)
  z <- as.matrix(z)
  N <- length(uid)
  n <- length(y)
  
  ##Initial value
  beta <- as.vector(solve(t(x)%*%x)%*%t(x)%*%y)
  g <- diag(rep(1, ncol(z)))
  sig2 <- 1
  residu <- as.vector(y - x%*%beta)
  
  ##Iteration
  for(j in 1:iteration){
    ##E-step
    P <- 0
    R <- 0
    C <- 0
    mu <- NULL
    u <- NULL
    for(i in uid){
      ##E-step
      xi <- x[i,]
      zi <- z[i,]
      residui <- residu[i]
      gammai <- solve(t(zi)%*%zi/sig2 + solve(g))
      mui <- (gammai%*%t(zi)%*%residui)/sig2
      mu <- c(mu, mui)
      u <- c(u, zi%*%mui)
      si <- gammai + mui %*% t(mui)
      R <- R + si
      P <- P + sum(diag(si%*%t(zi)%*%zi))
      C <- C + t(mui) %*% t(zi) %*% residui
    }
    ##M-step
    beta <- as.vector(solve(t(x) %*% x) %*% t(x) %*% (y-u))
    residu <- as.vector(y-x%*%beta)
    sig2 <- (sum(residu^2)-2*C[1] + P)/n
    g <- as.matrix(R/N)
  }
  return(list(beta=beta, g=g, sigma2 = sig2))
}


EM(y,x,z,uid, 10000)


########################
##3
########################
##data
postdf <- read.csv('output.csv')
clinic_df <- read.csv('data.csv')

#sk에 x로 포함되는 경우.
trape_list1 <- c()
simpson_list1 <- c()
for(i in 1:100){
  trape_list1[i] <- trape(clinic_df[i,1], clinic_df[i,2], 1000)
  simpson_list1[i] <- simpson(clinic_df[i,1], clinic_df[i,2], 1000)
}

list1_df <- data.frame(trape = trape_list1, simpson = simpson_list1)

#sk에 xi로 포함되는 경우.
trape_list2 <- c()
simpson_list2 <- c()
for(i in 1:100){
  trape_list2[i] <- trape2(clinic_df[i,4],clinic_df[i,1], clinic_df[i,2], 1000)
  simpson_list2[i] <- simpson2(clinic_df[i,4],clinic_df[i,1], clinic_df[i,2], 1000)
}

list2_df <- data.frame(trape = trape_list2, simpson = simpson_list2)

list1_df
list2_df


########################
##4
########################
##data
uscancer <- readLines('UScancer.txt')
uscancer <- as.data.frame(uscancer)
dim(uscancer)
length(strsplit(as.character(uscancer[1,]), split='')[[1]])

cancer <- matrix(NA, nrow=58, ncol = 66)
for(i in 1:58){
  charc <- strsplit(as.character(uscancer[i,]), split='')[[1]]
  num <- as.numeric(charc)
  cancer[i,] <- num
}

cancer[cancer==0] <- -1
cancer[cancer==2] <- 0


##bootstrapping
mple_df <- bootdf(1000, -0.3205, 0.1115, cancer)
dmh_df <- bootdf(1000, -0.3030, 0.1227, cancer)
aex_df <- bootdf(1000, -0.3017, 0.1224, cancer)


#bootstrap 결과 불러오기
mple_df <- read.csv('mple_df.csv', header=TRUE)
mple_df <- mple_df[,-1]
dmh_df <- read.csv('dmh_df.csv', header=TRUE)
dmh_df <- dmh_df[,-1]
aex_df <- read.csv('aex_df.csv', header=TRUE)
aex_df <- aex_df[,-1]



##RMSE 구하기
s1_obs <- sum(cancer)
s2_obs <- neigh_sum(cancer)

mple_s1_bias2 <- sum((mple_df$s1-s1_obs)^2)/1000
mple_s1_var <- sum((mple_df$s1-mean(mple_df$s1))^2)/1000
mple_s1 <- sqrt(mple_s1_bias2+mple_s1_var)

dmh_s1_bias2 <- sum((dmh_df$s1-s1_obs)^2)/1000
dmh_s1_var <- sum((dmh_df$s1-mean(dmh_df$s1))^2)/1000
dmh_s1 <- sqrt(dmh_s1_bias2+dmh_s1_var)

aex_s1_bias2 <- sum((aex_df$s1-s1_obs)^2)/1000
aex_s1_var <- sum((aex_df$s1-mean(aex_df$s1))^2)/1000
aex_s1 <- sqrt(aex_s1_bias2+aex_s1_var)

mple_s2_bias2 <- sum((mple_df$s2-s2_obs)^2)/1000
mple_s2_var <- sum((mple_df$s2-mean(mple_df$s2))^2)/1000
mple_s2 <- sqrt(mple_s2_bias2+mple_s2_var)

dmh_s2_bias2 <- sum((dmh_df$s2-s2_obs)^2)/1000
dmh_s2_var <- sum((dmh_df$s2-mean(dmh_df$s2))^2)/1000
dmh_s2 <- sqrt(dmh_s2_bias2+dmh_s2_var)

aex_s2_bias2 <- sum((aex_df$s2-s2_obs)^2)/1000
aex_s2_var <- sum((aex_df$s2-mean(aex_df$s2))^2)/1000
aex_s2 <- sqrt(aex_s2_bias2+aex_s2_var)


rmse_df <- data.frame(mple = c(mple_s1, mple_s2), dmh = c(dmh_s1, dmh_s2), aex = c(aex_s1, aex_s2))
rownames(rmse_df) <- c('S1', 'S2')
rmse_df

