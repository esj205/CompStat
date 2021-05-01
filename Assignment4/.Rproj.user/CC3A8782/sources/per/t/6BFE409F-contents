##5.3
library(Rcpp)
library(RcppArmadillo)
sourceCpp('HW4.cpp')


#data
x = c(6.52, 8.32, 0.31, 2.82, 9.96, 0.14, 9.64)
x_mean <- mean(x)
sd <- sqrt(9/7)
interval <- seq(0,10, length.out = 100)
norm_dens <- dnorm(interval, mean=x_mean, sd=sd)
plot(interval, norm_dens, type='l')

interval <- c(0,10)

1/riemann_cpp(interval, 100, mean(x), sd)
1/trapezoidal_cpp(interval, 100, mean(x), sd)
1/simpsons_cpp(interval, 100, mean(x), sd)


interval <- c(2,8)

riemann_df(interval, mean(x), sd, 0.0001)
trapezoidal_df(interval, mean(x), sd, 0.0001)
simpsons_df(interval, mean(x), sd, 0.0001)



#(c) 
#ignoring
interval <- c(exp(3)/(1+exp(3)), 1)
trapezoidal_cpp3(interval, 100, x_mean, sd)

#(d)
interval <- c(0, 1/3)
trapezoidal_cpp4(interval, 100, x_mean, sd)


#5.4
round(triangle(5, 6),7)
log(5)








