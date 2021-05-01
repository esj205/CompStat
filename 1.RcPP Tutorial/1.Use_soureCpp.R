library(Rcpp)
sourceCpp('example1.cpp')
a = 1:4
output1 = 2*a
output2 = timesTwo(a)