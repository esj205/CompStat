library("Rcpp")
gamma = 'NumericVector Ex_gamma(int n, int a, int b){
  NumericVector x = rgamma(n, a, b);
  return x;
}'
cppFunction(gamma)

gen = Ex_gamma(10000,1,1)
hist(gen, nclass = 100)
