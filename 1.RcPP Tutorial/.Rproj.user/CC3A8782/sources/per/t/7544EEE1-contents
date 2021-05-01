library('Rcpp')
src1 <- 'NumericVector example1(NumericVector x){
NumericVector y;
y = 2*x;
return(y);
}
'

cppFunction(src1)
example1(1:4)