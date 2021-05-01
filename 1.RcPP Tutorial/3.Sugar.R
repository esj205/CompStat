library('Rcpp')
sugar <- 'List sugar_Ex (NumericVector x, NumericVector y){
NumericVector soma = x + y, res = x - y ;
NumericVector prod = x * y, div = x / y ;
LogicalVector menor = x<y, maior = x>y ;
LogicalVector igual = x==y, dif = x!=y ;
return List::create(soma, div, menor, dif);
}
'

cppFunction(sugar)
x = c(8,9,7,6)
y = c(1,7,2,3)
sugar_Ex(x, y)


gamma <- 'NumericVector Ex_gamma(int n, int a, int b){
NumericVector x = rgamma(n,a,b);
return x;
}
'

cppFunction(gamma)
gen = Ex_gamma(1000, 1, 1)
hist(gen, nclass = 100)