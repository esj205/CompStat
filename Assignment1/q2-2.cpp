// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace std;
using namespace Rcpp;

double g(const double x, arma::dvec data){
  double value;
  value = sum(log(1-cos(data-x)) - log(2*datum::pi));
    return value;
}


double gprime1(double x, NumericVector data){
  double value;
  value = sum(-sin(data-x)/(1-cos(data-x)));
  return value;
}

double gprime2(double x, NumericVector data){
  double value;
  value = sum((cos(data-x)*(1-cos(data-x)) - sin(data-x)*sin(data-x))/
    ((1-cos(data-x))*(1-cos(data-x))));
  return value;
}


// [[Rcpp::export]]
double newton2(double x, double epsilon, NumericVector data, int max_iter){
  int i = 0;
  double diff = 1.0;
  
  while(fabs(diff) > epsilon){
    diff = -gprime1(x, data)/gprime2(x, data);
    x = x + diff;
    i = i + 1;
    if(i==max_iter){
      cout << "Iteration is over. Need more max_iter. \n"<<endl;
      break;
    }
  }
  
  return x;
}

