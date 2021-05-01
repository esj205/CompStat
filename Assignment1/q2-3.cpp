// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace std;
using namespace Rcpp;

double g(const double x){
  double value;
  value = log(x)/(1+x);
  return value;
}


double gprime1(double x){
  double value;
  value = (1+x-x*log(x))/(x*pow((1+x),2));
  return value;
}

double gprime2(double x){
  double value;
  value = (-log(x)*x*pow((1+x),2) - (1+x-x*log(x))*(3*pow(x,2)+4*x+1)  ) /
    (pow((x*pow((1+x),2)),2));
  return value;
}

// [[Rcpp::export]]
double newton3(double x, double epsilon, int max_iter){
  int i = 0;
  double diff = 1.0;
  
  while(fabs(diff) > epsilon){
    diff = -gprime1(x)/gprime2(x);
    x = x + diff;
    i = i + 1;
    if(i==max_iter){
      cout << "Iteration is over. Need more max_iter. \n"<<endl;
      break;
    }
  }
  
  return x;
}

// [[Rcpp::export]]
double bisec3(double a, double b, double epsilon, int max_iter){
  int i = 0;
  double diff = 1.0;
  double c = (a+b)/2; 
  
  while(fabs(diff) > epsilon){
    if (gprime1(b)*gprime1(c) < 0) 
      a = c ;
    else 
      b = c ;
    c=(a+b)/2;
    diff=fabs(b-a);
    i = i + 1;
    if(i==max_iter){
      cout << "Iteration is over. Need more max_iter. \n"<<endl; 
      break;
    }
  }
  
  return c;
}

// [[Rcpp::export]]
double fixed3(double x0, double alpha, double epsilon, int max_iter){
  int i = 0;
  double diff = 1.0;
  double x_n = 0.0;
  
  while(fabs(diff) > epsilon){
    x_n = x0 + alpha*gprime1(x0);
    diff = fabs(x_n-x0);
    x0 = x_n;
    i = i + 1;
    if(i==max_iter){
      cout << "Iteration is over. Need more max_iter. \n"<<endl;
      break;
    }
  }
  return x0;
}

// [[Rcpp::export]]
double secant3(double x0, double x1, double epsilon, int max_iter){
  int i = 0;
  double diff = 1.0;
  double x2 = 1.0;
  
  while(fabs(diff) > epsilon){
    x2 = x1 - gprime1(x1)*(x1-x0)/(gprime1(x1)-gprime1(x0));
    diff = fabs(x2-x1);
    x0 = x1;
    x1 = x2; 
    i = i + 1;
    if(i==max_iter){
      cout << "Iteration is over. Need more max_iter. \n"<<endl;
      break;
    }
  }
  
  return x1;
}

