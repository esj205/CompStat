// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


double g(const double x, arma::dvec data){
  double value;
  value = sum(-log(datum::pi)-log(pow((data-x),2) + 1));
  return value;
}

double gprime1(const double x, arma::dvec data){
  double value;
  value = 2*sum((data-x)/(1+arma::pow((data-x),2)));
  return value;
}

double gprime2(const double x, arma::dvec data){
  double value;
  value = 2*sum((arma::pow((data-x),2)-1)/arma::pow((1+arma::pow((data-x),2)),2));
  return value;
}

double g2(const double x, arma::dvec data){
  double value;
  value = -10*log(2*(datum::pi))-sum(data-x)/2
  return value;
}

double g2(const double x, arma::dvec data){
  double value;
  value = -10*log(2*(datum::pi))-sum(data-x)/2
  return value;
}



// [[Rcpp::export]]
double newton(double x, double epsilon, arma::dvec data, int max_iter){
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

// [[Rcpp::export]]
double bisec(double a, double b, double epsilon, arma::dvec data, int max_iter){
  int i = 0;
  double diff = 1.0;
  double c = (a+b)/2; 
  
  while(fabs(diff) > epsilon){
    if (gprime1(b, data)*gprime1(c, data) < 0) 
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
double fixed(double x0, double alpha, double epsilon, arma::dvec data, int max_iter){
  int i = 0;
  double diff = 1.0;
  double x_n = 0.0;
  
  while(fabs(diff) > epsilon){
    x_n = x0 + alpha*gprime1(x0, data);
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
double secant(double x0, double x1, double epsilon, arma::dvec data, int max_iter){
  int i = 0;
  double diff = 1.0;
  double x2 = 1.0;
  
  while(fabs(diff) > epsilon){
    x2 = x1 - gprime1(x1, data)*(x1-x0)/(gprime1(x1, data)-gprime1(x0, data));
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
