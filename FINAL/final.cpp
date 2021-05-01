// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;
using namespace std;
using namespace Rcpp;

double c1(double w, double a){
  double gam_11 = -0.1987954;
  double gam_12 = -0.6744738;
  double gam_13 = 0.1701579;
  double res = exp(gam_11 * w + gam_12 * a + gam_13 * w * a);
  return(res);
}

double c2(double w, double a){
  double gam_21 = -0.2518877;
  double gam_22 = 0.3991342;
  double gam_23 = 0.25158;
  double res = exp(gam_21 * w + gam_22 * a + gam_23 * w * a);
  return(res);
}

double f(double x, double eta1, double eta2, double w, double a){
  double res1 = eta1*c1(w,a) + eta2*c2(w,a);
  double res2 = exp(-res1*x);
  return(res2);
}


double trapezoidal_cpp(NumericVector interval, 
                       double n, 
                       double eta1, 
                       double eta2,
                       double w,
                       double a){
  vec x(n-1);
  vec f_list(n-1);
  double b = interval[0];
  double c = interval[1];
  double h = (c-b)/n;
  for(int i = 1; i < n; i++){
    x[i-1] = b + i*h;
  }
  double low = h/2*f(b, eta1, eta2, w, a);
  double high = h/2*f(c, eta1, eta2, w, a);
  for(int j = 1; j < n; j++){
    f_list[j-1] = f(x[j-1], eta1, eta2, w, a);
  }
  double middle = h*accu(f_list);
  double res = low + high + middle;
  return(res);
}

double first_function(double w, double a, double n){
  double eta_11 = 0.07116282;
  double eta_12 = 0.07766563;
  double eta_21 = 0.02380466;
  double eta_22 = 0.02865413;
  NumericVector inter1 = NumericVector::create(0.0, 1.5);
  NumericVector inter2 = NumericVector::create(1.5,3.0);
  double part1 = eta_11 * c1(w,a) * 
    trapezoidal_cpp(inter1, n, eta_11, eta_21, w, a);
  double part2 = eta_12 * c1(w,a) *
    trapezoidal_cpp(inter2, n, eta_12, eta_22, w, a);
  double res = (part1 + part2) * 0.5;
  return(res);
}

double second_function(double w, double a, double n){
  double eta_11 = 0.07116282;
  double eta_12 = 0.07766563;
  double eta_21 = 0.02380466;
  double eta_22 = 0.02865413;
  NumericVector inter1 = NumericVector::create(0.0, 1.5);
  NumericVector inter2 = NumericVector::create(1.5, 3.0);
  double part1 = eta_21 * c2(w,a) * 
    trapezoidal_cpp(inter1, n, eta_11, eta_21, w, a);
  double part2 = eta_22 * c2(w,a) *
    trapezoidal_cpp(inter2, n, eta_12, eta_22, w, a);
  double res = (part1 + part2) * 5.0;
  return(res);
}

double third_function(double w, double a, double n){
  double eta_13 = 0.1052774;
  double eta_14 = 0.1061366;
  double eta_23 = 0.03215047;
  double eta_24 = 0.03584044;
  NumericVector inter1 = NumericVector::create(3.0, 4.5);
  NumericVector inter2 = NumericVector::create(4.5, 6.0);
  double part1 = eta_13 * c1(w,a) * 
    trapezoidal_cpp(inter1, n, eta_13, eta_23, w, a);
  double part2 = eta_14 * c1(w,a) *
    trapezoidal_cpp(inter2, n, eta_14, eta_24, w, a);
  double res = (part1 + part2) * 10.0;
  return(res);
}

double fourth_function(double w, double a, double n){
  double eta_13 = 0.1052774;
  double eta_14 = 0.1061366;
  double eta_23 = 0.03215047;
  double eta_24 = 0.03584044;
  NumericVector inter1 = NumericVector::create(3.0, 4.5);
  NumericVector inter2 = NumericVector::create(4.5, 6.0);
  double part1 = eta_23 * c2(w,a) * 
    trapezoidal_cpp(inter1, n, eta_13, eta_23, w, a);
  double part2 = eta_24 * c2(w,a) *
    trapezoidal_cpp(inter2, n, eta_14, eta_24, w, a);
  double res = (part1 + part2) * 20.0;
  return(res);
}



// [[Rcpp::export]]
double trape(double w, double a, double n){
  double res1 = first_function(w, a, n);
  double res2 = second_function(w, a, n);
  double res3 = third_function(w, a, n);
  double res4 = fourth_function(w, a, n);
  double res = res1 + res2 + res3 + res4;
  return(res);
}



double simpsons_cpp(NumericVector interval, 
                    double n, 
                    double eta1, 
                    double eta2,
                    double w,
                    double a){
  vec x(n+1);
  vec out(n/2);
  double b = interval[0];
  double c = interval[1];
  double h = (c-b)/n;
  for(int i = 0; i < n+1; i++){
    x[i] = b + i*h;
  }
  for(int j = 0; j < n/2; j++){
    out[j] = h/3*f(x[2*j], eta1, eta2, w, a);
    out[j] = out[j] + 4*h/3*f(x[2*j+1], eta1, eta2, w, a);
    out[j] = out[j] + h/3*f(x[2*j+2], eta1, eta2, w, a);
  }
  double res = accu(out);
  return(res);
}

double first_function2(double w, double a, double n){
  double eta_11 = 0.07116282;
  double eta_12 = 0.07766563;
  double eta_21 = 0.02380466;
  double eta_22 = 0.02865413;
  NumericVector inter1 = NumericVector::create(0.0, 1.5);
  NumericVector inter2 = NumericVector::create(1.5,3.0);
  double part1 = eta_11 * c1(w,a) * 
    simpsons_cpp(inter1, n, eta_11, eta_21, w, a);
  double part2 = eta_12 * c1(w,a) *
    simpsons_cpp(inter2, n, eta_12, eta_22, w, a);
  double res = (part1 + part2) * 0.5;
  return(res);
}

double second_function2(double w, double a, double n){
  double eta_11 = 0.07116282;
  double eta_12 = 0.07766563;
  double eta_21 = 0.02380466;
  double eta_22 = 0.02865413;
  NumericVector inter1 = NumericVector::create(0.0, 1.5);
  NumericVector inter2 = NumericVector::create(1.5, 3.0);
  double part1 = eta_21 * c2(w,a) * 
    simpsons_cpp(inter1, n, eta_11, eta_21, w, a);
  double part2 = eta_22 * c2(w,a) *
    simpsons_cpp(inter2, n, eta_12, eta_22, w, a);
  double res = (part1 + part2) * 5.0;
  return(res);
}

double third_function2(double w, double a, double n){
  double eta_13 = 0.1052774;
  double eta_14 = 0.1061366;
  double eta_23 = 0.03215047;
  double eta_24 = 0.03584044;
  NumericVector inter1 = NumericVector::create(3.0, 4.5);
  NumericVector inter2 = NumericVector::create(4.5, 6.0);
  double part1 = eta_13 * c1(w,a) * 
    simpsons_cpp(inter1, n, eta_13, eta_23, w, a);
  double part2 = eta_14 * c1(w,a) *
    simpsons_cpp(inter2, n, eta_14, eta_24, w, a);
  double res = (part1 + part2) * 10.0;
  return(res);
}

double fourth_function2(double w, double a, double n){
  double eta_13 = 0.1052774;
  double eta_14 = 0.1061366;
  double eta_23 = 0.03215047;
  double eta_24 = 0.03584044;
  NumericVector inter1 = NumericVector::create(3.0, 4.5);
  NumericVector inter2 = NumericVector::create(4.5, 6.0);
  double part1 = eta_23 * c2(w,a) * 
    simpsons_cpp(inter1, n, eta_13, eta_23, w, a);
  double part2 = eta_24 * c2(w,a) *
    simpsons_cpp(inter2, n, eta_14, eta_24, w, a);
  double res = (part1 + part2) * 20.0;
  return(res);
}

// [[Rcpp::export]]
double simpson(double w, double a, double n){
  double res1 = first_function2(w, a, n);
  double res2 = second_function2(w, a, n);
  double res3 = third_function2(w, a, n);
  double res4 = fourth_function2(w, a, n);
  double res = res1 + res2 + res3 + res4;
  return(res);
}


