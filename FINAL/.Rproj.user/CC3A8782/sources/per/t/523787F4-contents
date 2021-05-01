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

double s(double x, double w, double a){
  double res;
  if(x>0 && x<=1.5){
    res = exp(-0.07116282*c1(w,a)*x -0.02380466*c2(w,a)*x);
  }
  else if(x>1.5 && x<=3){
    res = exp(-1.5*c1(w,a)*(0.07116282-0.07766563) - 1.5*c2(w,a)*(0.02380466-0.02865413) -
      0.07766563*c1(w,a)*x - 0.02865413*c2(w,a)*x);
  }
  else if(x>3 && x<=4.5){
    res = exp(-1.5*c1(w,a)*(0.07116282+0.07766563) - 1.5*c2(w,a)*(0.02380466+0.02865413) +
      3*(0.1052774*c1(w,a)+0.03215047*c2(w,a)) - x*(0.1052774*c1(w,a)+0.03215047*c2(w,a)));
  }
  else{
    res = exp(-1.5*c1(w,a)*(0.07116282+0.07766563+0.1052774) - 1.5*c2(w,a)*(0.02380466+0.02865413+0.03215047) +
      4.5*(0.1061366*c1(w,a) + 0.03584044*c2(w,a)) - x*(0.1061366*c1(w,a) + 0.03584044*c2(w,a)));
  }
  return(res);
}

double f1(double x, double t, double w, double a){
  double eta_11 = 0.07116282;
  double eta_21 = 0.02380466;
  double res1 = 0.5 * eta_11 * c1(w,a) * s(t,w,a);
  double res2 = 5.0 * eta_21 * c2(w,a) * s(t,w,a);
  double res = res1+res2;
  return(res);
}

double f2(double x, double t, double w, double a){
  double eta_12 = 0.07766563;
  double eta_22 = 0.02865413;
  double res1 = 0.5 * eta_12 * c1(w,a) * s(t,w,a);
  double res2 = 5.0 * eta_22 * c2(w,a) * s(t,w,a);
  double res = res1+res2;
  return(res);
}

double f3(double x, double t, double w, double a){
  double eta_13 = 0.1052774;
  double eta_23 = 0.03215047;
  double res1 = 10.0 * eta_13 * c1(w,a) * s(t,w,a);
  double res2 = 20.0 * eta_23 * c2(w,a) * s(t,w,a);
  double res = res1+res2;
  return(res);
}

double f4(double x, double t, double w, double a){
  double eta_14 = 0.1061366;
  double eta_24 = 0.03584044;
  double res1 = 10.0 * eta_14 * c1(w,a) * s(t,w,a);
  double res2 = 20.0 * eta_24 * c2(w,a) * s(t,w,a);
  double res = res1+res2;
  return(res);
}

// [[Rcpp::export]]
double trape2(double t, double w, double a, int n){
  double h = 1.5/n;
  double fsum;
  
  fsum = 0;
  for(int i = 1; i < n; i++){
    fsum = fsum + f1(i*h, t, w, a);
  }
  double res1 = h*(fsum + 0.5 * f1(0,t,w,a) + 0.5 * f1(1.5,t,w,a));
  
  fsum = 0;
  for(int i = 1; i < n; i++){
    fsum = fsum + f2((1.5+i*h),t, w, a);
  }
  double res2 = h*(fsum + 0.5 * f2(1.5,t,w,a) + 0.5 * f2(3.0,t,w,a));
  
  fsum = 0;
  for(int i = 1; i < n; i++){
    fsum = fsum + f3((3.0+i*h),t, w, a);
  }
  double res3 = h*(fsum + 0.5 * f3(3.0,t,w,a) + 0.5 * f3(4.5,t,w,a));
  
  fsum = 0;
  for(int i = 1; i < n; i++){
    fsum = fsum + f4((4.5+i*h),t, w, a);
  }
  double res4 = h*(fsum + 0.5 * f4(4.5,t,w,a) + 0.5 * f3(6.0,t,w,a));
  
  double res = res1+res2+res3+res4;
  return(res);
}

// [[Rcpp::export]]
double simpson2(double t, double w, double a, int n){
  double h = 1.5/n;
  double fsum1;
  double fsum2;
  double fsum3;
  
  fsum1 = 0;
  fsum2 = 0;
  fsum3 = 0;
  for(int i = 1; i < (n/2); i++){
    fsum1 = fsum1 + f1((2*i-2)*h, t, w, a);
    fsum2 = fsum2 + 4.0 * f1((2*i-1)*h, t, w, a);
    fsum3 = fsum3 + f1((2*i)*h, t, w, a);
  }
  double res1 = h*(fsum1 + fsum2 + fsum3)/3.0;
  
  fsum1 = 0;
  fsum2 = 0;
  fsum3 = 0;
  for(int i = 1; i < (n/2); i++){
    fsum1 = fsum1 + f2(1.5+(2*i-2)*h, t, w, a);
    fsum2 = fsum2 + 4.0 * f2(1.5+(2*i-1)*h, t, w, a);
    fsum3 = fsum3 + f2(1.5+(2*i)*h, t, w, a);
  }
  double res2 = h*(fsum1 + fsum2 + fsum3)/3.0;
  
  fsum1 = 0;
  fsum2 = 0;
  fsum3 = 0;
  for(int i = 1; i < (n/2); i++){
    fsum1 = fsum1 + f3(3.0+(2*i-2)*h, t, w, a);
    fsum2 = fsum2 + 4.0 * f3(3.0+(2*i-1)*h, t, w, a);
    fsum3 = fsum3 + f3(3.0+(2*i)*h, t, w, a);
  }
  double res3 = h*(fsum1 + fsum2 + fsum3)/3.0;
  
  fsum1 = 0;
  fsum2 = 0;
  fsum3 = 0;
  for(int i = 1; i < (n/2); i++){
    fsum1 = fsum1 + f4(4.5+(2*i-2)*h, t, w, a);
    fsum2 = fsum2 + 4.0 * f4(4.5+(2*i-1)*h, t, w, a);
    fsum3 = fsum3 + f4(4.5+(2*i)*h, t, w, a);
  }
  double res4 = h*(fsum1 + fsum2 + fsum3)/3.0;
  
  double res = res1+res2+res3+res4;
  return(res);
}
