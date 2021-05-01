// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace std;
using namespace Rcpp;

// Multinomial Example
NumericVector multi_e(NumericVector x, NumericVector p){
  double n1;
  double n2;
  double n3;
  n1 = 38;
  n2 = 34;
  n3 = 125*(p[2]/(p[3]+p[2]));
  NumericVector n = {n1,n2,n3};
  return n;
}

NumericVector multi_m(NumericVector n){
  double theta = (34+n[2])/(72+n[2]);
  NumericVector p = {1.0/2.0-theta/2, theta/4, theta/4, 1.0/2.0};
  return p;
}

// [[Rcpp::export]]
NumericVector multi_out(NumericVector x, NumericVector p, int itr){
  
  for(int i = 0; i<itr; i++){
    NumericVector n = multi_e(x, p);
    p = multi_m(n);
  }
  
  return(p);
  
}

// Peppered moths Example
NumericVector pep_e1(NumericVector x, NumericVector p){
  double n_cc = (x[0]*pow(p[0],2))/(pow(p[0],2)+2*p[0]*p[1]+2*p[0]*p[2]);
  double n_ci = (2*x[0]*p[0]*p[1])/(pow(p[0],2)+2*p[0]*p[1]+2*p[0]*p[2]);
  double n_ct = (2*x[0]*p[0]*p[2])/(pow(p[0],2)+2*p[0]*p[1]+2*p[0]*p[2]);
  double n_ii = (x[1]*pow(p[1],2))/(pow(p[1],2)+2*p[1]*p[2]);
  double n_it = (2*x[1]*p[1]*p[2])/(pow(p[1],2)+2*p[1]*p[2]);
  NumericVector n = {n_cc,n_ci,n_ct,n_ii,n_it,x[2]};
  return n;
}

NumericVector pep_m1(NumericVector x, NumericVector n){
  double p_c = (2*n[0]+n[1]+n[2])/(2*sum(x));
  double p_i = (2*n[3]+n[4]+n[1])/(2*sum(x));
  double p_t = (2*n[5]+n[2]+n[4])/(2*sum(x));
  NumericVector p = {p_c,p_i,p_t};
  return p;
}

// [[Rcpp::export]]
NumericVector pep_out1(NumericVector x, NumericVector p, int itr){
  
  for(int i = 0; i<itr; i++){
    NumericVector n = pep_e1(x, p);
    p = pep_m1(x, n);
  }
  
  return(p);
  
}


// Exercise 4.1
NumericVector pep_e2(NumericVector x, NumericVector p){
  double n_cc = (x[0]*pow(p[0],2))/(pow(p[0],2)+2*p[0]*p[1]+2*p[0]*p[2]);
  double n_ci = (2*x[0]*p[0]*p[1])/(pow(p[0],2)+2*p[0]*p[1]+2*p[0]*p[2]);
  double n_ct = (2*x[0]*p[0]*p[2])/(pow(p[0],2)+2*p[0]*p[1]+2*p[0]*p[2]);
  double n_ii = (x[1]*pow(p[1],2))/(pow(p[1],2)+2*p[1]*p[2]) + ((x[3])*pow(p[1],2))/(pow(p[1],2)+2*p[1]*p[2]+pow(p[2],2));
  double n_it = (2*x[1]*p[1]*p[2])/(pow(p[1],2)+2*p[1]*p[2]) + (2*x[3]*p[1]*p[2])/(pow(p[1],2)+2*p[1]*p[2]+pow(p[2],2));
  double n_tt = x[2] + ((x[3])*pow(p[2],2))/(pow(p[1],2)+2*p[1]*p[2]+pow(p[2],2));
  NumericVector n = {n_cc,n_ci,n_ct,n_ii,n_it,n_tt};
  return n;
}

NumericVector pep_m2(NumericVector x, NumericVector n){
  double p_c = (2*n[0]+n[1]+n[2])/(2*sum(x));
  double p_i = (2*n[3]+n[4]+n[1])/(2*sum(x));
  double p_t = (2*n[5]+n[2]+n[4])/(2*sum(x));
  NumericVector p = {p_c,p_i,p_t};
  return p;
}

// [[Rcpp::export]]
NumericVector pep_out2(NumericVector x, NumericVector p, int itr){
  
  for(int i = 0; i<itr; i++){
    NumericVector n = pep_e2(x, p);
    p = pep_m2(x, n);
  }
  
  return(p);
  
}


// Exercise 4.2
// [[Rcpp::export]]
vec hiv_out(vec param, vec obs, int itr){
  vec new_param = zeros(4);
  double len = obs.size();
  double N = sum(obs);
  double z_0, t_sum, p_sum, z_sum;
  vec t(len);
  vec p(len);
  vec pi_f(len);
  
  for(int i = 0; i<itr; i++){
    //E-step
    //statistics update
    //0 group setting
    pi_f[0] = param[0] + param[1]*exp(-param[2]) + (1-param[0]-param[1])*exp(-param[3]);
    z_0 = param[0] / pi_f[0];
    t[0] = (param[1]*exp(-param[2])) / pi_f[0];
    p[0] = ( (1-param[0]-param[1])*exp(-param[3]) ) / pi_f[0];
    
    //1~17 group update
    for(int j = 1; j < len; j++){
      pi_f[j] = param[1]*pow(param[2], j)*exp(-param[2]) + (1-param[0]-param[1])*pow(param[3], j)*exp(-param[3]);
      t[j] = (param[1]*pow(param[2], j)*exp(-param[2])) / pi_f[j];
      p[j] = ((1-param[0]-param[1])*pow(param[3], j)*exp(-param[3])) / pi_f[j];
    }
    
    t_sum = 0;
    p_sum = 0;
    z_sum = 0;
    
    z_sum = obs[0] * z_0;
    
    for(int k=0; k < len; k++){
      t_sum = t_sum + obs[k]*t[k];
      p_sum = p_sum + obs[k]*p[k];
    }
    
    // M-step
    new_param = {0, 0, 0, 0};
    new_param[0] = z_sum / N ;
    for(int l=0; l<len; l++){
      new_param[1] = new_param[1] + (obs[l]*t[l]) / N ;
      new_param[2] = new_param[2] + (l*obs[l]*t[l]) / t_sum ;
      new_param[3] = new_param[3] + (l*obs[l]*p[l]) / p_sum ;
      
    }
    
    param = new_param;
    
  }
  
  return(param);
  
}

