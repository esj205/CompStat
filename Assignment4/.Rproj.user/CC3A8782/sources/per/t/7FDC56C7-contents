// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;
using namespace std;
using namespace Rcpp;



double cauchy_pdf(double x, double a, double b){
  double y = (x - a) / b;
  double res = 1.0 / (datum::pi * b * (1.0 + y * y));
  return(res);
}


vec cauchy_pdf_vec(vec x, double a, double b){
  vec res;
  vec y;
  vec a_vec(x.n_rows);
  vec b_vec(x.n_rows);
  vec pi_vec(x.n_rows);
  vec ones(x.n_rows, fill::ones);
  a_vec.fill(a);
  b_vec.fill(b);
  pi_vec.fill(datum::pi);
  y = ( x - a_vec ) / b_vec;
  res = ones / ( pi_vec % b_vec % ( ones + y % y ) );
  return res;
}


double f(double x, double mean, double sd){
  double res = normpdf(x, mean, sd) * cauchy_pdf(x, 5.0, 2.0);
  return res;
}

double f2(double x, double mean, double sd){
  double res = 7.84654*normpdf(x, mean, sd) * cauchy_pdf(x, 5.0, 2.0);
  return res;
}

vec f_vec(vec x, double mean, double sd){
  vec res = normpdf(x, mean, sd) % cauchy_pdf_vec(x, 5.0, 2.0);
  return res;
}

vec f2_vec(vec x, double mean, double sd){
  vec constant(x.n_rows);
  constant.fill(7.84654);
  vec res = constant % normpdf(x, mean, sd) % cauchy_pdf_vec(x, 5.0, 2.0);
  return res;
}

double g(double x, double mean, double sd){
  double res = f2(log(x/(1-x)), mean, sd) / ( x*(1-x) );
  return res;
}

double g2(double x, double mean, double sd){
  double res = f2(1/x, mean, sd)/pow(x,2);
  return res;
}



// [[Rcpp::export]]
double riemann_cpp(NumericVector interval, double n, double mean, double sd){
  vec x(n);
  double h = (interval[1]-interval[0])/n;
  for(int i = 0; i<n; i++){
    x[i] = interval[0] + i*h;
  }
  double res = h*accu(f_vec(x, mean, sd));
  return(res);
}


// [[Rcpp::export]]
double trapezoidal_cpp(NumericVector interval, double n, double mean, double sd){
  vec x(n-1);
  double a = interval[0];
  double b = interval[1];
  double h = (b-a)/n;
  for(int i = 1; i < n; i++){
    x[i-1] = a + i*h;
  }
  double low = h/2*f(a, mean, sd);
  double high = h/2*f(b, mean, sd);
  double middle = h*accu(f_vec(x, mean, sd));
  double res = low + high + middle;
  return(res);
}


// [[Rcpp::export]]
double simpsons_cpp(NumericVector interval, double n, double mean, double sd){
  vec x(n+1);
  vec out(n/2);
  double a = interval[0];
  double b = interval[1];
  double h = (b-a)/n;
  for(int i = 0; i < n+1; i++){
    x[i] = a + i*h;
  }
  for(int j = 0; j < n/2; j++){
    out[j] = h/3*f(x[2*j], mean, sd);
    out[j] = out[j] + 4*h/3*f(x[2*j+1], mean, sd);
    out[j] = out[j] + h/3*f(x[2*j+2], mean, sd);
  }
  double res = accu(out);
  return(res);
}








double riemann_cpp2(NumericVector interval, double n, double mean, double sd){
  vec x(n);
  double h = (interval[1]-interval[0])/n;
  for(int i = 0; i<n; i++){
    x[i] = interval[0] + i*h;
  }
  double res = h*accu(f2_vec(x, mean, sd));
  return(res);
}


double trapezoidal_cpp2(NumericVector interval, double n, double mean, double sd){
  vec x(n-1);
  double a = interval[0];
  double b = interval[1];
  double h = (b-a)/n;
  for(int i = 1; i < n; i++){
    x[i-1] = a + i*h;
  }
  double low = h/2*f2(a, mean, sd);
  double high = h/2*f2(b, mean, sd);
  double middle = h*accu(f2_vec(x, mean, sd));
  double res = low + high + middle;
  return(res);
}

double simpsons_cpp2(NumericVector interval, double n, double mean, double sd){
  vec x(n+1);
  vec out(n/2);
  double a = interval[0];
  double b = interval[1];
  double h = (b-a)/n;
  for(int i = 0; i < n+1; i++){
    x[i] = a + i*h;
  }
  for(int j = 0; j < n/2; j++){
    out[j] = h/3*f2(x[2*j], mean, sd);
    out[j] = out[j] + 4*h/3*f2(x[2*j+1], mean, sd);
    out[j] = out[j] + h/3*f2(x[2*j+2], mean, sd);
  }
  double res = accu(out);
  return(res);
}




// [[Rcpp::export]]
double trapezoidal_cpp3(NumericVector interval, double n, double mean, double sd){
  vec x(n-1);
  double a = interval[0];
  double b = interval[1];
  double h = (b-a)/n;
  double middle = 0;
  for(int i = 1; i < n; i++){
    x[i-1] = a + i*h;
  }
  double low = h/2*g(a, mean, sd);
  for(int j = 0; j < n-1; j++){
    middle = middle + g(x[j], mean, sd);
  }
  middle = h*middle;
  double res = low + middle;
  return(res);
}



// [[Rcpp::export]]
double trapezoidal_cpp4(NumericVector interval, double n, double mean, double sd){
  vec x(n-1);
  double a = interval[0];
  double b = interval[1];
  double h = (b-a)/n;
  double middle = 0;
  for(int i = 1; i < n; i++){
    x[i-1] = a + i*h;
  }
  double high = h/2*g2(b, mean, sd);
  for(int j = 0; j < n-1; j++){
    middle = middle + g2(x[j], mean, sd);
  }
  middle = h*middle;
  double res = high + middle;
  return(res);
}




// [[Rcpp::export]]
DataFrame riemann_df(NumericVector interval, 
                  double mean, 
                  double sd,
                  double tol){
  vector<double> estimation;
  vector<double> rel_err;
  vector<double> sub_interval;
  int iter = 1;
  int n = 2;
  double diff = 1.0;
  double old_est;
  double new_est;
  while(diff > tol){
    if(iter==1){
      double first_est = riemann_cpp2(interval, n, mean, sd);
      estimation.push_back(first_est);
      sub_interval.push_back(n);
      rel_err.push_back(0);
      n = n+2;
      iter = iter+1;
      old_est = first_est;
    }
    else{
      new_est = riemann_cpp2(interval, n, mean, sd);
      diff = abs(new_est - old_est)/old_est;
      estimation.push_back(new_est);
      sub_interval.push_back(n);
      rel_err.push_back((new_est-old_est)/old_est);
      n = n+2;
      iter = iter +1;
      old_est = new_est;
    }
  }
  DataFrame res = DataFrame::create(Named("est") = estimation,
                                    Named("sub_int") = sub_interval,
                                    Named("rel_err") = rel_err);
  return(res);
}

// [[Rcpp::export]]
DataFrame trapezoidal_df(NumericVector interval, 
                     double mean, 
                     double sd,
                     double tol){
  vector<double> estimation;
  vector<double> rel_err;
  vector<double> sub_interval;
  int iter = 1;
  int n = 2;
  double diff = 1.0;
  double old_est;
  double new_est;
  while(diff > tol){
    if(iter==1){
      double first_est = trapezoidal_cpp2(interval, n, mean, sd);
      estimation.push_back(first_est);
      sub_interval.push_back(n);
      rel_err.push_back(0);
      n = n+2;
      iter = iter+1;
      old_est = first_est;
    }
    else{
      new_est = trapezoidal_cpp2(interval, n, mean, sd);
      diff = abs(new_est - old_est)/old_est;
      estimation.push_back(new_est);
      sub_interval.push_back(n);
      rel_err.push_back((new_est-old_est)/old_est);
      n = n+2;
      iter = iter +1;
      old_est = new_est;
    }
  }
  DataFrame res = DataFrame::create(Named("est") = estimation,
                                    Named("sub_int") = sub_interval,
                                    Named("rel_err") = rel_err);
  return(res);
}

// [[Rcpp::export]]
DataFrame simpsons_df(NumericVector interval, 
                     double mean, 
                     double sd,
                     double tol){
  vector<double> estimation;
  vector<double> rel_err;
  vector<double> sub_interval;
  int iter = 1;
  int n = 2;
  double diff = 1.0;
  double old_est;
  double new_est;
  while(diff > tol){
    if(iter==1){
      double first_est = simpsons_cpp2(interval, n, mean, sd);
      estimation.push_back(first_est);
      sub_interval.push_back(n);
      rel_err.push_back(0);
      n = n+2;
      iter = iter+1;
      old_est = first_est;
    }
    else{
      new_est = simpsons_cpp2(interval, n, mean, sd);
      diff = abs(new_est - old_est)/old_est;
      estimation.push_back(new_est);
      sub_interval.push_back(n);
      rel_err.push_back((new_est-old_est)/old_est);
      n = n+2;
      iter = iter +1;
      old_est = new_est;
    }
  }
  DataFrame res = DataFrame::create(Named("est") = estimation,
                                    Named("sub_int") = sub_interval,
                                    Named("rel_err") = rel_err);
  return(res);
}








double f3(double x){
  return 1/x;
}

double f3_sum(double a, int n){
  double s = 0;
  double h = (a-1)/n;
  for(int i = 1; i < n+1; i++){
    s = s + f3(1 + (i-0.5)*h);
  }
  return s;
}




// [[Rcpp::export]]
mat triangle(double a, int m){
  mat T(m+1, m+1);
  
  T(0,0) = 0.5*(f3(1)+f3(a));
  for(int i = 1; i < m+1; i++){
    T(i,0) = 0.5 * T(i-1,0) + f3_sum(a, pow(2, i-1))*(a-1)/pow(2.0,i);
    for(int j = 1; j < i+1; j++){
      T(i,j) = ( pow(4,j)*T(i,j-1) - T(i-1,j-1) ) / ( pow(4,j) -1 ); 
    }
  }
  return T;
}

