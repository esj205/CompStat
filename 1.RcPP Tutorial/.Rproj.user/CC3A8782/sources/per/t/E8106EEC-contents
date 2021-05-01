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

// [[Rcpp::export]]
Rcpp::List gibbs_cpp(int N, int thin){
  arma::dmat mat(N, 2, fill::zeros);
  double x = 0.0, y = 0.0;
  //int x = 0, y = 0;
  
  for(int i = 0; i < N; i++){
    for(int j = 0; j < thin; j++){
      x = R::rgamma(3, 1 / (y * y + 4));
      y = R::rnorm(1 / (x + 1), 1 / sqrt(2 * (x + 1)));
    }
    mat(i,0) = x;
    mat(i,1) = y;
  }
  
  Rcpp::List output;
  output["mat"] = mat;
  
  return(output);
}

// [[Rcpp::export]]
Rcpp::List gibbs_cpp2(int N, int thin){
  arma::dvec x_result(N, fill::zeros);
  arma::dvec z_result(N, fill::zeros);
  double x = 0.0, y = 0.0;
  //int x = 0, y = 0;
  
  for(int i = 0; i < N; i++){
    for(int j = 0; j < thin; j++){
      x = R::rgamma(3, 1 / (y * y + 4));
      y = R::rnorm(1 / (x + 1), 1 / sqrt(2 * (x + 1)));
    }
    x_result(i) = x;
    z_result(i) = y;
  }
  
  Rcpp::List output;
  output["x"] = x_result;
  output["y"] = z_result;
  
  return(output);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
