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
Rcpp::List gibbs_cpp2(int N, int thin){
  arma::dmat mat(2, N, fill::zeros); //dmat은 double matrix라는 뜻
  double x = 0.0, y = 0.0 ;
  for(int i =0; i < N; i++){
    for(int j=0; j < thin; j++){
      x = R::rgamma(3, 1 / (y*y+4));
      y = R::rnorm(1 / (x+1), 1 / sqrt(2*(x+1)));
    }
    mat(0, i) = x;
    mat(1, i) = y;
  }
  
  Rcpp::List output;
  output["mat"] = mat;
  return(output);
  
}

// [[Rcpp::export]]
Rcpp::List gibbs_cpp3(int N, int thin){
  arma::dvec x_results(N, fill::zeros);
  arma::dvec y_results(N, fill::zeros);
  double x = 0.0, y = 0.0;
  
  for(int i = 0; i<N; i++){
    for(int j = 0; j<thin; j++){
      x = R::rgamma(3, 1 / (y*y+4));
      y = R::rnorm(1 / (x+1), 1 / sqrt(2*(x+1)));
    }
    x_results(i) = x;
    y_results(i) = y;
  }
  
  Rcpp::List output;
  output["x"] = x_results;
  output["y"] = y_results;
  
  return(output);
}
