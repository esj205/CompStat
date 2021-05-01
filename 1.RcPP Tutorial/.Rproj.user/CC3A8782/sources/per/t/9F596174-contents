#include <Rcpp.h>
using namespace Rcpp;

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
NumericMatrix gibbs_cpp(int N, int thin) {
  NumericMatrix mat(2, N);
  double x = 0.0, y = 0.0;
  
  for(int i =0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rgamma(1, 3, 1 / (y*y +4))[0];
      y = rnorm(1, 1 / (x+1), 1 / sqrt(2*(x+1)))[0];
    }
    mat(0,i) = x;
    mat(1,i) = y;
  }
  return(mat);
}
