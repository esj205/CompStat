// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;
using namespace std;
using namespace Rcpp;


// [[Rcpp::export]]
mat zero_pad(mat dat){
  mat zero_mat = dat;
  zero_mat.insert_rows(0,1);
  zero_mat.insert_rows(dat.n_rows+1, 1);
  zero_mat.inset_cols(0,1);
  zero_mat.inset_cols(dat.n_cols+1, 1);
  return zero_mat;
}
