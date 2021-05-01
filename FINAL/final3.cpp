// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;
using namespace std;
using namespace Rcpp;


// [[Rcpp::export]]
mat zero_pad(mat dat){
  mat zero_mat = dat;
  zero_mat.insert_cols(0,1);
  zero_mat.insert_cols(dat.n_cols+1, 1);
  zero_mat.insert_rows(0,1);
  zero_mat.insert_rows(dat.n_rows+1, 1);
  return zero_mat;
}

// [[Rcpp::export]]
double neigh_sum(mat dat){
  mat zero_mat = zero_pad(dat);
  mat neigh_mat = zero_pad(dat);
  for(int i=1; i<(neigh_mat.n_rows-1); i++){
    for(int j=1; j<(neigh_mat.n_cols-1); j++){
      neigh_mat(i,j) = zero_mat(i,j)*(zero_mat(i-1,j) +
        zero_mat(i+1,j) + zero_mat(i,j-1) + zero_mat(i, j+1));
    }
  }
  double res = accu(neigh_mat);
  return(res);
}

// [[Rcpp::export]]
mat boot(double alpha, double beta, mat dat){
  mat temp_mat = dat;
  mat u_mat(1,1);
  double u;
  for(int i=1; i<11; i++){
    for(int j=0; j<dat.n_rows; j++){
      for(int k=0; k<dat.n_cols; k++){
        if(temp_mat(j,k)==0){
          temp_mat(j,k) = 0;
        }
        else{
          temp_mat(j,k) = 1;
          double dens1 = exp(alpha * accu(temp_mat) + 0.5*beta*neigh_sum(temp_mat));
          temp_mat(j,k) = -1;
          double dens2 = exp(alpha * accu(temp_mat) + 0.5*beta*neigh_sum(temp_mat));
          double r = dens1/(dens1+dens2);
          u_mat.randu();
          u = u_mat(0,0);
          if(u < r){
            temp_mat(j,k)=1;
          }
          else{
            temp_mat(j,k)=-1;
          }
        }
      }
    }
  }
  return(temp_mat);
}

// [[Rcpp::export]]
DataFrame bootdf(int iteration, double alpha, double beta, mat dat){
  vector<double> s1;
  vector<double> s2;
  mat temp_mat;
  for(int i = 0; i < iteration; i++){
    temp_mat = boot(alpha, beta, dat);
    s1.push_back(accu(temp_mat));
    s2.push_back(neigh_sum(temp_mat));
  }
  DataFrame res = DataFrame::create(Named("s1") = s1,
                                    Named("s2") = s2);
  return(res);
}
