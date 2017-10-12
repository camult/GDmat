#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Calculates the realized dominance relationship matrix
//
// [[Rcpp::export]]
arma::mat DRM(arma::mat X){
  arma::vec p = arma::sum(X,1) / (2*X.n_cols);
  arma::vec q = 1 - p;
  arma::vec pq2 = (2*p%q);
  double varH = accu(pq2%(1-pq2));
  arma::mat H = 1 - abs(X-1);
  H.each_col() -= pq2;
  arma::mat D = (H.t()*H)/varH;
  return D;
}


// Calculates the realized additive relationship matrix
//
// [[Rcpp::export]]
arma::mat GRM(arma::mat X, double smallValue){
  arma::vec p = arma::sum(X,1) / (2*X.n_cols);
  arma::vec q = 1 - p;
  double sum2pq = accu(2*p%q);
  arma::vec P = 2*p;
  X.each_col() -= P;
  arma::mat G = (X.t()*X)/sum2pq;
  G.diag() += smallValue;
  return G;
}

