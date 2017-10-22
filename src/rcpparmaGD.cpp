#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]


//' @title Calculates the realized dominance relationship matrix
//' 
//' @description Calculates realized dominance relationship matrix.
//'
//' @param X is a (numeric, \eqn{n}{m}) matrix of genotypes with a number of rows identical to the number of SNPs (m) and a number of columns identical to the number of genotyped individuals (n). Values in the matrix are 0 and 2 for homozygous and 1 for heterozygous.
//' 
//' @examples
//' ## Not to run ##
//' 
//' ## DRM(X)
//'
//' ## End(Not run)
//' 
//' @references VanRaden, Paul M. Efficient methods to compute genomic predictions. Journal of dairy science 91.11 (2008): 4414-4423.
//' 
//' @references Su G, Christensen OF, Ostersen T, Henryon M, Lund MS. 2012. Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293
//' 
//' 
//' @useDynLib GDmat
//'
//' @export DRM
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


//' @title Calculates the realized additive relationship matrix
//' 
//' @description Calculates realized additive (genomic) relationship matrix.
//'
//' @param X is a (numeric, \eqn{n}{m}) matrix of genotypes with a number of rows identical to the number of SNPs (m) and a number of columns identical to the number of genotyped individuals (n). Values in the matrix are 0 and 2 for homozygous and 1 for heterozygous.
//' @param smallValue is a numeric valeu that can be added to the diagonal of the genomic relationship matrix (GRM).
//' 
//' @examples
//' ## Not to run ##
//' 
//' ## GRM(X, smallValue)
//'
//' ## End(Not run)
//' 
//' @references VanRaden, Paul M. Efficient methods to compute genomic predictions. Journal of dairy science 91.11 (2008): 4414-4423.
//' 
//' @references Su G, Christensen OF, Ostersen T, Henryon M, Lund MS. 2012. Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293
//' 
//' 
//' @useDynLib GDmat
//'
//' @export GRM
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

 