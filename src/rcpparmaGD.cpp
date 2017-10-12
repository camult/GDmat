#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 

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


// Calculates the numerator relationship matrix. Code from optiSel R package
//
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_makeA(const arma::ivec& numSire, const arma::ivec& numDam, const arma::mat& AFounder, const arma::ivec& numFounder, const Rcpp::CharacterVector IndivName){
  int i, j, nSire, nDam, ni, nj;
  int N = numSire.n_elem;
  int NFounder = numFounder.n_elem;
  double tmp;
  Rcpp::NumericMatrix rpedKin(N, N);
  arma::mat pedKin(rpedKin.begin(), rpedKin.nrow(), rpedKin.ncol(), false);
  pedKin.eye();
  
  for(i=0;i<NFounder;i++){
    ni = numFounder.at(i);
    for(j=0;j<=i;j++){
      nj  = numFounder.at(j);
      tmp = AFounder.at(i, j);
      pedKin.at(ni, nj) = tmp;
      pedKin.at(nj, ni) = tmp;
    }
  }
  
  for(i=0;i<N;i++){
    nSire = numSire.at(i);
    nDam  = numDam.at(i);
    if((nSire>0) & (nDam>0)){pedKin.at(i, i) = 1.0 + 0.5*pedKin.at(nSire-1, nDam-1);}
    if((nSire>0) | (nDam>0)){
      for(j=0; j<i; j++){
        tmp = ((nSire>0)?(0.5*pedKin.at(nSire-1, j)):(0.0)) + ((nDam>0)?(0.5*pedKin.at(nDam-1, j)):(0.0));
        pedKin.at(i, j) = tmp;
        pedKin.at(j, i) = tmp;
      }
    }
  }
  
  rpedKin.attr("dimnames") = Rcpp::List::create(IndivName, IndivName);
  return rpedKin;
}


// Calculates the numerator relationship matrix. Code from optiSel R package
//
// [[Rcpp::export]]

Rcpp::NumericMatrix rcpp_makeA_lowMem(const arma::ivec& numSire, const arma::ivec& numDam, const arma::mat& AFounder, const arma::ivec& numFounder, const Rcpp::CharacterVector IndivName,  const arma::ivec& numKeep, const arma::ivec& ainKeep, const arma::ivec& anOff){
  int i, j, nSire, nDam, ni, nj, nextFounder, iF;
  int N = numSire.n_elem;
  int NFounder = numFounder.n_elem;
  int NKeep    = numKeep.n_elem;
  
  double** Amat  = (double**)calloc(N, sizeof(double*));  
  int*  nOff     = (int*)calloc(N, sizeof(int));  
  int*  inKeep   = (int*)calloc(N, sizeof(int));  
  int*  isAlloc  = (int*)calloc(N, sizeof(int));  
  
  iF = 0;
  nextFounder =  numFounder.at(iF);
  
  for(i=0;i<N;i++){
    nOff[i]    = anOff.at(i);
    inKeep[i]  = ainKeep.at(i);
    Amat[i]    = (double*) calloc(i+1, sizeof(double));
    isAlloc[i] = 1;
    if(i==nextFounder){
      for(j=0;j<=iF;j++){
        Amat[i][numFounder.at(j)] = AFounder.at(iF, j);
      }
      if(iF<NFounder-1){iF=iF+1;}
      nextFounder =  numFounder.at(iF);
    }else{
      Amat[i][i] = 1.0;
    }
    nSire = numSire.at(i);
    nDam  = numDam.at(i);
    if((nSire>0) & (nDam>0)){
      Amat[i][i] = 1.0 + 0.5*((nDam<nSire)?(Amat[nSire-1][nDam-1]):(Amat[nDam-1][nSire-1]));
    }
    if(nSire>0){
      nOff[nSire-1] = nOff[nSire-1]-1;
      for(j=0; j<nSire-1; j++){
        Amat[i][j]  = 0.5*Amat[nSire-1][j];
      }
      for(j=nSire-1; j<i; j++){
        if(isAlloc[j]){
          Amat[i][j]  = 0.5*Amat[j][nSire-1];
        }
      }
      if((!nOff[nSire-1]) & (!inKeep[nSire-1])){
        free(Amat[nSire-1]);
        isAlloc[nSire-1] = 0;
      }
    }
    if(nDam>0){
      nOff[nDam-1] = nOff[nDam-1]-1;
      for(j=0; j<nDam-1; j++){
        Amat[i][j]  += 0.5*Amat[nDam-1][j];
      }
      for(j=nDam-1; j<i; j++){
        if(isAlloc[j]){
          Amat[i][j] += 0.5*Amat[j][nDam-1];
        }
      }
      if((!nOff[nDam-1]) & (!inKeep[nDam-1])){
        free(Amat[nDam-1]);
        isAlloc[nDam-1] = 0;
      }
    }
  }
  
  Rcpp::NumericMatrix rAmat(NKeep, NKeep);
  arma::mat aAmat(rAmat.begin(), rAmat.nrow(), rAmat.ncol(), false);
  for(i=0;i<NKeep;i++){
    ni = numKeep.at(i);
    for(j=0;j<=i;j++){
      nj = numKeep.at(j);
      aAmat.at(i, j) = Amat[ni][nj];
      aAmat.at(j, i) = Amat[ni][nj];
    }
  }
  for(i=0;i<N;i++){
    if(isAlloc[i]){
      free(Amat[i]);
    }
  }
  free(Amat);
  free(isAlloc);
  free(inKeep);
  free(nOff);
  
  rAmat.attr("dimnames") = Rcpp::List::create(IndivName, IndivName);
  return rAmat;
}
