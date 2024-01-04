#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec getClusters(arma::mat X, arma::vec weights, arma::mat Mean, arma::mat Sigma) {
  
  int n = X.n_rows;
  int p = Sigma.n_cols;
  arma::vec clusterlabels(n);
  int K = weights.n_elem;
  
  arma::mat INV(K*p,p);
  arma::vec NORM(K);
  for (int k = 0; k < K; ++k) {
    arma::mat S = Sigma.submat(k*p,0,k*p + (p-1),p-1);
    float DET = arma::det(S);
    NORM[k] = std::pow(std::sqrt(2 * arma::datum::pi),p) * std::sqrt(DET);
    INV.submat(k*p,0,k*p + (p-1),p-1) = inv(S); 
  }
  
  arma::vec probs(K);
  
  for (int i = 0; i < n; ++i) {
    
    arma::colvec y = trans(X.row(i));
    
    for (int k = 0; k < K; ++k) {
      
      arma::colvec M = arma::trans(Mean.row(k));
      
      arma::mat inverse =  INV.submat(k*p,0,k*p + (p-1),p-1);
      arma::mat MAH = arma::trans(y - M) * inverse * (y - M);
      float normalization = NORM[k];
      float GAUSS = std::exp(-0.5 * MAH[0]) / normalization;
      probs[k] = weights[k] * GAUSS;
    }
    
    clusterlabels[i] = index_max(probs) + 1;
  }
  return clusterlabels;
}
