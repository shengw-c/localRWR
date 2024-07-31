#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppEigen.h> // For efficient matrix operations
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using Eigen::Map;       // Eigen mapping
using Eigen::MatrixXd;   // Eigen matrix
using Eigen::VectorXd;   // Eigen vector
using Eigen::SelfAdjointEigenSolver;

// [[Rcpp::export]]
NumericVector Cpp_rowSums(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nr);
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      ans[i] += x(i, j);
    }
  }
  return ans;
}

// [[Rcpp::export]]
NumericVector Cpp_colSums(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nc);
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      ans[i] += x(i, j);
    }
  }
  return ans;
}

// [[Rcpp::export]]
arma::mat Arma_Inv(const arma::mat & x) {
  return arma::inv(x); 
}

// [[Rcpp::export]]
List rwr_cpp(NumericMatrix adj_mat, 
             NumericVector init_scores, // init scores
             Nullable<NumericVector> int_scores = R_NilValue, //optional, only necessary when manually run iteration and serve intermediate scores
             double restart_prop = 0.15,
             int num_iter = 1000, 
             double delta = 1e-6
             ) {
  const Map<MatrixXd> AdjMat(as<Map<MatrixXd> >(adj_mat)); 
  const Map<VectorXd> Init_Scores(as<Map<VectorXd> >(init_scores));
  // const Map<VectorXd> Int_Scores(as<Map<VectorXd> >(int_scores));
  
  // Scale scores to sum to 1
  VectorXd Init_Scores_norm;
  if (std::abs(Init_Scores.sum() - 1.0) < 1e-10) { 
    Init_Scores_norm = Init_Scores; // Init_Scores already sum to 1, no need to scale
  } else {
    warning("Normalizing scores to sum as 1...");
    Init_Scores_norm = Init_Scores.array() / Init_Scores.sum(); // Scale scores to sum to 1
  }
  
  // Progress bar setup
  Progress pb(num_iter, true);
  
  MatrixXd W_p = AdjMat * (1.0 - restart_prop); 
  W_p.diagonal().setZero();
  
  VectorXd restart_term = restart_prop * Init_Scores_norm;
  
  VectorXd F;
  if (int_scores.isNotNull()) {
    warning("Using provided int scores to continue propagation...");
    F = as<VectorXd>(int_scores);
  }else{
    F = Init_Scores_norm;
  }
  
  double prevNorm = F.norm();
  bool converged = false;
  
  MatrixXd W_p_t = W_p.transpose();
  for (int i = 0; i < num_iter; ++i) {
    VectorXd newF = (W_p_t * F) + restart_term; 
    double currentNorm = newF.norm(); 
    
    if (std::abs(currentNorm - prevNorm) < delta) {
      converged = true;
      Rcout << "Convergence at iteration " << ++i << std::endl; // message
      break;  // Exit loop early if converged
    }
    
    F = newF;  
    prevNorm = currentNorm; 
    
    pb.increment(); // Update the progress bar
  }
  if (!converged) {
    warning("RWR did not converge within the specified number of iterations.");
  }
  
  Rcpp::List result;
  result["vector"] = wrap(F); // Wrap Eigen vector as R vector
  result["isconverged"] = converged;
  
  return result; 
}
