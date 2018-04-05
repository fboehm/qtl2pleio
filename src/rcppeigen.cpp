// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

//find the inverse of a matrix with RcppEigen
//[[Rcpp::export]]
Eigen::MatrixXd rcppeigen_invert_matrix(const Eigen::Map<Eigen::MatrixXd> & A){
  Eigen::MatrixXd Ainv = A.inverse();
  return Ainv;
}
//find the log-determinant of the diagonal of a matrix with RcppEigen -- useful for dmvnorm
//[[Rcpp::export]]
Eigen::MatrixXd rcppeigen_get_diag(const Eigen::Map<Eigen::MatrixXd> & A){
  Eigen::MatrixXd temp = A.llt().matrixL();
  return temp.diagonal();
}
//find the determinant of a matrix with RcppEigen
//[[Rcpp::export]]
double rcppeigen_get_det(const Eigen::Map<Eigen::MatrixXd> & A){
 return A.determinant();
}
//find the Cholesky decomposition of a matrix with RcppEigen
//[[Rcpp::export]]
  Eigen::MatrixXd rcppeigen_get_chol(const Eigen::Map<Eigen::MatrixXd> & A){
  Eigen::MatrixXd chol = A.llt().matrixL();
  return chol;
}
//find the Cholesky decomposition of a matrix with RcppEigen using a more stable ldlt() function
//[[Rcpp::export]]
  Eigen::MatrixXd rcppeigen_get_chol_stable(const Eigen::Map<Eigen::MatrixXd> & A){
  Eigen::MatrixXd chol = A.ldlt().matrixL();
  return chol;
}
//find the Cholesky decomposition of a matrix with RcppEigen using a more stable ldlt() function
//[[Rcpp::export]]
  Eigen::VectorXd rcppeigen_get_chol_diag(const Eigen::Map<Eigen::MatrixXd> & A){
  Eigen::VectorXd diag = A.ldlt().vectorD();
  return diag;
}
