// Matrix utilities
#ifndef MATRIX_H
#define MATRIX_H

#include <RcppEigen.h>

// find columns that exactly match previous columns
// returns numeric vector with -1 indicating no match to an earlier column and
//                             >0 indicating matches that earlier column
//                                (indexes starting at 1)
Rcpp::NumericVector find_matching_cols(const Rcpp::NumericMatrix& mat,
                                       const double tol);

// find set of linearly independent columns in a matrix
// returns a vector of column indices (starting at 1)
Rcpp::IntegerVector find_lin_indep_cols(const Rcpp::NumericMatrix& mat,
                                        const double tol);

// form X matrix with intcovar
// has_intercept = true indicates that addcovar has an intercept
//                 and so probs reduced by one column
//               = false means probs has full set of columns
//
// This is maybe a bit confusing.
//
// In the has_intercept=true case, an intercept is included in the
// addcovar matrix, and probs is missing the first column.
// The matrix formed is [A P (P.I)] where A=addcovar, P=probs, I=intcovar
//
// In the has_intercept=false case, no intercept is included in the
// addcovar matrix, and the probs have all columns (so each row sums
// to 1). The matrix formed is [P A (P*.I)] where in the P*.I bit we
// drop the first column of the probs when getting interactions with
// intercovar.
Rcpp::NumericMatrix formX_intcovar(const Rcpp::NumericVector& probs,
                                   const Rcpp::NumericMatrix& addcovar,
                                   const Rcpp::NumericMatrix& intcovar,
                                   const int position,
                                   const bool has_intercept);

// multiply each column of a matrix by a set of weights
Rcpp::NumericMatrix weighted_matrix(const Rcpp::NumericMatrix& mat,
                                    const Rcpp::NumericVector& weights);

// multiply each element of a vector by the corresponding weight
Rcpp::NumericVector weighted_3darray(const Rcpp::NumericVector& array,
                                     const Rcpp::NumericVector& weights);

// expand genotype probabilities with intcovar
Rcpp::NumericVector expand_genoprobs_intcovar(const Rcpp::NumericVector& probs, // 3d array ind x prob x pos
                                              const Rcpp::NumericMatrix& intcovar);


// matrix multiplication
Rcpp::NumericMatrix matrix_x_matrix(const Rcpp::NumericMatrix& X,
                                    const Rcpp::NumericMatrix& Y);

// multiply matrix by vector
Rcpp::NumericVector matrix_x_vector(const Rcpp::NumericMatrix& X,
                                    const Rcpp::NumericVector& y);

// multiply matrix by array
Rcpp::NumericVector matrix_x_3darray(const Rcpp::NumericMatrix& X,
                                     Rcpp::NumericVector& A);

#endif // MATRIX_H
