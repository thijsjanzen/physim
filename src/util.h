#pragma once

#include <Rcpp.h>
#include <vector>
#include <array>

inline void vector_to_numericmatrix(const std::vector< std::array< float, 4 >>& v,
                                    Rcpp::NumericMatrix& m) {
  int n_rows = v.size();
  m = Rcpp::NumericMatrix(n_rows, 4);
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < 4; ++j) {
      m(i, j) = v[i][j];
    }
  }
  return;
}

inline void vector_to_numericmatrix(const std::vector< std::array< double, 4 >>& v,
                                    Rcpp::NumericMatrix& m) {
  int n_rows = v.size();
  m = Rcpp::NumericMatrix(n_rows, 4);
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < 4; ++j) {
      m(i, j) = v[i][j];
    }
  }
  return;
}

inline void particle_to_numericmatrix(const std::vector< std::array<double, 10>>& v,
                                      Rcpp::NumericMatrix& m) {
  int n_rows = v.size();
  m = Rcpp::NumericMatrix(n_rows, 10);
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < 10; ++j) {
      m(i, j) = v[i][j];
    }
  }
  return;
}
