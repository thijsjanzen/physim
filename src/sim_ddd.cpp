#include <Rcpp.h>
#include "sim_dd.h"

void vector_to_numericmatrix(const std::vector< std::array< float, 4 >>& v,
                             Rcpp::NumericMatrix* m) {
  int n_rows = v.size();
  (*m) = Rcpp::NumericMatrix(n_rows, 4);
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < 4; ++j) {
      (*m)(i, j) = v[i][j];
    }
  }
  return;
}

// [[Rcpp::export]]
Rcpp::List sim_ddd_cpp(double la,
                            double mu,
                            double K,
             double max_t,
             double num_species,
             int seed = -1) {

  dd_sim sim(la, mu, K, max_t, num_species, seed);
  while (true) {
    sim.run();
    if (sim.run_info != extinct) {
      auto sim_num_spec = sim.get_num_species();
      if (num_species >= 0) { // conditioning on number of species
          if (sim_num_spec == num_species) {
            break;
          }
      } else { // conditioning on time
          break;
      }
    }
  }
  Rcpp::NumericMatrix ltable_for_r;
  vector_to_numericmatrix(sim.L, &ltable_for_r);
  return Rcpp::List::create(Rcpp::Named("ltable") = ltable_for_r,
                            Rcpp::Named("crown_age") = sim.t);
}
