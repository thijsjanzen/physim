#include <Rcpp.h>
#include "sim_pbd.h"
#include "util.h"

// [[Rcpp::export]]
Rcpp::List sim_pbd_cpp(double la0,
                       double mu0,
                       double la1,
                       double mu1,
                       double trans_rate,
                       double max_t,
                       double num_species,
                       int seed = -1) {

  pbd_sim sim(la0, mu0, la1, mu1, trans_rate, max_t, num_species, seed);
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
