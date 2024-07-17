#include "sim_pbd.h"

// [[Rcpp::export]]
Rcpp::List sim_pbd_cpp(double la0,
                       double mu0,
                       double la1,
                       double mu1,
                       double trans_rate,
                       double max_t,
                       double max_num_species,
                       int num_tries = 100) {

   sim_pbd sim(la0, la1, mu0, mu1, trans_rate, max_t, max_num_species);

  for (size_t i = 0; i < num_tries; ++i) {
    sim.run();
    if (sim.status == "success") {

      // check if number of alive species is > 1
      size_t num_species = 0;
      for (const auto& i : sim.L) {
        if (i[species_property::death_time] == -1) num_species++;

        if (num_species > 2) break;
      }

      if (num_species > 2) break;
    }
  }

  Rcpp::NumericMatrix ltable_for_r;
  vector_to_numericmatrix(sim.L, ltable_for_r);

  return Rcpp::List::create(Rcpp::Named("ltable") = ltable_for_r,
                            Rcpp::Named("status") = sim.status,
                            Rcpp::Named("Ng") = sim.get_num_good_species(),
                            Rcpp::Named("Ni") = sim.get_num_incipient_species());
}

// [[Rcpp::export]]
Rcpp::List sim_pbd_conditional_cpp(double la0,
                       double mu0,
                       double la1,
                       double mu1,
                       double trans_rate,
                       double max_t,
                       double min_num_species,
                       double max_num_species,
                       int num_tries = 100) {

  sim_pbd sim(la0, la1, mu0, mu1, trans_rate, max_t, 100000);
  size_t num_species = 0;
  bool found_num_species = false;
  for (size_t i = 0; i < num_tries; ++i) {
    sim.run();
    if (sim.status == "success") {

      // check if number of alive species is required number of species
      num_species = 0;
      for (const auto& i : sim.L) {
        if (i[species_property::death_time] == -1) num_species++;
      }

     // std::cerr << num_species << " ";

      if (num_species >= min_num_species &&
          num_species <= max_num_species) {
        found_num_species = true;
        break;
      }
    }
  }

  if (found_num_species) {

  Rcpp::NumericMatrix ltable_for_r;
  vector_to_numericmatrix(sim.L, ltable_for_r);

  return Rcpp::List::create(Rcpp::Named("ltable") = ltable_for_r,
                            Rcpp::Named("status") = sim.status);
  } else {
    return(Rcpp::List::create(Rcpp::Named("status") = "failure"));
  }
}




