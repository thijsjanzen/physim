// Copyright 2022 - 2023 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
//
#include <Rcpp.h>

#include "bisse_sim.h"   // NOLINT [build/include_subdir]
#include "util.h"

// [[Rcpp::export]]
Rcpp::List bisse_sim_cpp(const std::vector<float>& pars,
                         float crown_age,
                         int num_species,
                         float init_state,
                         bool verbose,
                         int max_tries) {

    bisse_sim sim( {pars[2], pars[3]}, // mu0, mu1
                   {pars[0], pars[1]}, // lambda0, lambda1
                   {pars[4], pars[5]}, // q01, q10
                    crown_age,
                    num_species * 1.05,
                    init_state);
   for (size_t cnt = 0; cnt < max_tries; ++cnt) {
    sim.run();
    if (sim.get_num_species() == num_species) {
      break;
    }
   }

    // extract and return
    Rcpp::NumericMatrix ltable_for_r;
    vector_to_numericmatrix(sim.extract_ltable(), ltable_for_r);

    auto traits = sim.get_traits();
    auto init   = sim.get_initial_state();

    Rcpp::List output = Rcpp::List::create(Rcpp::Named("ltable") = ltable_for_r,
                                           Rcpp::Named("traits") = traits,
                                           Rcpp::Named("initial_state") = init);
    return output;
}
