#include <Rcpp.h>
#include "abc.h"
#include "util.h"

void update_output(std::vector< std::array<double, 10>>& out,
                   const std::vector< particle >& gen,
                   int iteration) {

  for (const auto& i : gen) {

    std::array<double, 10> to_add;
    to_add[0] = static_cast<double>(iteration);
    for (size_t j = 0; j < i.params_.size(); ++j) {
      to_add[j + 1] = i.params_[j];
    }
    to_add[6] = i.gamma;
    to_add[7] = i.colless;
    to_add[8] = static_cast<double>(i.num_lin);
    to_add[9] = i.weight;

    out.push_back(to_add);
  }
  return;
}


//' function to do abc using rcpp
//' @param num_particles number of particles
//' @param num_iterations number of iterations
//' @param crown_age crown age
//' @param min_lin minimum number of lineages from the prior
//' @param max_lin maximum number oflineages from the prior
//' @param lambdas vector of lambdas for exponential priors (5)
//' @param s perturbation standard deviation. Perturbations are made on a log scale,
//' e.g. new_param = exp( log(s) + N(0, s))
//' @param obs_gamma observed gamma value to fit on
//' @param obs_colless observed colless value to fit on
//' @param obs_num_lin observed number of lineages to fit on
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix perform_abc_rcpp(int num_particles,
                            int num_iterations,
                            double crown_age,
                            double min_lin,
                            double max_lin,
                            std::vector<double> lambdas,
                            double s,
                            double obs_gamma,
                            double obs_colless,
                            double obs_num_lin) {

    analysis focal_analysis(num_particles,
                            num_iterations,
                            crown_age,
                            min_lin,
                            max_lin,
                            lambdas,
                            s,
                            obs_gamma,
                            obs_colless,
                            obs_num_lin);

    std::vector< std::array<double, 10>> res;

    focal_analysis.iterate_first();
    update_output(res, focal_analysis.current_sample, 0);

    for (size_t i = 1; i < num_iterations; ++i) {
      focal_analysis.iterate(i);
      update_output(res, focal_analysis.current_sample, i);
    }

    Rcpp::NumericMatrix out;
    particle_to_numericmatrix(res, out);

    return out;
}

