#include <mutex>
#include <atomic>
#include <tbb/tbb.h>

#include <cmath>
#include <string>

#include <random>

#include <Rcpp.h>

#include "sim_pbd.h"
#include "L2newick.h"
#include "rnd_thijs.h"

size_t get_rcpp_num_threads() {
  auto* nt_env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
  return (nullptr == nt_env)
    ? tbb::task_arena::automatic  // -1
  : static_cast<size_t>(std::atoi(nt_env));
}

std::array<double, 5> draw_from_prior(rnd_t& rndgen,
                                      const std::vector<double>& prior_means) {
  std::array<double, 5> params;
  for (size_t i = 0; i < 5; ++i) {
    params[i] = rndgen.exp(prior_means[i]);
  } // order: la_0, la_1, mu_0, mu_1, comp_time

  return params;
}

//' function to generate 5 parameters from exponential prior, given mean
//' values
//' @param prior_means vector of prior means
//' @return vector of values
//' @export
// [[Rcpp::export]]
std::vector<double> draw_from_prior_rcpp(const std::vector<double>& prior_means) {
  rnd_t rndgen;

  std::vector<double> answ(5);
  for (size_t i = 0; i < 5; ++i) {
    answ[i] = rndgen.exp(prior_means[i]);
  } // order: la_0, la_1, mu_0, mu_1, comp_time

  return answ;
}

//' function to calculate density of parameters
//' @param prior_means vector of prior means
//' @param x vector of x values for which to calculate the density
//' @return vector of values
//' @export
// [[Rcpp::export]]
double prior_dens_rcpp(const std::vector<double>& prior_means,
                  const std::vector<double>& x) {
  double answ = 0.0;
  for (size_t i = 0; i < prior_means.size(); ++i) {
    answ += log(prior_means[i]) - prior_means[i] * x[i];
  }
  return exp(answ);
}

//' simulate many trees drawing from the prior
//' @param num_repl a vector that indicates the time points of water level changes
//' @param crown_age crown age of the tree to be simulated
//' @param min_lin minimum number of lineages
//' @param max_lin maximum number of lineages
//' @return list with two entries: 1) list of newick strings, 2) matrix of
//' used parameter values
//' @export
// [[Rcpp::export]]
 Rcpp::List create_ref_table_tbb_par(int num_repl,
                                     std::vector<double> prior_means,
                                     double crown_age,
                                     int min_lin,
                                     int max_lin) {

   auto num_threads = get_rcpp_num_threads();
   auto global_control = tbb::global_control(tbb::global_control::max_allowed_parallelism, num_threads);


   std::vector< std::string > trees;
   std::vector< std::array<double, 5>> parameter_list;

   auto T0 = std::chrono::high_resolution_clock::now();
   int loop_size = num_repl - trees.size();

   Rcpp::Rcout << "0--------25--------50--------75--------100\n";
   Rcpp::Rcout << "*";

   int updateFreq = num_repl / 20;
   if(updateFreq < 1) updateFreq = 1;

   int prev_update = 0;

   while(trees.size() < num_repl) {
     // yes, this is silly code. Thank you ;)
     for (size_t t = prev_update; t < trees.size(); ++t) {
       if (t % updateFreq == 0) {
         Rcpp::Rcout << "**";
       }
     }
     prev_update = trees.size();




     // loop size can be optimized further, depending on the average success rate
     // e.g. loop_size = loop_size * 1.0f / success_rate
     // this is especially interesting once only a few are left.
     loop_size = num_repl - trees.size();

     std::vector< std::string > add(loop_size);
     std::vector< std::array<double, 5> > add_params(loop_size);
     std::vector< bool > add_flag(loop_size, false);

     tbb::parallel_for(
       tbb::blocked_range<unsigned>(0, loop_size),
       [&](const tbb::blocked_range<unsigned>& r) {

         rnd_t rndgen;

         bool success = false;
         int num_lin = 0;
         for (unsigned i = r.begin(); i < r.end(); ++i) {
           auto parameters = draw_from_prior(rndgen, prior_means);

           auto l_table = sim_once(parameters,
                                   crown_age,
                                   max_lin * 10,
                                   &success,
                                   &num_lin);

           if (success && num_lin >= min_lin && num_lin <= max_lin) {
             auto x = ltable_to_newick(l_table, true);
             add[i]      = x;
             add_params[i] = parameters;
             add_flag[i] = true;
           }
         }
       });

     for(int j = 0; j < add_flag.size(); ++j) {
       if(add_flag[j]) {
         trees.push_back(add[j]);
         parameter_list.push_back(add_params[j]);
       }
     }
   }

   Rcpp::List output(num_repl);
   for(int k = 0; k < trees.size(); ++k) {
     output[k] = trees[k];
   }

   Rcpp::NumericMatrix parameter_matrix(num_repl, 5);
   for(int k = 0; k < parameter_list.size(); ++k) {
     for(int j = 0; j < parameter_list[0].size(); ++j) {
       parameter_matrix(k, j) = parameter_list[k][j];
     }
   }


   auto T1 = std::chrono::high_resolution_clock::now();
   auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::seconds>(T1 - T0).count());
   Rcpp::Rcout << "\ntrees simulated in: " << elapsed << "seconds\n";
   return Rcpp::List::create(Rcpp::Named("trees") = trees,
                             Rcpp::Named("parameters") = parameter_matrix);
 }







//' simulate many trees drawing from the prior serially
//' @param num_repl a vector that indicates the time points of water level changes
//' @param prior_means means of 5 priors
//' @param crown_age crown age of the tree to be simulated
//' @param min_lin minimum number of lineages
//' @param max_lin maximum number of lineages
//' @param num_threads number of threads
//' @return list with two entries: 1) list of newick strings, 2) matrix of
//' used parameter values
//' @export
// [[Rcpp::export]]
Rcpp::List create_ref_table_serial(int num_repl,
                                   std::vector<double> prior_means,
                                    double crown_age,
                                    int min_lin,
                                    int max_lin,
                                    int num_threads) {

  std::vector< std::string > trees;
  std::vector< std::array< double, 5 > > parameter_list;

  auto T0 = std::chrono::high_resolution_clock::now();
  int loop_size = num_repl - trees.size();

  while(trees.size() < num_repl) {
    // loop size can be optimized further, depending on the average success rate
    // e.g. loop_size = loop_size * 1.0f / success_rate
    // this is especially interesting once only a few are left.
    loop_size = num_repl - trees.size();

    std::vector< std::string > add(loop_size);

    std::vector< std::array< double, 5 > > add_params(loop_size);
    std::vector< bool > add_flag(loop_size, false);

    std::random_device rd;
    rnd_t rndgen;
    bool success = false;
    int num_lin = 0;
    for (unsigned i = 0; i < loop_size; ++i) {
      auto parameters = draw_from_prior(rndgen, prior_means);

      auto l_table = sim_once(parameters,
                              crown_age,
                              max_lin * 10,
                              &success,
                              &num_lin);

      if (success && num_lin >= min_lin && num_lin <= max_lin) {
        auto x = ltable_to_newick(l_table, true);
        add[i]      = x;
        add_params[i] = parameters;
        add_flag[i] = true;
      }
    }
    for(int j = 0; j < add_flag.size(); ++j) {
      if(add_flag[j]) {
        trees.push_back(add[j]);
        parameter_list.push_back(add_params[j]);
      }
    }
  }


  Rcpp::List output(num_repl);
  for(int k = 0; k < trees.size(); ++k) {
    output[k] = trees[k];
  }

  Rcpp::NumericMatrix parameter_matrix(num_repl, 5);
  for(int k = 0; k < parameter_list.size(); ++k) {
    for(int j = 0; j < parameter_list[0].size(); ++j) {
      parameter_matrix(k, j) = parameter_list[k][j];
    }
  }


  auto T1 = std::chrono::high_resolution_clock::now();
  auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::seconds>(T1 - T0).count());
  Rcpp::Rcout << "trees simulated in: " << elapsed << "seconds\n";
  return Rcpp::List::create(Rcpp::Named("trees") = trees,
                            Rcpp::Named("parameters") = parameter_matrix);
}
