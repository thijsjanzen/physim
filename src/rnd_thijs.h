#pragma once
#include <random>

struct rnd_t {
  std::mt19937_64 rndgen_;

  rnd_t() {
    std::random_device rd;
    rndgen_ = std::mt19937_64(rd());
    unif_dist = std::uniform_real_distribution<>(0, 1.0);
  }

  rnd_t(size_t seed) {
    rndgen_ = std::mt19937_64(seed);
    unif_dist = std::uniform_real_distribution<>(0, 1.0);
  }

  rnd_t(const std::vector<double>& lambda_vals,
        double sigma) :
    lambdas(lambda_vals) {
    std::random_device rd;
    rndgen_ = std::mt19937_64(rd());
    unif_dist = std::uniform_real_distribution<>(0, 1.0);
    perturb = std::normal_distribution<double>(0.0, sigma);
  }

  std::uniform_real_distribution<> unif_dist;

  double uniform() {
    return unif_dist(rndgen_);
  }

  int random_number(unsigned int n) {
    return std::uniform_int_distribution<> (0, n-1)(rndgen_);
  }

  double exp(double lambda) {
    if (lambda == 0.0) return 1e20f;
    return std::exponential_distribution<double>(lambda)(rndgen_);
  }

  double normal(double m, double s) {
    std::normal_distribution<double> d(m, s);
    return(d(rndgen_));
  }

  double perturb_particle_val(double m) {
    return m + perturb(rndgen_);
  }

  std::array<double, 5> draw_from_prior() {
    std::array<double, 5> out;
    for (size_t i = 0; i < out.size(); ++i) {
      out[i] = exp(lambdas[i]);
    }
    return out;
  }

  std::vector<double> lambdas;
  std::normal_distribution<double> perturb;
};
