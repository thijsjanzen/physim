#pragma once

#include <vector>
#include <array>

#include "colless.h"
#include "gamma.h"
#include "rnd_thijs.h"
#include "sim_pbd.h"

using ltable = std::vector< std::array<double, 4>>;

std::vector<double> brts_from_ltable(const ltable& L) {
  std::vector<double> out;
  out.reserve(L.size());
  for (const auto&i : L) {
    out.push_back(i[3]);
  }
  return out;
}

ltable drop_extinct(const ltable& L) {
  ltable out;
  out.reserve(L.size());
  for (const auto& j : L) {
    if (j[3] < 0) {
      out.push_back(j);
    }
  }

  // now we have a broken Ltable, we need to fill in parents of parents
  for (size_t i = 0; i < out.size(); ++i) {
    auto parent = out[i][species_property::parent];
    if (parent != 0) {
      // find parent:
      bool parent_found = false;
      while (!parent_found) {
        for (size_t j = 0; j < i; ++j) {
          if (out[j][species_property::id] == parent) {
            parent_found = true;
            break;
          }
        }
        if (!parent_found) {
          for (size_t j = 0; j < L.size(); ++j) {
            if (L[j][species_property::id] == parent) {
              out[j][species_property::parent] = L[j][species_property::parent];
              out[j][species_property::birth_time] = L[j][species_property::birth_time];
              parent = out[j][species_property::parent];
              break;
            }
          }
        }
      }
    }
  }

  // now, sort by branching time:
  std::sort(out.begin(), out.end(), [&](const auto& a, const auto& b) {
    return a[0] > b[0];
  });

std::cerr << "\n";
  for (const auto& i : L) {
    for (const auto& j : i) {
      std::cerr << j << " ";
    } std::cerr << "\n";
  }

  std::cerr << "\n";
  for (const auto& i : out) {
    for (const auto& j : i) {
      std::cerr << j << " ";
    } std::cerr << "\n";
  }

  // and now we need to renumber to ensure consistent numbering
  for (int i = 0; i < out.size(); ++i) {
    auto id = std::abs(out[i][species_property::id]);
    if (id != (i + 1)) {
      out[i][species_property::id] = out[i][species_property::id] < 0 ? -(i+1) : i+1;
      for (size_t j = i; j < out.size(); ++j) {
        if (out[j][species_property::parent] == id) {
          out[j][species_property::parent] = out[i][species_property::id];
        }
      }
    }
  }

  std::cerr << "\n";
  std::cerr << "\n";
  for (const auto& i : out) {
    for (const auto& j : i) {
      std::cerr << j << " ";
    } std::cerr << "\n";
  }

  return out;
}

struct particle {
  std::array<double, 5> params_;
  double gamma;
  double colless;
  int num_lin;

  double weight;
  double sigma;
  ltable ltable_;


  particle(rnd_t& rndgen,
           double s) : sigma(s) {
    params_ = rndgen.draw_from_prior();
  }

  void perturb(rnd_t& rndgen) {
      size_t index = rndgen.random_number(params_.size());
      double new_val = rndgen.perturb_particle_val(params_[index]);
      params_[index] = exp(new_val);
      return;
  }

  double prob_perturb(const particle& other) {
    static double prefactor = -log(sigma) - 0.5 * log(2 * 3.141592653589793238);
    double s = other.params_.size() * prefactor;

    for (size_t i = 0; i < other.params_.size(); ++i) {
      double d = (params_[i] - other.params_[i]) * 1.0 / sigma;
      s += -0.5 * d * d;
    }
    return exp(s);
  }

  void update_weight(const std::vector<particle>& other) {
    weight = 0.0;
    for (const auto& i : other) {
      weight += prob_perturb(i) * i.weight;
    }
  }

  void sim(double crown_age,
           int min_lin,
           int max_lin) {

    sim_pbd sim(params_, crown_age, max_lin * 10);
    sim.run();
    num_lin = 0;
    for (const auto& i : sim.L) {
      if (i[3] < 0) num_lin++;
    }
    bool crowns_alive = sim.L[0][species_property::death_time] < 0 &&
                        sim.L[1][species_property::death_time] < 0;

    std::cerr << sim.status << " " << num_lin << " ";

    if (sim.status == "success" &&
        num_lin >= min_lin &&
        num_lin <= max_lin &&
        crowns_alive) {

      ltable_ = drop_extinct(sim.L);
      std::cerr << "ltable dropped ";
      // calc gamma and colless


      std::vector<double> brts = brts_from_ltable(ltable_);
      gamma = calc_gamma(brts);
      std::cerr << "gamma calculated\n";

      for (const auto& i : ltable_) {
        for (const auto& j : i) {
          std::cerr << j << " ";
        } std::cerr << "\n";
      }

      colless_stat_ltable s(ltable_);
      colless = static_cast<double>(s.colless());
      std::cerr << "colless calculated ";

    } else {
      colless = 1e6;
      gamma = 1e6;
      num_lin = 1e6;
    }
    return;
  }
};

struct analysis {
  std::vector<particle> current_sample;

  std::vector<particle> new_sample;

  rnd_t rndgen_;

  const double ref_gamma;
  const double ref_colless;
  const double ref_num_lin;
  const double crown_age;
  const double min_lin;
  const double max_lin;
  const int num_particles;
  const double sigma;

  std::vector<double> threshold;

  analysis(int n,
           int num_iterations,
           double ca,
           double minimum_lineages,
           double maximum_lineages,
           std::vector<double> lambdas,
           double s,
           double obs_gamma,
           double obs_colless,
           double obs_num_lin) :
    ref_gamma(obs_gamma),
    ref_colless(obs_colless),
    ref_num_lin(obs_num_lin),
    crown_age(ca),
    min_lin(minimum_lineages),
    max_lin(maximum_lineages),
    num_particles(n),
    sigma(s) {
    rndgen_ = rnd_t(lambdas, sigma);
    for (size_t i = 0; i < num_iterations; ++i) {
      threshold.push_back(1000 * std::exp(-0.5 * (i - 1)));
    }
  }

  void iterate_first() {

    std::cerr << "0--------25--------50--------75--------100\n";
    std::cerr << "*";
    int updateFreq = num_particles / 20;
    if(updateFreq < 1) updateFreq = 1;

    int prev_print = 0;
    std::cerr << "\n";
    while(current_sample.size() < num_particles) {
      auto new_particle = particle(rndgen_, sigma);
      for (auto i : new_particle.params_) {
        std::cerr << i << " ";
      }

      new_particle.sim(crown_age, min_lin, max_lin);

      std::cerr << new_particle.gamma << " " <<
                   new_particle.colless << " " <<
                   new_particle.num_lin << "\n";

      if (new_particle.num_lin >= min_lin &&
          new_particle.num_lin <= max_lin) {
        current_sample.push_back(new_particle);
      }

      if (current_sample.size() % updateFreq == 0 &&
          current_sample.size() != prev_print) {
        std::cerr << "**";
        prev_print = current_sample.size();
      }
    }
  }

  void iterate(int iteration) {
    new_sample.clear();
    std::vector<double> weights (current_sample.size());
    for (size_t i = 0; i < current_sample.size(); ++i) {
      weights[i] = current_sample[i].weight;
    }

    std::discrete_distribution pick_particle(weights.begin(), weights.end());

    std::cerr << "0--------25--------50--------75--------100\n";
    std::cerr << "*";
    int updateFreq = num_particles / 20;
    if(updateFreq < 1) updateFreq = 1;
    size_t prev_print = 0;

    while(new_sample.size() < num_particles) {
        auto new_particle = current_sample[pick_particle(rndgen_.rndgen_)];
        new_particle.perturb(rndgen_);
        new_particle.sim(crown_age, min_lin, max_lin);
        double dist = calc_dist(new_particle);
        if (dist < threshold[iteration]) {
          new_particle.update_weight(current_sample);
          new_sample.push_back(new_particle);
        }

        if (new_sample.size() % updateFreq == 0 &&
            new_sample.size() != prev_print) {
          std::cerr << "**";
          prev_print = new_sample.size();
        }
    }

    current_sample = new_sample;
  }


  double calc_dist(const particle& p) {
    double d1 = (p.gamma - ref_gamma);
    double d2 = (p.colless - ref_colless);
    double d3 = (p.num_lin- ref_num_lin);

    std::array<double, 3> diff;
    diff[0] = (d1 * d1) / ref_gamma;
    diff[1] = (d2 * d2) / ref_colless;
    diff[2] = (d3 * d3) * 1.0 / ref_num_lin;

    return *std::max_element(diff.begin(), diff.end());
  }
};
