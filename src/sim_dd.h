#pragma once

#include <vector>
#include <array>
#include <random>

enum event_type {speciation, extinction, max_num};

enum finish_type {done, extinct, overshoot, not_run_yet,
                  max_types};

using num_mat = std::vector< std::vector<float >>;

using ltable = std::vector< std::array<float, 4> >;


struct dd_sim {
  ltable L;
  float t;
  std::mt19937 rndgen_;
  std::vector< float > pop;
  std::array<int, 2> track_crowns;
  finish_type run_info;

  const float lambda_;
  const float mu_;
  const float K_;

  float spec_rate;
  float ext_rate;

  const int max_spec_;
  const float max_t_;
  const int seed_;

  dd_sim(float lambda,
         float mu,
         float K,
         float max_t,
         int max_spec,
         int s) :
    lambda_(lambda),
    mu_(mu),
    K_(K),
    max_spec_(max_spec),
    max_t_(max_t),
    seed_(s) {

    init_randomizer();
  }

  void run() {
    t = 0.0;
    pop.clear();
    L.clear();
    run_info = not_run_yet;
    pop.push_back(-1.f);
    pop.push_back(2.f);
    L.push_back({0.0, 0, -1, -1});
    L.push_back({0.0, -1, 2, -1});

    track_crowns = {1, 1};

    ext_rate = mu_;

    while (true) {
      spec_rate = calc_spec_rate(pop.size());
      ext_rate = mu_;

      float dt = draw_dt();
      t += dt;
      if (t > max_t_ && max_t_ > 0.0)  {
        run_info = done; break;
      }
      event_type event = draw_event();
      apply_event(event);

      if (track_crowns[0] < 1 || track_crowns[1] < 1) {
        run_info = extinct;
        break;
      }
      if (pop.size() > max_spec_ && max_spec_ > 0) {
        run_info = overshoot; break;
      }
    }
  }

  float draw_dt() {
    std::exponential_distribution<float> exp_dist(pop.size() * (spec_rate + ext_rate));
    return exp_dist(rndgen_);
  }

  event_type draw_event() {
    auto unif_dist = std::uniform_real_distribution<float>(0.0, spec_rate + ext_rate);
    float r = unif_dist(rndgen_);
    if (r <= spec_rate) {
      return speciation;
    } else {
      return extinction;
    }
  }

  float calc_spec_rate(size_t N) {
    float answ =  lambda_ - (lambda_ - mu_) * N / K_;
    return answ  < 0 ? 0 : answ;
  }

  void init_randomizer() {
    if (seed_ >= 0) {
      std::mt19937 rndgen_t(seed_);
      rndgen_= rndgen_t;
    } else {
      std::random_device rd;
      std::mt19937 rndgen_t(rd());
      rndgen_ = rndgen_t;
    }
  }

  void apply_event(const event_type event) {
    if (event == speciation) {
      do_speciation();
    }
    if (event == extinction) {
      do_extinction();
    }
  }

  void do_speciation() {
    auto dist = std::uniform_int_distribution<size_t>(0, pop.size() - 1);
    auto mother = dist(rndgen_);
    auto mother_id = pop[mother];
    float daughter_id = static_cast<float>(L.size() + 1);
    if (mother_id < 0) {
      daughter_id *= -1.f;
      track_crowns[0]++;
    } else {
      track_crowns[1]++;
    }
    L.push_back({t, mother_id, daughter_id, -1.f});
    pop.push_back(daughter_id);
  }

  void do_extinction() {
    auto dist = std::uniform_int_distribution<size_t>(0, pop.size() - 1);
    auto unlucky_spec = dist(rndgen_);
    auto dead_id = pop[unlucky_spec];
    int index = static_cast<int>(dead_id);
    if (index < 0) {
      track_crowns[0]--;
      index *= -1;
    } else {
      track_crowns[1]--;
    }
    index = index - 1; // starting indexing at 1
    L[index][3] = t;
    pop[unlucky_spec] = pop.back();
    pop.pop_back();
  }

  int get_num_species() {
    return pop.size();
  }
};
