#pragma once

#include <vector>
#include <array>
#include <random>

enum event_type {speciation_g, extinction_g, speciation_i, extinction_i,
                 completion, max_num};

enum finish_type {done, extinct, overshoot, not_run_yet,
                  max_types};

using num_mat = std::vector< std::vector<float >>;

using ltable = std::vector< std::array<float, 4> >;


struct pbd_sim {
  ltable L;
  float t;
  std::mt19937 rndgen_;
  std::vector< float > pop_g;
  std::vector< float > pop_i;
  std::array<int, 2> track_crowns;
  finish_type run_info;

  const float lambda_g_;
  const float lambda_i_;
  const float mu_g_;
  const float mu_i_;
  const float compl_rate;

  float sum_rate;
  std::array<float, 5> rates;

  const int max_spec_;
  const float max_t_;
  const int seed_;

  pbd_sim(float lambda_0,
         float mu_0,
         float lambda_1,
         float mu_1,
         float trans_rate,
         float max_t,
         int max_spec,
         int s) :
    lambda_g_(lambda_0),
    mu_g_(mu_0),
    lambda_i_(lambda_1),
    mu_i_(mu_1),
    compl_rate(trans_rate),
    max_spec_(max_spec),
    max_t_(max_t),
    seed_(s) {

    init_randomizer();
  }

  void run() {
    t = 0.0;
    pop_g.clear();
    pop_i.clear();
    L.clear();
    run_info = not_run_yet;
    pop_g.push_back(-1.f);
    pop_g.push_back(2.f);
    L.push_back({0.0, 0, -1, -1});
    L.push_back({0.0, -1, 2, -1});

    track_crowns = {1, 1};



    while (true) {
      update_rates();
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
    }
  }

  float update_rates() {
    size_t N_g = pop_g.size();
    size_t N_i = pop_i.size();
    rates[speciation_g] = lambda_g_ * N_g;
    rates[extinction_g] = mu_g_ * N_g;
    rates[speciation_i] = lambda_i_ * N_i;
    rates[extinction_i] = mu_i_ * N_i;
    rates[completion]   = compl_rate * N_i;
    sum_rate = std::accumulate(rates.begin(), rates.end(), 0.0);
  }

  float draw_dt() {
    std::exponential_distribution<float> exp_dist(sum_rate);
    return exp_dist(rndgen_);
  }

  event_type draw_event() {
    auto unif_dist = std::uniform_real_distribution<float>(0.0, sum_rate);
    float r = unif_dist(rndgen_);
    size_t i = 0;
    for (; i < event_type::max_num; ++i) {
      r -= rates[i];
      if (r <= 0) break;
    }
    return static_cast<event_type>(i);
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
    if (event == speciation_g) {
      do_speciation(pop_g);
    }
    if (event == speciation_i) {
      do_speciation(pop_i);
    }
    if (event == extinction_g) {
      do_extinction(pop_g);
    }
    if (event == extinction_i) {
      do_extinction(pop_i);
    }
    if (event == completion) {
      do_completion();
    }
    return;
  }

  void do_speciation(std::vector<float>& focal_pop) {
    auto dist = std::uniform_int_distribution<size_t>(0, focal_pop.size() - 1);
    auto parent_index = dist(rndgen_);

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
