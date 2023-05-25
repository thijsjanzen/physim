#pragma once

#include <vector>
#include <array>
#include <random>
#include <cmath>

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

  std::vector< std::vector< int > > track_parents;


  const float lambda_g_;
  const float mu_g_;
  const float lambda_i_;
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

  ~pbd_sim() = default;

  void run() {
   // std::cerr << "start\n";
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
      if (track_crowns[0] + track_crowns[1] > max_spec_ && max_spec_ > 0) {
        run_info = overshoot;
        break;
      }
  //    std::cerr << pop_g.size() << " " << pop_i.size() << "\n";
    }

  //  std::cerr << "reducing incipient\n";
    purge_incipient_random();
 //   std::cerr << "simulation done\n";
  }

  void update_rates() {
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
      do_speciation(pop_g, pop_i);
    }
    if (event == speciation_i) {
      do_speciation(pop_i, pop_i);
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

  void do_speciation(const std::vector<float>& parent_pop,
                     std::vector<float>& daughter_pop) {
    auto dist = std::uniform_int_distribution<size_t>(0, parent_pop.size() - 1);

    auto mother = dist(rndgen_);
    auto mother_id = parent_pop[mother];
    float daughter_id = static_cast<float>(L.size() + 1);
    if (mother_id < 0) {
      daughter_id *= -1.f;
      track_crowns[0]++;
    } else {
      track_crowns[1]++;
    }
    L.push_back({t, mother_id, daughter_id, -1.f});
    daughter_pop.push_back(daughter_id);
  }

  void do_extinction(std::vector<float>& focal_pop) {
    auto dist = std::uniform_int_distribution<size_t>(0, focal_pop.size() - 1);
    auto unlucky_spec = dist(rndgen_);
    auto dead_id = focal_pop[unlucky_spec];
    int index = static_cast<int>(dead_id);
    if (index < 0) {
      track_crowns[0]--;
      index *= -1;
    } else {
      track_crowns[1]--;
    }
    index = index - 1; // starting indexing at 1
    L[index][3] = t;
    focal_pop[unlucky_spec] = focal_pop.back();
    focal_pop.pop_back();
  }

  void do_completion() {
    auto dist = std::uniform_int_distribution<size_t>(0, pop_i.size() - 1);
    auto index = dist(rndgen_);
    pop_g.push_back(pop_i[index]);
    pop_i[index] = pop_i.back();
    pop_i.pop_back();
  }

  void purge_incipient_random2() {
    // we want to select one incipient species per good species
    // and make all others extinct.
    std::vector< std::vector< float > > local_parents(L.size() + 10, std::vector<float>());

    for (const auto& i : pop_g) { // these are the good species themselves
      int parent_id = static_cast<int>(std::fabs(i));
      local_parents[parent_id].push_back(parent_id);
    }

    for (const auto& i : pop_i) { // add all incipient daughters
      int incipient_id = static_cast<int>(std::fabs(i));
      int parent_id = static_cast<int>(L[incipient_id - 1][1]);
      parent_id = std::fabs(parent_id);
      local_parents[parent_id].push_back(incipient_id);
    }

    for (size_t j = 0; j < local_parents.size(); ++j) {
      if (local_parents[j].size() > 1) {
        auto dist = std::uniform_int_distribution<size_t>(0, local_parents[j].size() - 1);
        size_t extant_species = dist(rndgen_);
        for (size_t i = 0; i < local_parents[j].size(); ++i) {
          if (i != extant_species) {
            int focal_index = static_cast<int>(std::fabs(local_parents[j][i])) - 1;
            if (focal_index < 0) throw "index < 0";
            L[focal_index][4] = t; // fake them being extinct
          }
        }
      }
    }
    return;
  }

  void purge_incipient_random() {
    track_parents.clear();
    track_parents.resize(L.size() + 1);
    for (size_t i = 0; i < pop_g.size(); ++i) {
      int parent_id = std::abs(static_cast<int>(pop_g[i]));
      track_parents[parent_id].push_back(parent_id);
    }
    for (size_t i = 0; i < pop_i.size(); ++i) {
      int focal_id = std::abs(static_cast<int>(pop_i[i]));
      // find parent:
      int parent_id = std::abs(static_cast<int>(L[focal_id - 1][1]));
      track_parents[parent_id].push_back(focal_id);
    }
    // now we go over all collections and make all extinct except one
    for (size_t i = 0; i < track_parents.size(); ++i) {
      if (track_parents[i].size() > 1) { // if there is only one, we retain that one
        std::uniform_int_distribution<int> d(0, track_parents[i].size() - 1);
        int lucky_one = d(rndgen_);
        for (size_t j = 0; j < track_parents[i].size(); ++j) {
          if (j != lucky_one) {
            L[track_parents[i][j] - 1][3] = t;
          }
        }
      }
    }

    return;
  }




  bool is_present(const std::vector<float>& present, float index) {
    for (const auto& i : present) {
      if (i == index) return true;
    }
    return false;
  }

  int get_num_species() {
    int cnt = 0;
    for (const auto& i : L) {
      if (i[3] < 0) cnt++;
    }
    return cnt;
  }

  size_t get_num_good_species() {
    return pop_g.size();
  }
  size_t get_num_incipient_species() {
    return pop_i.size();
  }
};
