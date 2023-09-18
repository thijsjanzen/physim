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

#pragma once

#include <memory>
#include <vector>
#include <array>
#include <random>

enum event_type {shift, speciation, extinction, max_num};

enum finish_type {done, extinct, overshoot, conditioning, not_run_yet,
                  max_types};

using num_mat = std::vector< std::array<float, 4 >>;

struct ltab_species {
  enum info_index {time, p_id, self_id, extinct_time};

  ltab_species(float brts, int parent, int ID, float death) {
    data_[time] = brts;
    data_[p_id] = static_cast<float>(parent);
    data_[self_id] = static_cast<float>(ID);
    data_[extinct_time] = death;
  }

  float get_id() const {
    return(data_[self_id]);
  }

  float get_parent() const  {
    return(data_[p_id]);
  }

  void set_death(float d) {
    data_[extinct_time] = d;
  }

  bool is_dead() {
    if (data_[extinct_time] > -1) return true;
    return false;
  }

  ltab_species() {
    data_[time] = -1e6;
    data_[p_id] = -1e6;
    data_[self_id] = -1e6;
    data_[extinct_time] = -1e6;
  }

  std::array<float, 4>& get_data() {
    return data_;
  }

private:
  std::array<float, 4> data_;
};

struct ltable {
  std::vector< ltab_species > data_;

  ltable() {
    data_.emplace_back(ltab_species(0.0,  0, -1, -1));
    data_.emplace_back(ltab_species(0.0, -1,  2, -1));
  }

  void clear() {
    data_.clear();
  }
};

struct species_info {
  species_info(const std::array<float, 2>& m,
               const std::array<float, 2>& l,
               const std::array<float, 2>& s) :
    trait_mu(m), trait_lambda(l), trait_qs(s) {
    max_mu = m[0] > m[1] ? m[0] : m[1];
    max_la = l[0] > l[1] ? l[0] : l[1];
    max_qs = s[0] > s[1] ? s[0] : s[1];
  }

  float mu(size_t trait) const {
    return trait_mu[trait];
  }

  float lambda(size_t trait) const {
    return trait_lambda[trait];
  }

  float shift(size_t trait) const {
    return trait_qs[trait];
  }

  float max_ext()const  {
    return max_mu;
  }

  float max_spec() const {
    return max_la;
  }

  float max_shift() const {
    return max_qs;
  }

private:
  const std::array<float, 2> trait_mu;
  const std::array<float, 2> trait_lambda;
  const std::array<float, 2> trait_qs;
  float max_mu = 0.0;
  float max_la = 0.0;
  float max_qs = 0.0;
};

struct species {
private:
  size_t trait_;

public:
  int id_;
  float mu_;
  float lambda_;
  float shiftprob_;
  float sum_rate_;

  species(size_t trait, int ID, const species_info& info) :
    trait_(trait),  id_(ID), mu_(info.mu(trait)), lambda_(info.lambda(trait)),
    shiftprob_(info.shift(trait)) {
    sum_rate_ = mu_ + lambda_ + shiftprob_;
  }

  void change_trait(const species_info& info) {
    trait_ = 1 - trait_;
    mu_ = info.mu(trait_);
    lambda_ = info.lambda(trait_);
    shiftprob_ = info.shift(trait_);
    sum_rate_ = mu_ + lambda_ + shiftprob_;
  }

  size_t get_trait() const {
    return trait_;
  }
};


struct population {
  std::vector<species> pop;
  std::array<float, 3> rates;

  population() {
    rates = {0.0, 0.0, 0.0};
  }

  void add(const species& s) {
    rates[shift]      += s.shiftprob_;
    rates[extinction] += s.mu_;
    rates[speciation] += s.lambda_;
    pop.push_back(s);
  }

  void remove(size_t index) {
    rates[shift]      -= pop[index].shiftprob_;
    rates[extinction] -= pop[index].mu_;
    rates[speciation] -= pop[index].lambda_;

    pop[index] = pop.back();
    pop.pop_back();
  }

  void change_trait(size_t index, const species_info& info) {
    auto old_mu    = pop[index].mu_;
    auto old_la    = pop[index].lambda_;
    auto old_shift = pop[index].shiftprob_;
    pop[index].change_trait(info);

    rates[shift]      += pop[index].shiftprob_ - old_shift;
    rates[extinction] += pop[index].mu_ - old_mu;
    rates[speciation] += pop[index].lambda_ - old_la;
  }

  bool empty() const {
    return pop.empty();
  }

  size_t size() const {
    return pop.size();
  }

  size_t get_trait(size_t index) const {
    return pop[index].get_trait();
  }

  int get_id(size_t index) const {
    return pop[index].id_;
  }

  void clear() {
    pop.clear();
    rates = {0.0, 0.0, 0.0};
  }
};

struct bisse_sim {
  std::mt19937_64 rndgen_;

  ltable L;
  float t;

  finish_type run_info;

  int init_state;

  population pop;

  std::array<int, 2> track_crowns;

  // external data:
  const species_info trait_info;
  const size_t num_states;
  const float max_t;
  const size_t max_spec;
  const int init_states;

  bisse_sim(const std::array<float, 2>& m,
            const std::array<float, 2>& l,
            const std::array<float, 2>& q,
            float mt,
            size_t max_s,
            int init) :
    trait_info(m, l, q),
    num_states(m.size()), max_t(mt),
    max_spec(max_s),
    init_states(init) {
    // randomize randomizer
    std::random_device rd;
    std::mt19937_64 rndgen_t(rd());
    rndgen_ = rndgen_t;
    run_info = not_run_yet;
    t = 0.0;
    init_state = 0;
  }

  void run() {
    t = 0.0;

    init_state = init_states;
    if (init_state < 0) { // randomly draw initial trait
      std::uniform_int_distribution<size_t> d(0, 2);
      init_state = d(rndgen_);
    }

    run_info = not_run_yet;

    pop.clear();
    L.clear();

    pop.add(species(init_state, -1, trait_info));
    pop.add(species(init_state,  2, trait_info));

    track_crowns = {1, 1};

    L = ltable();

    while (true) {
      float dt = draw_dt();
      t += dt;

      if (t > max_t)  {
        run_info = done; break;
      }

      event_type event = draw_event();
      apply_event(event);

      if (track_crowns[0] < 1 || track_crowns[1] < 1) {
        run_info = extinct;
        break;
      }
      if (pop.size() > max_spec) {
        run_info = overshoot; break;
      }
    }
  }

  void apply_event(const event_type event) {
    switch (event) {
      case shift: {
        event_traitshift();
        break;
      }
      case speciation: {
        event_speciation();
        break;
      }
      case extinction: {
        event_extinction();
        break;
      }
      default: break;
    }
    return;
  }

  void event_extinction() {
    size_t dying = 0;
    if (pop.size() > 1) {
      // sample one at randomly following mus
      auto get_val = [](const species& s) { return s.mu_;};
      dying = sample_from_pop(get_val,
                              trait_info.max_ext());
    }
    auto dying_id = pop.get_id(dying);

    for (auto& i : L.data_) {
      if (std::abs(i.get_id()) == std::abs(dying_id)) {
        if (i.get_id() < 0) {
          track_crowns[0]--;
        } else {
          track_crowns[1]--;
        }

        i.set_death(t);
        break;
      }
    }
    pop.remove(dying);
  }

  void event_speciation() {
    size_t mother = 0;
    if (pop.size() > 1) {
      // sample one at randomly following lambdas
      auto get_val = [](const species& s) { return s.lambda_;};
      mother = sample_from_pop(get_val,
                               trait_info.max_spec());
    }
    auto mother_trait = pop.get_trait(mother);

    int new_id = static_cast<int>(L.data_.size()) + 1;
    if (pop.get_id(mother) < 0) {
      track_crowns[0]++;
      new_id *= -1;
    } else {
      track_crowns[1]++;
    }

    pop.add(species(mother_trait, new_id, trait_info));
    L.data_.emplace_back(ltab_species(t, pop.get_id(mother), new_id, -1));
  }

  void event_traitshift() {
    size_t index_chosen_species = 0;
    if (pop.size() > 1) {
      // sample one at randomly following shiftprob
      auto get_val = [](const species& s) { return s.shiftprob_;};
      index_chosen_species = sample_from_pop(get_val,
                                             trait_info.max_shift());
    }

    pop.change_trait(index_chosen_species, trait_info);
    return;
  }

  event_type draw_event() {
    float total_rate = pop.rates[shift] +
      pop.rates[extinction] +
      pop.rates[speciation];
    std::uniform_real_distribution<float> unif_dist(0.0, total_rate);
    float r = unif_dist(rndgen_);

    // ordering of rates is:
    // {shift, speciation, extinction, max_num};
    if (r < pop.rates[shift]) return shift;
    if (r < pop.rates[shift] + pop.rates[speciation]) return speciation;

    return extinction;
  }

  float draw_dt() {
    float total_rate = pop.rates[shift] +
      pop.rates[extinction] +
      pop.rates[speciation];

    std::exponential_distribution<float> exp_dist(total_rate);
    return exp_dist(rndgen_);
  }

  size_t sample_from_pop(float (*getvalfrom_species)(const species&),
                         float max_val) {

    std::uniform_int_distribution<> d(0, static_cast<int>(pop.size()) - 1);
    std::uniform_real_distribution<float> r(0.0, 1.0);
    int index;
    float mult = 1.0 / max_val;
    float ulim = 1.0 - 1e-9;
    while (true) {
      index = d(rndgen_);
      float rel_prob = getvalfrom_species(pop.pop[index]) * mult;
      if (rel_prob > 0.0) {
        if (rel_prob >= (ulim)) break;

        if (r(rndgen_) < rel_prob) {
          break;
        }
      }
    }
    return index;
  }

  std::vector<int> get_traits() {
    std::vector<int> traits(pop.size() * 2);
    for (size_t i = 0; i < pop.size(); ++i) {
      auto index = i * 2;
      traits[index] = pop.get_trait(i);
      traits[index + 1] = pop.pop[i].id_;
    }
    return traits;
  }

  size_t get_initial_state() {
    return init_state;
  }

  size_t get_num_species() {
    return pop.size();
  }

  num_mat extract_ltable() {
    num_mat extracted_ltable(L.data_.size(), std::array<float, 4>());
    for (int i = 0; i < L.data_.size(); ++i) {
      extracted_ltable[i] = L.data_[i].get_data();
    }
    return extracted_ltable;
  }
};

