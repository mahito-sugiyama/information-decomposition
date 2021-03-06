/*
Information decomposition on structured space
Copyright (C) 2016 Mahito Sugiyama

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

Contact:
Mahito Sugiyama
ISIR, Osaka University,
8-1, Mihogaoka, Ibaraki-shi, Osaka, 567-0047, Japan
E-mail: mahito@ar.sanken.osaka-u.ac.jp
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <list>
#include <numeric>
#include <random>
#include <functional>
#include <utility>
#include <iomanip>

using namespace std;

// node structure
class node {
public:
  bool is_lower;
  bool is_upper;
  int32_t id, id_org;
  int32_t check;
  int32_t supp;
  double p, p_init, p_prev;
  double theta, theta_init, theta_prev;
  double eta, eta_init, eta_prev;
  double score, pvalue;
  vector<reference_wrapper<node>> from;
  vector<reference_wrapper<node>> to;
};

// output a poset "s"
ostream &operator<<(ostream& out, const vector<node>& s) {
  int32_t width = 8;
  out << setw(4) << right << "id" << setw(width) << right << "prob" << setw(width) << right << "theta" << setw(width) << right << "eta" << endl;
  for (auto&& x : s) {
    out << setw(4) << right << setprecision(4) << fixed << x.id << setw(width) << right << x.p << setw(width) << right << x.theta << setw(width) << right << x.eta << endl;
  }
  return out;
}

// for sorting
struct lessId { bool operator()(const node& left, const node& right) const { return left.id < right.id; } };
struct greaterId { bool operator()(const node& left, const node& right) const { return left.id > right.id; } };

// for "reverse" in range-based loop
template<class Cont> class const_reverse_wrapper {
  const Cont& container;
public:
  const_reverse_wrapper(const Cont& cont) : container(cont){ }
  decltype(container.rbegin()) begin() const { return container.rbegin(); }
  decltype(container.rend()) end() const { return container.rend(); }
};
template<class Cont> class reverse_wrapper {
  Cont& container;
public:
  reverse_wrapper(Cont& cont) : container(cont){ }
  decltype(container.rbegin()) begin() { return container.rbegin(); }
  decltype(container.rend()) end() { return container.rend(); }
};
template<class Cont> const_reverse_wrapper<Cont> reverse(const Cont& cont) {
  return const_reverse_wrapper<Cont>(cont);
}
template<class Cont> reverse_wrapper<Cont> reverse(Cont& cont) {
  return reverse_wrapper<Cont>(cont);
}

// get upper and lower elements
void getUpperNext(node& x) {
  if (x.is_upper)
    return;
  x.is_upper = true;
  for (auto&& ptr : x.to)
    getUpperNext(ptr.get());
}
void getUpper(node& x_target, vector<node>& s) {
  for (auto&& x : s)
    x.is_upper = false;
  x_target.is_upper = true;
  for (auto&& ptr : x_target.to)
    getUpperNext(ptr.get());
}
void getLowerNext(node& x) {
  if (x.is_lower)
    return;
  x.is_lower = true;
  for (auto&& ptr : x.from)
    getLowerNext(ptr.get());
}
void getLower(node& x_target, vector<node>& s) {
  for (auto&& x : s)
    x.is_lower = false;
  x_target.is_lower = true;
  for (auto&& ptr : x_target.from)
    getLowerNext(ptr.get());
}

void revertCheckTo(node& x) {
  if (x.check == -1)
    return;
  x.check = -1;
  for (auto&& ptr : x.to)
    revertCheckTo(ptr.get());
}
void revertCheckFrom(node& x) {
  if (x.check == -1)
    return;
  x.check = -1;
  for (auto&& ptr : x.from)
    revertCheckFrom(ptr.get());
}
// compute theta for all elements in "s"
double aggregateTheta(node& x, node& y) {
  if (y.check == x.id)
    return 0.0;
  y.check = x.id;

  double theta_sum = y.theta;
  for (auto&& ptr : y.from)
    theta_sum += aggregateTheta(x, ptr.get());
  return theta_sum;
}
void computeTheta(node& x) {
  x.theta = log(x.p);
  x.check = x.id;
  for (auto&& ptr : x.from)
    x.theta -= aggregateTheta(x, ptr.get());
}
void computeThetaAll(vector<node>& s) {
  for (auto&& x : s)
    x.check = -1;
  for (auto&& x : s)
    computeTheta(x);
}

// compute eta for all elements in "s"
double aggregateEta(node& x, node& y) {
  if (y.check == x.id)
    return 0.0;
  y.check = x.id;

  double eta_sum = y.p;
  for (auto&& ptr : y.to)
    eta_sum += aggregateEta(x, ptr.get());
  return eta_sum;
}
void computeEta(node& x) {
  x.eta = x.p;
  x.check = x.id;
  for (auto&& ptr : x.to)
    x.eta += aggregateEta(x, ptr.get());
}
void computeEtaAll(vector<node>& s) {
  for (auto&& x : s)
    x.check = -1;
  for (auto&& x : reverse(s))
    computeEta(x);
}

// poset manipulations
void initializePoset(vector<node>& s) {
  for (auto&& x : s) {
    x.check = -1;
    x.is_lower = false;
    x.is_upper = false;
    x.p_init = x.p;
    x.p_prev = x.p;
    x.theta_init = x.theta;
    x.theta_prev = x.theta;
    x.eta_init = x.eta;
    x.eta_prev = x.eta;
  }
}
void revertPoset(vector<node>& s) {
  for (auto&& x : s) {
    x.check = -1;
    x.is_lower = false;
    x.is_upper = false;
    x.p = x.p_init;
    x.p_prev = x.p_init;
    x.theta = x.theta_init;
    x.theta_prev = x.theta_init;
    x.eta = x.eta_init;
    x.eta_prev = x.eta_init;
  }
}
void copyPoset(vector<node>& s, vector<node>& t, vector<double> q) {
  t.resize(s.size());
  auto iter1 = s.begin();
  auto iter2 = t.begin();
  for(; iter1 != s.end() && iter2 != t.end(); ++iter1, ++iter2) {
    iter2->id = iter1->id;
    for (auto&& p : iter1->from)
      (iter2->from).push_back(ref(t[p.get().id]));
    for (auto&& p : iter1->to)
      (iter2->to).push_back(ref(t[p.get().id]));
  }
  for (int32_t i = 0; i < (int32_t)t.size(); i++)
    t[i].p = q[i];
  computeThetaAll(t);
  computeEtaAll(t);
  initializePoset(t);
}
// compute a null distribution by copying "s"
void copyToNullPoset(vector<node>& s, vector<node>& s_null) {
  s_null.resize(s.size());
  auto iter1 = s.begin();
  auto iter2 = s_null.begin();
  for(; iter1 != s.end() && iter2 != s_null.end(); ++iter1, ++iter2) {
    iter2->id_org = iter1->id_org;
    iter2->id = iter1->id;
    iter2->p = 1.0 / (double)s.size();
    for (auto&& p : iter1->from)
      (iter2->from).push_back(ref(s_null[p.get().id]));
    for (auto&& p : iter1->to)
      (iter2->to).push_back(ref(s_null[p.get().id]));
  }
  for (auto&& x : s_null) x.theta = 0;
  s_null[0].theta = log(s_null[0].p);
  computeEtaAll(s_null);
  initializePoset(s_null);
}


// compute the mixed distribution of "s" and "t"
double aggregateDiff(node& x, node& y) {
  if (y.check == x.id || !y.is_lower)
    return 0.0;
  y.check = x.id;

  double p_diff = y.p_prev - y.p;
  for (auto&& ptr : y.to)
    p_diff += aggregateDiff(x, ptr.get());
  return p_diff;
}
void computePSingle(node& x) {
  x.p = x.p_prev;
  for (auto&& ptr : x.to)
    x.p += aggregateDiff(x, ptr.get());
}
bool computeThetaSingle(node& x_target, vector<reference_wrapper<node>>& x_target_lower, double p_new) {
  // initialization
  for (auto&& ptr : x_target_lower) {
    ptr.get().is_lower = true;
    ptr.get().check = -1;
    ptr.get().p = ptr.get().p_prev;
    ptr.get().theta = ptr.get().theta_prev;
    ptr.get().eta = ptr.get().eta_prev;
  }

  x_target.p = p_new; // set p
  // update p in the lower set (without x_target)
  for (auto&& ptr : x_target_lower)
    if (ptr.get().id != x_target.id)
      computePSingle(ptr.get());

  // check whether there is a negative probability
  bool positive = true;
  for (auto&& ptr : x_target_lower) {
    if (ptr.get().p < 0) {
      positive = false;
      break;
    }
  }

  // update theta if there is no negative probability
  if (positive) {
    for (auto&& ptr : x_target_lower)
      ptr.get().check = -1;
    for (auto&& ptr : reverse(x_target_lower))
      computeTheta(ptr.get());
  }
  // clean up
  for (auto&& ptr : x_target_lower) {
    ptr.get().check = -1;
    ptr.get().is_lower = false;
  }

  return positive;
}
double computeMixedSingle(node& x_target, double t_theta, double eps, bool revert) {
  if (fabs(x_target.theta - t_theta) < eps) return 0;

  // get the lower set by breadth first search
  vector<reference_wrapper<node>> x_target_lower;
  vector<reference_wrapper<node>> queue;
  x_target.check = x_target.id;
  queue.push_back(ref(x_target));
  while (queue.size() > 0) {
    reference_wrapper<node> tmp = queue.back();
    queue.pop_back();
    x_target_lower.push_back(tmp);
    for (auto&& ptr : tmp.get().from) {
      if (ptr.get().check == -1) {
	ptr.get().check = x_target.id;
	queue.push_back(ref(ptr.get()));
      }
    }
  }
  sort(x_target_lower.begin(), x_target_lower.end(), greaterId());

  // Find the starting probabilities for the bisection method
  double theta_b1 = x_target.theta - t_theta;
  double diff = 1.0;
  double logp_ub = log(x_target.p_prev) + diff;
  while (logp_ub > 0 || !computeThetaSingle(x_target, x_target_lower, exp(logp_ub))) {
    diff *= 0.5;
    logp_ub = log(x_target.p_prev) + diff;
  }
  double theta_b2 = x_target.theta - t_theta;
  if (theta_b1 * (theta_b2 - theta_b1) > 0) diff *= -1.0;
  while (theta_b1 * theta_b2 > 0) {
    logp_ub += diff;
    while (logp_ub > 0 || !computeThetaSingle(x_target, x_target_lower, exp(logp_ub))) {
      diff *= 0.5;
      logp_ub -= diff;
    }
    theta_b2 = x_target.theta - t_theta;
  }
  double p_lb = x_target.p_prev;
  double p_ub = exp(logp_ub);
  double theta_diff_lb = theta_b1;
  double theta_diff_ub = theta_b2;
  if (p_ub < p_lb) {
    swap(p_lb, p_ub);
    swap(theta_diff_lb, theta_diff_ub);
  }

  // Perform the bisection method
  // The solusion should exist between p_lb and p_ub
  double p_mid, theta_diff_mid;
  do {
    p_mid = (p_lb + p_ub) / 2.0;
    computeThetaSingle(x_target, x_target_lower, p_mid); // update theta
    theta_diff_mid = x_target.theta - t_theta;
    if (theta_diff_lb * theta_diff_mid > 0) {
      p_lb = p_mid;
      theta_diff_lb = theta_diff_mid;
    } else {
      p_ub = p_mid;
      theta_diff_ub = theta_diff_mid;
    }
  } while (fabs(theta_diff_mid) > eps);

  // compute the KL divergence
  double kl = 0;
  for (auto&& ptr : x_target_lower)
    kl += ptr.get().p_init * log(ptr.get().p_init / ptr.get().p);


  // clean up
  for (auto&& ptr : x_target_lower) {
    ptr.get().p_prev = ptr.get().p;
    ptr.get().theta_prev = ptr.get().theta;
    ptr.get().eta_prev = ptr.get().eta;
  }

  // revert to the initial state
  if (revert) {
    for (auto&& ptr : x_target_lower) {
      ptr.get().check = -1;
      ptr.get().is_lower = false;
      ptr.get().is_upper = false;
      ptr.get().p = ptr.get().p_init;
      ptr.get().p_prev = ptr.get().p_init;
      ptr.get().theta = ptr.get().theta_init;
      ptr.get().theta_prev = ptr.get().theta_init;
      ptr.get().eta = ptr.get().eta_init;
      ptr.get().eta_prev = ptr.get().eta_init;
    }
  }

  return kl;
}
void computeMixed(vector<reference_wrapper<node>>& s_nodes, vector<node>& s, vector<node>& t, double eps) {
  double theta_prev, theta = s[0].theta;
  do {
    theta_prev = theta;
    for (auto&& ptr : s_nodes) {
      computeTheta(ptr.get());
      revertCheckFrom(ptr.get());
      computeMixedSingle(ptr.get(), t[ptr.get().id].theta, eps, false);
    }
    theta = s[0].theta;
  } while (fabs(theta - theta_prev) > eps);

  computeEtaAll(s);
}
