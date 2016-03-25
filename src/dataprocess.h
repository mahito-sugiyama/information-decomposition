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

#include <boost/math/distributions/chi_squared.hpp>

// for index sort
int32_t sort_index;

// sort a 2D-vector according to the size of vectors
struct lessSize { bool operator()(const vector<int32_t>& left, const vector<int32_t>& right) const { return left.size() < right.size(); } };
struct lessIndex { bool operator()(const vector<int32_t>& left, const vector<int32_t>& right) const { return left[sort_index] < right[sort_index]; } };

// show a progress bar
void progress(double progress) {
  int32_t barWidth = 100;

  cerr << "  [";
  int32_t pos = barWidth * progress;
  for (int32_t i = 0; i < barWidth; i++) {
    if (i < pos) cerr << "=";
    else if (i == pos) cerr << ">";
    else cerr << " ";
  }
  cerr << "] " << ceil(progress * 100.0) << " %\r";
  cerr.flush();
}

// TRUE if [trans1] is included in [trans2]
// Both [trans1] and [trans2] should be sorted beforehand
bool included(vector<int32_t>& trans1, vector<int32_t>& trans2) {
  if (trans1.size() > trans2.size()) return false;

  int32_t counter = 0;
  bool match;
  for (auto&& x : trans1) {
    match = false;
    for (; counter < (int32_t)trans2.size(); counter++) {
      if (x == trans2[counter]) {
	match = true;
	break;
      } else if (x < trans2[counter]) {
	break;
      }
    }
    if (!match)
      break;
  }
  return match;
}

// read a database file
void readDatabase(ifstream& ifs, vector<vector<int32_t>>& data, vector<int32_t>& item_list) {
  int32_t buf;
  string str;
  while (getline(ifs, str)) {
    stringstream sstr(str);
    vector<int32_t> tmp;
    while (sstr >> buf) {
      tmp.push_back(buf);
      item_list.push_back(buf);
    }
    sort(tmp.begin(), tmp.end());
    data.push_back(tmp);
  }
  // get the list of items
  sort(item_list.begin(), item_list.end());
  item_list.erase(unique(item_list.begin(), item_list.end()), item_list.end());
}

// compute the support for each transaction
void computeSupport(vector<vector<int32_t>>& data, double sigma, vector<node>& s) {
  // sort a dataset according to the length of transactions
  sort(data.begin(), data.end(), lessSize());
  vector<int32_t> idx_length_change;
  idx_length_change.push_back(0);
  // get indices where the length of a transaction changes
  for (int32_t i = 1; i < (int32_t)data.size(); i++) {
    if (data[i].size() != data[i - 1].size())
      idx_length_change.push_back(i);
  }
  idx_length_change.push_back(data.size());
  // sort a dataset
  for (int32_t i = 1; i < (int32_t)idx_length_change.size(); i++) {
    int32_t idx_str = idx_length_change[i - 1];
    int32_t idx_end = idx_length_change[i];
    int32_t current_size = data[idx_str].size();
    for (sort_index = 0; sort_index < current_size; sort_index++)
      stable_sort(data.begin() + idx_str, data.begin() + idx_end, lessIndex());
    progress((double)idx_length_change[i] / (double)data.size());
  }
  cerr << endl;
  reverse(data.begin(), data.end());


  // compute the support
  s.push_back(node());
  s.back().id_org = -1; s.back().supp = -1;
  bool checking = false;
  double thr = sigma * (double)data.size();
  for (int32_t i = 0; i < (int32_t)data.size() - 1; i++) {
    if (!checking) {
      if ((double)(s.back()).supp < thr) s.pop_back();
      s.push_back(node());
      s.back().id_org = i;
      s.back().supp = 1;
    }
    checking = false;
    if (data[i].size() == data[i + 1].size()) {
      int32_t counter;
      for (counter = 0; counter < (int32_t)data[i].size(); counter++) {
	if (data[i][counter] != data[i + 1][counter])
	  break;
      }
      if (counter == (int32_t)data[i].size()) { // match
	checking = true;
	(s.back().supp)++;
      }
    }
  }
  // the last transaction
  if ((double)(s.back()).supp < thr) s.pop_back();
  if (!checking && 1.0 >= thr) {
    s.push_back(node());
    s.back().id_org = (int32_t)data.size() - 1;
    s.back().supp = 1;
  }
  // if there is no bottom
  if (data.back().size() > 0) {
    // computing the bottom
    int32_t bottom_supp = 0;
    for (auto&& x : s)
      bottom_supp += x.supp;
    bottom_supp = data.size() - bottom_supp; // the support of the bottom
    if (bottom_supp == 0) {
      cout << "end" << endl << flush;
      cerr << "> ERROR: The threshold " << sigma << " is too small!" << endl;
      exit(1);
    }
    // add the bottom
    vector<int32_t> bottom;
    data.push_back(bottom);
    s.push_back(node());
    s.back().id_org = (int32_t)data.size() - 1;
    s.back().supp = bottom_supp;
  }
  // compute probability with normalization
  int32_t sum_supp = 0;
  for (auto&& x : s)
    sum_supp += x.supp;
  for (auto&& x : s)
    x.p = (double)x.supp / (double)sum_supp;

  reverse(s.begin(), s.end());
  for (int32_t i = 0; i < (int32_t)s.size(); i++)
    s[i].id = i;
}

// compute the structure of a poset
int32_t computePoset(vector<node>& s, vector<vector<int32_t>>& data) {
  for (int32_t i = 0; i < (int32_t)s.size(); i++) {
    for (int32_t j = i - 1; j >= 0; j--) s[j].is_lower = false;
    for (int32_t j = i - 1; j >= 0; j--) {
      if (!s[j].is_lower && included(data[s[j].id_org], data[s[i].id_org])) {	
	s[i].from.push_back(ref(s[j]));
	s[j].to.push_back(ref(s[i]));
	getLowerNext(s[j]);
      }
    }
  }
  // compute the number of edges
  int32_t num_edges = 0;
  for (auto&& x : s) num_edges += x.to.size();
  return num_edges;
}

// compute the KL divergence for each node
void computeFeatureScore(vector<node>& s, double eps) {
  for (auto&& x : s) {
    x.score = x.id > 0 ? computeMixedSingle(x, 0, eps, true) : 0;
    progress((double)x.id / (double)s.size());
  }
  cerr << endl;
}

// compute p-values
void computePvalues(vector<node>& s, int32_t data_size) {
  boost::math::chi_squared chisq_dist((int32_t)s.size() - 1);
  for (auto&& x : s) {
    if (x.score <= pow(10, -8)) x.pvalue = 1;
    else x.pvalue = (1 - boost::math::cdf(chisq_dist, 2 * data_size * x.score)) * s.size();
    if (x.pvalue > 1) x.pvalue = 1;
  }
}

// compute the number of significant combinations for each size
void computeEachSignum(vector<node>& s, vector<vector<int32_t>>& data, double alpha, vector<int32_t>& sig_num) {
  int32_t sig_size_max = 0;
  for (auto && x : s) {
    if (x.pvalue > alpha) break;
    if ((int32_t)data[x.id_org].size() > sig_size_max) sig_size_max = data[x.id_org].size();
  }
  sig_num.resize(sig_size_max + 1);
  fill(sig_num.begin(), sig_num.end(), 0);
  for (auto && x : s) {
    if (x.pvalue > alpha) break;
    sig_num[data[x.id_org].size()]++;
  }
}

// write significant combinations to "ofs"
void outputSig(vector<node>& s, vector<vector<int32_t>>& data, const double alpha, const int32_t lb, const int32_t topk, const int32_t sp, ostream &ofs) {
  int32_t counter = 1;
  for (auto && x : s) {
    if (x.pvalue > alpha) break;
    if (x.id > 0 && (int32_t)data[x.id_org].size() >= lb) {
      for (int32_t i = 0; i < sp; i++)
	ofs << " ";
      for (int32_t i = 0; i < (int32_t)data[x.id_org].size() - 1; i++)
	ofs << data[x.id_org][i] << " ";
      ofs << data[x.id_org].back() << ", " << x.score << ", " << x.pvalue << endl;
      if (topk > 0 && ++counter > topk)
	break;
    }
  }
}
