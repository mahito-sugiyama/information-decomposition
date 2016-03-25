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

#include "mix.h"
#include "dataprocess.h"
#include <time.h>
#include <unistd.h>

using namespace std;

struct greaterScore { bool operator()(const node& left, const node& right) const { return left.score > right.score; } };

int main(int argc, char *argv[]) {
  char *input_file = NULL;
  char *output_file = NULL;
  bool verbose = false;
  bool flag_in = false, flag_out = false;
  int32_t item_length_lb = 1;
  int32_t topk = 10;
  double sigma = 0.1;
  double eps_exponent = 5.0; // error
  double alpha = 0.05;

  clock_t time_total_start = clock();

  // get arguments
  char opt;
  while ((opt = getopt(argc, argv, "i:o:s:e:l:k:a:v")) != -1) {
    switch ( opt ) {
    case 'i': input_file = optarg; flag_in = true; break;
    case 'o': output_file = optarg; flag_out = true; break;
    case 's': sigma = atof(optarg); break;
    case 'e': eps_exponent = atof(optarg); break;
    case 'l': item_length_lb = atoi(optarg); break;
    case 'k': topk = atoi(optarg); break;
    case 'a': alpha = atof(optarg); break;
    case 'v': verbose = true; break;
    }
  }
  double eps = pow(10.0, -1 * eps_exponent);

  vector<vector<int32_t>> data;
  vector<int32_t> item_list;

  // read a database
  // ifstream ifs("T10I4D100K.dat");
  if (!flag_in) {
    cerr << "> ERROR: Input file (-i [input_file]) is not specified!" << endl;
    exit(1);
  }
  cout << "> Reading a database file \"" << input_file << "\" ... " << flush;
  ifstream ifs(input_file);
  readDatabase(ifs, data, item_list);
  cout << "end" << endl << flush;
  cout << "  Information:" << endl << flush;
  cout << "  Number of transactions: " << data.size() << endl << flush;
  cout << "  Number of items:        " << item_list.size() << endl << flush;

  // check duplication and compute the support for each transaction
  // keep transactions whose support is larger than sigma
  cout << "> Computing the support" << endl << flush;
  vector<node> s;
  computeSupport(data, sigma, s);
  cout << "  Information:" << endl << flush;
  cout << "  Number of frequent combinations: " << s.size() << endl << flush;

  if (s.size() <= 1) {
    cerr << "> ERROR: The threshold " << sigma << " is too large!" << endl;
    exit(1);
  }

  if (verbose) {
    for (auto&& x : s) {
      cout << "  " << x.id_org << ": ";
      for (auto&& t : data[x.id_org]) cout << t << " ";
      cout << "(" << x.supp << ")" << endl;
    }
  }

  // make a poset
  cout << "> Constructing a poset ... " << flush;
  int32_t num_edges = computePoset(s, data);
  cout << "end" << endl << flush;

  cout << "  Information:" << endl << flush;
  cout << "  Number of edges: " << num_edges << endl << flush;

  if (verbose) {
    cout << "  Structure:" << endl;
    for (auto&& x : s) {
      for (auto&& p : x.to) {
	cout << "  " << x.id_org << " -> " << p.get().id_org << endl;
      }
    }
  }

  cout << "> Conputing theta ... " << flush;
  computeThetaAll(s);
  cout << "end" << endl << flush;
  cout << "> Conputing eta ... " << flush;
  computeEtaAll(s);
  cout << "end" << endl << flush;
  cout << "> Initializing poset ... " << flush;
  initializePoset(s);
  cout << "end" << endl << flush;

  if (verbose) {
    cout << "  Distribution:" << endl << s;
  }

  // compute feature scores (KL divergence)
  cout << "> Conputing feature scores (KL divergence)" << endl << flush;
  clock_t time_score_start = clock();
  computeFeatureScore(s, eps);
  clock_t time_score_end = clock();
  cout << "  Runtime for computing scores: " << (float)(time_score_end - time_score_start) / CLOCKS_PER_SEC << endl;
  // compute p-values
  cout << "> Conputing p-values ... " << flush;
  computePvalues(s, (int32_t)data.size());
  cout << "end" << endl << flush;

  // compute the number of significant combinations for each size
  sort(s.begin(), s.end(), greaterScore()); // Note: Sorting will destroy the structure of "s"!!
  vector<int32_t> sig_num;
  computeEachSignum(s, data, alpha, sig_num);
  cout << "  Information:" << endl << flush;
  cout << "  Number of significant patterns: " << accumulate(sig_num.begin(), sig_num.end(), 0) << endl;
  for (int32_t i = 1; i < (int32_t)sig_num.size(); i++) {
    if (sig_num[i] > 0) cout << "  size " << i << ": " << sig_num[i] << endl;
  }

  // write all significant combinations to "output_file"
  if (flag_out) {
    cout << "> Writing all significant combinations to \"" << output_file << "\" ... " << flush;
    ofstream ofs(output_file);
    outputSig(s, data, alpha, item_length_lb, -1, 0, ofs);
    ofs.close();
    cout << "end" << endl << flush;
  }

  // output top-k combinations
  cout << "> Top-" << topk << " combinations (KL & Bonferroni-corrected p-value):" << endl;
  outputSig(s, data, alpha, item_length_lb, topk, 2, cout);

  clock_t time_total_end = clock();
  cout << "> Total runtime: " << (float)(time_total_end - time_total_start) / CLOCKS_PER_SEC << endl;
}
