#include "mix.h"

using namespace std;

int main(int argc, char *argv[]) {
  /*
    make a poset "s" such that
      3
     / \
    1   2
     \ /
      0
    with probabilities (0.1, 0.4, 0.2, 0.3)
  */
  vector<node> s;
  s.resize(4);
  // set "id" and "p" (probabilities)
  vector<double> p = {0.1, 0.4, 0.2, 0.3};
  for (int32_t i = 0; i < (int32_t)s.size(); i++) {
    s[i].id = i;
    s[i].p = p[i];
  }
  // set orders (edges)
  s[0].to.push_back(ref(s[1]));
  s[0].to.push_back(ref(s[2]));
  s[1].to.push_back(ref(s[3]));
  s[2].to.push_back(ref(s[3]));
  s[1].from.push_back(ref(s[0]));
  s[2].from.push_back(ref(s[0]));
  s[3].from.push_back(ref(s[1]));
  s[3].from.push_back(ref(s[2]));
  // compute theta
  computeThetaAll(s);
  // compute eta
  computeEtaAll(s);
  // initialization
  initializePoset(s);

  // copy a poset "s" to "t"
  vector<node> t;
  vector<double> q = {0.01, 0.1, 0.8, 0.09};
  copyPoset(s, t, q);

  cout << "Distribution s:" << endl;
  cout << s << endl;

  cout << "Distibution t:" << endl;
  cout << t << endl;

  // compute the KL divergence between "s" and "t"
  double kl_s_t = 0;
  for (int32_t i = 0; i < (int32_t)s.size(); i++)
    kl_s_t += s[i].p * log(s[i].p / t[i].p);

  // mix "s" and "t" with respect to {1, 2}
  vector<reference_wrapper<node>> s_nodes = {ref(s[1]), ref(s[2])};
  double eps = pow(10.0, -1 * 10.0); // error bound for numerical optimization
  computeMixed(s_nodes, s, t, eps);

  // output the mixed distribution
  cout << "Mixed distibution of (s, t) w.r.t. {1, 2}:" << endl;
  cout << s << endl;

  // compute the KL divergence of "s -> s_mixed" and "s_mixed -> t"
  double kl_s_smixed = 0;
  for (int32_t i = 0; i < (int32_t)s.size(); i++)
    kl_s_smixed += s[i].p_init * log(s[i].p_init / s[i].p);
  double kl_smixed_t = 0;
  for (int32_t i = 0; i < (int32_t)s.size(); i++)
    kl_smixed_t += s[i].p * log(s[i].p / t[i].p);

  // check the decomposition
  cout << "KL(s, t) = " << kl_s_t << endl;
  cout << "KL(s, mixed) + KL(mixed, t) = " << kl_s_smixed << " + " << kl_smixed_t << " = " << kl_s_smixed + kl_smixed_t << endl;
}
