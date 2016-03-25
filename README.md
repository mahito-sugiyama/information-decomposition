# Information Decomposition
An implementation of information decomposition algorithms on structured space.  
Please see the following paper for more details:
* Sugiyama, M., Nakahara, H., Tsuda, K.: **Information Decomposition on Structured Space**, *[arXiv:1601.05533](http://arxiv.org/abs/1601.05533)* (2016)

## Usage
The file "mix.h" contains the source code of mixing two distributions on a poset that achieves the KL divergence decomposition.
Please check "example.cc", where you can see how to use the functions in "mix.h".

We also provide data processing functions in "dataprocess.h" that construct a poset from a transaction database and compute the score (KL divergence) and the *p*-value for each combination of features (attributes).

The code is written in C++11.
To compile the code including "dataprocess.h", the Boost library is needed, which is used to compute the *p*-values from the *&chi;*<sup>2</sup> distribution.
Please modify the "Makefile" according to your environment.

The following is an example in the terminal (in the directory "src"):
```
$ make
$ ./decomp -i ../data/T10I4D100K.data -o output -s 0.0001
> Reading a database file "../data/T10I4D100K.data" ... end
  Information:
  Number of transactions: 100000
  Number of items:        870
> Computing the support
  [====================================================================================================] 100 %
  Information:
  Number of frequent combinations: 364
> Constructing a poset ... end
  Information:
  Number of edges: 363
> Conputing theta ... end
> Conputing eta ... end
> Initializing poset ... end
> Conputing feature scores (KL divergence)
  [===================================================================================================>] 100 %
  Runtime for computing scores: 0.0496
> Conputing p-values ... end
  Information:
  Number of significant patterns: 359
  size 2: 1
  size 3: 20
  size 4: 58
  size 5: 93
  size 6: 91
  size 7: 55
  size 8: 24
  size 9: 9
  size 10: 6
  size 11: 1
  size 13: 1
> Writing all significant combinations to "output" ... end
> Top-10 combinations (KL & Bonferroni-corrected p-value):
  351 477 766 820, 0.640463, 0
  275 305 368 471, 0.640463, 0
  72 234 343 427 528 897, 0.640463, 0
  0 140 239 257 368 871, 0.640463, 0
  58 140 190 354 460 859, 0.640463, 0
  38 319 522 588, 0.640463, 0
  39 170 438 500 575 738, 0.640463, 0
  120 181 333 798 998, 0.640463, 0
  73 334 474 638 989, 0.640463, 0
  95 159 517 722 966, 0.640463, 0
> Total runtime: 0.490335
```

#### Command-line arguments

  `-i <input_file>` : An input file of a transaction database  
  `-o <output_file>` : An output file of the full list of significant combinations of features  
  `-s <prob_threshold>` : A threshold for probability of combinations that are included in the resulting poset    
  `-l <lower_bound>` : The lower bound of the size of combinations (default: 1)  
  `-a <alpha>` : The significance level (default: 0.05)  

## Contact
Author: Mahito Sugiyama  
Affiliation: ISIR, Osaka University, Japan  
E-mail: mahito@ar.sanken.osaka-u.ac.jp
