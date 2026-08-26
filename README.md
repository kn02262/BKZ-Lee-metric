# BKZ-Lee-metric
BKZ basis reduction for linear codes in Lee metric.

* Generic BKZ algorithm for linear codes in non-Hamming metric.
* Statistics on the low Lee-weight codewords obtained using reduction.

System requirements:
- Sagemath 10.7 (currently does not work on later versions)
- gcc, cython

Before usage, please compile .so modules by running:
> make

usage:
> LeeBKZ.sage.py [-h] -q Q -b B -s S -n N [--seed SEED] [--ignore_potential_growth | --no-ignore_potential_growth]
                      [--verbose | --no-verbose] [--algorithm ALGORITHM]

This script compares efficiency of different BKZ modes for codes over $F_q$.

Options:
```
  -h, --help            show this help message and exit
  -q Q                  Finite field cardinality, prime
  -b B                  Beta block size in BKZ
  -s S                  Number of samples
  -n N                  Length of the code
  --seed SEED           Random seed: default is 0=do not set
  --verbose, --no-verbose
                        Verbose output, print raw data
  --algorithm ALGORITHM
                        Algorithm: 0=Classic BKZ, 1=Classic bkz, selecting lower LeeWt Vec., 2=bkz with
                        minimum_lee_wt_codeword+primitive, 3=LeeBKZ, minimum_lee_wt_codeword may be non-primitive
                        4=LeeLLL, LeeBKZ with no potential growth
```

Script generates -s samples of random (full-rank) matrices of size $k \times n$ over finite field $F_q$ with $k = n/2$.

Matrices are being BKZ-reduced (with block size -b) in 4 modes (select with --algorithm):

1. HammingBKZ.
2. HammingBKZ + Least Lee-weight codeword selection.
3. LeeBKZ + Proper bases.
4. LeeBKZ.
5. LeeLLL = LeeBKZ for block size beta=2, polynomial complexity.

For selected mode the script computes mean value of the found low Lee-weight codewords.

See paper "Block basis reduction for linear codes in non-Hamming metrics" for details.

To compute Gilbert-Varshamov bounds for Lee and Hamming metrics use:
> sage bounds.sage

To reproduce experiments from the paper, run:
> parallel -j <NUMBER_OF_AVAILABLE_CORES> :::: experiment_tasks_all.txt
