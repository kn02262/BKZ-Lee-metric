# BKZ-Lee-metric
BKZ basis reduction for linear codes in Lee metric.

* Generic BKZ algorithm for linear codes in non-Hamming metric.
* Statistics on the low Lee-weight codewords obtained using reduction.

Usage:
> sage LeeBKZ.sage -q 127 -b 2 -s 10

Script generates -s samples of random (full-rank) matrices of size $k \times n$ over finite field $F_q$ with $n = 32,64,128,256,512,1024$, $k = n/2$.

Matrices are being BKZ-reduced (with block size -b) in 4 modes:

1. HammingBKZ.
2. HammingBKZ + Least Lee-weight codeword selection.
3. LeeBKZ + Proper bases.
4. LeeBKZ.

For each mode the script computes mean value of the found low Lee-weight codewords.

See paper "Block basis reduction for linear codes in non-Hamming metrics" for details.

To compute Gilbert-Varshamov bounds for Lee and Hamming metrics use:
> sage bounds.sage
