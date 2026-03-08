# BKZ-Lee-metric
BKZ basis reduction for linear codes in Lee metric

Generic BKZ algorithm for linear codes in non-Hamming metric.
Statistics on the lowest Lee Weight codeword generated with reduction.

Useage:
> sage LeeBKZ.sage -q 127 -b 2 -s 10

Script generates -s samples of random (full rank) matrices of size k*n over finite field F_q
with n = 32,64,128,256,512,1024
k = n/2
Matrices are being BKZ-reduced (with block size -b) in 4 modes:

(1) BKZ
(2) BKZ + Lower Lee wt.
(3) LeeBKZ + Primitivity
(4) LeeBKZ

For each mode computes mean value of the lowest Lee weight.

See paper "Block basis reduction for linear codes in non-Hamming metrics" for details.

Compute Gilbert-Varshamov bounds for Lee and Hamming metrics:
> sage bounds.sage
