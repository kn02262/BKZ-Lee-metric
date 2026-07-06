#load("../LEE.sage")
from LEE import *
from bkz import *
import sys
import argparse
import numpy
import time

# Read input parameters
parser = argparse.ArgumentParser(description='This script compares efficiency of different BKZ modes for codes over F_q')
parser.add_argument("-q", required=True, type=int, default=0, help='Finite field cardinality, prime')
parser.add_argument("-b", required=True, type=int, default=0, help='Beta block size in BKZ')
parser.add_argument("-s", required=True, type=int, default=10, help='Number of samples')
parser.add_argument("-n", required=True, type=int, default=0, help='Length of the code')
parser.add_argument("--seed", dest="seed", required=False, type=int, default=0, help="Random seed: default is 0=do not set")
parser.add_argument("--ignore_potential_growth", dest="ignore_potential_growth", required=False, action=argparse.BooleanOptionalAction, default=False, help="Ignore potential function growth in LLL")
parser.add_argument("--verbose", dest="verbose", required=False, action=argparse.BooleanOptionalAction, default=False, help="Verbose output, print raw data")
parser.add_argument("--algorithm", dest="algorithm", required=False, type=int, default=0, help="Algorithm: 0=Classic BKZ, 1=Classic bkz, selecting lower LeeWt Vec., 2=bkz with minimum_lee_wt_codeword+primitive, 3=LeeBKZ, minimum_lee_wt_codeword may be non-primitive")

args = parser.parse_args()
beta = args.b
q = args.q
samples = args.s
if args.algorithm == 0:
    alg_caption_string = "\n--- BKZ with beta=" + str(beta) + ", classic ---"
elif args.algorithm == 1:
    alg_caption_string = "\n--- BKZ with beta=" + str(beta) + ", selecting lowest lee wt. vec ---"
elif args.algorithm == 2:
    alg_caption_string = "\n--- LeeBKZ with beta=" + str(beta) + ", bkz with minimum_lee_wt_codeword+primitive ---"
elif args.algorithm == 3:
    alg_caption_string = "\n--- LeeBKZ with beta=" + str(beta) + ", bkz with minimum_lee_wt_codeword (may be non-primitive) ---"
else:
    print("ERROR: Invalid algorithm selected, run LeeBKZ.sage --help")
    assert false

F=GF(q)
n=args.n
R=0.5
if args.seed != 0:
    print(f"RANDOM_SEED={args.seed}")
    set_random_seed(args.seed)
else:
    print("RANDOM_SEED is not set")

allstat = []
wtH0 = [] # B[0].hamming_wt (proper basis)
wtL0 = [] # B[0].lee_wt (proper basis)
wtH = [] # B[0].hamming_wt (reduced basis)
wtL = [] # B[0].lee_wt (reduced basis)
T = []
MCalls = []
l1_fix_fails = []

k=int(n*R)
for s in range(samples):
    print()
    while true:
        B = random_matrix(F, k, n)
        if B[:,:k].rank() == k:
            break
    assert B.rank() == B.nrows()
    B = proper_basis(B)
    if args.verbose:
        print(B, flush=true)
    B_profile = ell_profile(B)
    print(f"profile: {B_profile}")
    print(f"{B[0].hamming_weight()}<=B[0].wtL<={B[0].hamming_weight()*floor(F.cardinality()/2)}")
    print(f"B[0].wtH={B[0].hamming_weight()}, B[0].wtL={WtLeeVec(B[0])}")
    wtH0.append(B[0].hamming_weight())
    wtL0.append(WtLeeVec(B[0]))

    # Start algorithm
    print(alg_caption_string)
    STAT=[0,0]
    t=time.time()
    if (args.algorithm == 0) or (args.algorithm == 1):
        B_red=bkz(B, beta, args.algorithm, verbose=args.verbose, STAT=STAT)
    elif (args.algorithm == 2):
        B_red=LEEbkz(B, beta, primitive=True, verbose=args.verbose, ignore_potential=args.ignore_potential_growth, STAT=STAT)
    elif (args.algorithm == 3):
        B_red=LEEbkz(B, beta, primitive=False, verbose=args.verbose, ignore_potential=args.ignore_potential_growth, STAT=STAT)
    T.append(time.time()-t)
    B_red_profile = ell_profile(B_red)
    print(f"profile after BKZ: {B_red_profile}")
    print(f"{B_red[0].hamming_weight()}<=B_red[0].wtL<={B_red[0].hamming_weight()*floor(F.cardinality()/2)}")
    print(f"B_red[0].wtH={B_red[0].hamming_weight()}, B_red[0].wtL={WtLeeVec(B_red[0])}")
    wtH.append(B_red[0].hamming_weight())
    wtL.append(WtLeeVec(B_red[0]))
    MCalls.append(STAT[0])
    l1_fix_fails.append(STAT[1])


print("\n\n****************************************************")
print(f"******** STATISTICS OVER N={args.s} samples ********")
print("****************************************************")
print(f"Reduced {samples} matrices of size {k}x{n} over F{q}")
print(f"Proper basis wtH[0]={numpy.mean(wtH0)} wtL[0]={numpy.mean(wtL0)}")
print(f"Reduced basis wtH[0]={numpy.mean(wtH)} wtL[0]={numpy.mean(wtL)}")
print(f"Avg.time={numpy.mean(T)}")
print(f"Avg.Mw_codeword.calls={numpy.mean(MCalls)}")
if beta == 2:
    print(f"Avg.potential_growth_fails={numpy.mean(l1_fix_fails)}")
