load("LEE.sage")
import sys
import argparse
import numpy
import time

def WtLee(x): # Coordinate Lee weight
    F=x.base_ring()
    b=floor(F.cardinality()/4)
    c=floor(F.cardinality()/2)
    if x>c:
        x=F.cardinality()-x
    return x.lift()

def WtLeeVec(x): # Vector Lee weight
    return sum(WtLee(s) for s in x)

def proj_orthogonal_lee_vector_q(a,b):
     # Projection of b orthogonal to a (Lee)
     res = []
     for i in range(len(a)):
         if WtLee(a[i])+WtLee(b[i]) == WtLee(a[i]+b[i]):
             res.append(b[i])
         else:
             res.append(0)
     return vector(res)

def matrix_weight(B):
    c = 0
    for i in range(B.nrows()):
        c += B[i].hamming_weight()
    return c

def intersection_basis(S, B):
    # Implementation of the method from Lemma 1 of the paper:
    # Micciancio D. - Efficient reductions among lattice problems (2007)
    #
    S = matrix(S)
    H = matrix(S.transpose().left_kernel().basis())
    assert H.ncols() == B.ncols(), f"{H.ncols()} != {B.ncols()}"
    C = H * B.transpose()
    #C = C.stack(zero_matrix(ZZ, B.nrows() - H.nrows(), B.nrows()))
    sD,sU,sV = C.smith_form(transformation=True)
    assert sD == sU * C * sV
    B_t = B.transpose() * sV
    return B_t.transpose()[::-1,:]

def intersection_basis_reduced(S, B):
    S = matrix(S)
    assert S.nrows() <= B.nrows()
    B0 = intersection_basis(S, B)
    if S.nrows() == B.nrows():
        return B0
    for i in range(S.nrows(), B.nrows()):
        Bt = B0[i:]
        C = LinearCode(Bt)
        p = C._minimum_weight_codeword()
        S = S.stack(matrix(p))
        B0 = intersection_basis(S, B)
    return B0

def proper_basis(B):
    k = B.nrows()
    A = B[:,:k]
    assert A.ncols() == k and A.nrows() == k
    return A^(-1) * B

def or_vector(a,b):
    return vector([a[i] or b[i] for i in range(len(a))])

def proj_orthogonal_vector(a,b):
    return vector([(a[i]+1)*b[i] for i in range(len(a))])

def proj_orthogonal_vector_q(a,b):
     # Projection of b orthogonal to a
     res = []
     for i in range(len(a)):
         if a[i] != 0:
             res.append(0)
         else:
             res.append(b[i])
     return vector(res)

def orthogonal_projections(B, i, j, proj_vector=proj_orthogonal_vector_q):
    k = B.nrows()
    assert i < k and i >= 0 and i <= j
    assert j >= 0
    if i == 0:
        return matrix(B.base_ring(), B[i:j+1,:])
    B_res = []
    for ii in [i..j]:
        p = B[ii]
        for jj in [0..i-1]:
            p = proj_vector(B[jj], p)
        B_res.append(p)
    return matrix(B_res)

def insert_primitive(B, p, proj_vector=proj_orthogonal_vector_q):
    k = B.nrows()
    B = copy(B)
    B0 = copy(B)
    a = B.solve_left(p)
    assert a * B == p
    m = infinity
    for i in range(0,B.nrows()):
        if a[i] != 0:
            m = i
            break
    B[m] = p
    B.swap_rows(0, m)
    B_proj = orthogonal_projections(B, 1, k - 1, proj_vector)
    B_proj_ech = B_proj.echelon_form()
    T = B_proj.solve_left(B_proj_ech)
    assert T.nrows() == k-1
    assert T.ncols() == k-1
    T1 = block_diagonal_matrix(identity_matrix(1), T)
    B = T1 * B
    A = B0.solve_left(B)
    #print("Transformation A:", A.ncols(), A.nrows(), A.rank())
    assert (A*B0).rank() == k
    assert (A*B0)[0] == p
    return A

def insert_vector(B, p):
    F=B.base_ring()
    k = B.nrows()
    B = copy(B)
    a = B.solve_left(p)
    assert a * B == p
    m = infinity
    for i in range(0,B.nrows()):
        if a[i] != 0:
            m = i
            break
    A = identity_matrix(F, k);
    b0=A[0]
    A[0]=a
    if m>0:
        A[m]=b0
    #print("Transformation: ", A.ncols(), A.nrows(), A.rank(), A)
    #print("---")
    return A

def insert_primitive_global(B, p):
    A = insert_primitive(B, p)
    # Apply transformation A to matrix B
    #T = block_diagonal_matrix(identity_matrix(i).change_ring(GF(2)), A.change_ring(GF(2)), identity_matrix(k-j-1).change_ring(GF(2)))
    B = A * B
    return B

def MakePrimitive(B, c): # Definition: Claim 4.7 [More basis reduction for linear codes.]
    S=[] # Supp(c)
    Sb=[] # [n]-Supp(c)
    for i in [0..len(c)-1]:
        if c[i] == 0:
            Sb.append(i)
        else:
            S.append(i)
    B_sub = B[ range(B.nrows()), [ k for k in [0..B.ncols()-1] if k in Sb ] ]
    B_sub_pivots = B_sub.pivots()
    m=len(B_sub_pivots)
    B1_pivots = []
    for i in B_sub_pivots:
        B1_pivots.append(Sb[i]+1)
    B1_pivots = B1_pivots + [k for k in [1..B.ncols()] if k not in B1_pivots]
    # Permute columns of B, B1_pivots come first
    P=Permutation(B1_pivots).to_matrix().change_ring(B.base_ring())
    B1=(B*P).echelon_form()
    cp = matrix(B1[m])*P^(-1)
    return vector(cp[0])

def MW_Codeword_Hamming_fast(B, lee=false):
    zero_pos = []
    Bn=B.ncols()
    for j in range(Bn):
        col = B.column(j)
        if col.is_zero():
            zero_pos.append(j)
    B=B.delete_columns(zero_pos)
    C=LinearCode(B)
    if lee == false:
        p = C._minimum_weight_codeword(algorithm='guava')
    else:
        p=C[0]
        dw=C.length()
        dl=C.length()*floor(B.base_ring().cardinality()/2)
        for c in C:
            if c.hamming_weight() != 0 and c.hamming_weight()<dw:
                dw=c.hamming_weight()
                dl=WtLeeVec(c)
                p=c
            elif c.hamming_weight != 0 and c.hamming_weight()==dw:
                if WtLeeVec(c) < dl:
                    dl=WtLeeVec(c)
                    p=c
    res = []
    ptr = 0
    for i in range(Bn):
        if i in zero_pos:
            res.append(0)
        else:
            res.append(p[ptr])
            ptr = ptr + 1
    return vector(res)


def MW_Codeword(B, lee=0):
    if lee == 0:
        p = MW_Codeword_Hamming_fast(B,false)
    elif lee == 1:
        p = MW_Codeword_Hamming_fast(B,true)
    elif lee == 2:
        if B.rank() == 0:
            p = zero_vector(F, B.ncols())
        else:
            C=LinearLeeMetricCode(B)
            p = C.minimum_lee_weight_codeword_fast(primitive=true)
    elif lee == 3:
        if B.rank() == 0:
            p = zero_vector(F, B.ncols())
        else:
            C=LinearLeeMetricCode(B)
            p = C.minimum_lee_weight_codeword_fast()
    return p

def bkz(B, beta, preproc=[], lee=0, proj_vector=proj_orthogonal_vector_q):
    # lee=0: hamming bkz
    # lee=1: hamming bkz + lower WtLee codeword selection
    # lee=2: lee bkz (inserting primitive codeword)
    # lee=3: lee bkz (inserting codeword may be non-primitive)

    # proj_vector=proj_orthogonal_vector_q: Projection with non-intersecting supports
    # proj_vector=proj_orthogonal_lee_vector_q: Projection with orthogonal elements
    k = B.nrows()
    assert 2 <= beta and beta <= k
    i = 0
    while i < k - 1:
        j = min(i+beta-1, k-1)
        B_proj = orthogonal_projections(B, i, j, proj_vector)
        p = MW_Codeword(B_proj, lee)
        if (lee>=2 and WtLeeVec(B_proj[0]) == WtLeeVec(p)) or (lee<=1 and B_proj[0].hamming_weight() == p.hamming_weight()):
            i = i + 1
        else:
            if (lee <= 2) and (proj_vector == proj_orthogonal_vector_q):
                A = insert_primitive(B_proj, p, proj_vector)
            else:
                A = insert_vector(B_proj, p)
            T = block_diagonal_matrix(identity_matrix(i), A, identity_matrix(k-j-1))
            E = T * B
            if ((WtLeeVec(E[i]) > WtLeeVec(B[i])) and (proj_vector != proj_orthogonal_vector_q)):
                i=i+1
                continue
            B = E
            i = max(0, i-beta+1)
    return B

def epipodal_vector(B, i):
    # Build a proj system
    p = vector([0 for x in range(B.ncols())])
    for i0 in range(i):
        p = or_vector(p, B[i0])
    return proj_orthogonal_vector_q(p,B[i])

def epipodal_matrix(B, proj_vector=proj_orthogonal_vector_q):
    Bp = []
    Bp.append(B[0])
    for i in [1..B.nrows()-1]:
        p = B[i]
        for j in [0..i-1]:
            pp = B[j]
            p = proj_vector(pp, p)
        Bp.append(p)
    return matrix(Bp)

def ell(B, i):
    return epipodal_vector(B, i).hamming_weight()

def ell_profile(B):
    return [epipodal_vector(B, i).hamming_weight() for i in range(B.nrows())]

def get_k1(ells):
    return len([x for x in ells if x != 1])


# Read input parameters
parser = argparse.ArgumentParser(description='This script compares efficiency of different BKZ modes for codes over F_q')
parser.add_argument("-q", required=True, type=int, default=0, help='Finite field cardinality, prime')
parser.add_argument("-b", required=True, type=int, default=0, help='Beta block size in BKZ')
parser.add_argument("-s", required=True, type=int, default=10, help='Number of samples')
args = parser.parse_args()
beta = args.b
q = args.q
samples = args.s

F=GF(q)
n_set = [32,64,128,256,512,1024]
R=0.5

allstat = []
for n in n_set:
    wtH0 = [] # B[0].hamming_wt (proper basis)
    wtH1 = []
    wtH2 = []
    wtH3 = []
    wtH4 = []
    wtL0 = [] # B[0].lee_wt (proper basis)
    wtL1 = []
    wtL2 = []
    wtL3 = []
    wtL4 = []
    T1 = []
    T2 = []
    T3 = []
    T4 = []
    k=int(n*R)
    for s in range(samples):
        while true:
            B = random_matrix(F, k, n)
            if B[:,:k].rank() == k:
                break
        assert B.rank() == B.nrows()
        B = proper_basis(B)
        print(B, flush=true)
        B_profile = ell_profile(B)
        print(f"profile: {B_profile}")
        print(f"{B[0].hamming_weight()}<=B[0].wtL<={B[0].hamming_weight()*floor(F.cardinality()/2)}")
        print(f"B[0].wtH={B[0].hamming_weight()}, B[0].wtL={WtLeeVec(B[0])}")
        wtH0.append(B[0].hamming_weight())
        wtL0.append(WtLeeVec(B[0]))

        print(f"\n--- BKZ with beta={beta}, classic ---")
        t=time.time()
        B_red=bkz(B, beta, [], 0)
        T1.append(time.time()-t)
        B_red_profile = ell_profile(B_red)
        print(f"profile after BKZ: {B_red_profile}")
        print(f"{B_red[0].hamming_weight()}<=B_red[0].wtL<={B_red[0].hamming_weight()*floor(F.cardinality()/2)}")
        print(f"B_red[0].wtH={B_red[0].hamming_weight()}, B_red[0].wtL={WtLeeVec(B_red[0])}")
        wtH1.append(B_red[0].hamming_weight())
        wtL1.append(WtLeeVec(B_red[0]))

        print(f"\n--- BKZ with beta={beta}, selecting lowest lee wt. vec ---")
        t=time.time()
        B_red1=bkz(B, beta, [], 1)
        T2.append(time.time()-t)
        B_red1_profile = ell_profile(B_red1)
        print(f"profile after BKZ: {B_red1_profile}")
        print(f"{B_red1[0].hamming_weight()}<=B_red1[0].wtL<={B_red1[0].hamming_weight()*floor(F.cardinality()/2)}")
        print(f"B_red1[0].wtH={B_red1[0].hamming_weight()}, B_red1[0].wtL={WtLeeVec(B_red1[0])}")
        wtH2.append(B_red1[0].hamming_weight())
        wtL2.append(WtLeeVec(B_red1[0]))

        print(f"\n--- LeeBKZ with beta={beta}, bkz with minimum_lee_wt_codeword+primitive ---")
        t=time.time()
        B_red2=bkz(B, beta, [], 2)
        T3.append(time.time()-t)
        B_red2_profile = ell_profile(B_red2)
        print(f"profile after BKZ: {B_red2_profile}")
        print(f"{B_red2[0].hamming_weight()}<=B_red2[0].wtL<={B_red2[0].hamming_weight()*floor(F.cardinality()/2)}")
        print(f"B_red2[0].wtH={B_red2[0].hamming_weight()}, B_red2[0].wtL={WtLeeVec(B_red2[0])}")
        wtH3.append(B_red2[0].hamming_weight())
        wtL3.append(WtLeeVec(B_red2[0]))

        print(f"\n--- LeeBKZ with beta={beta}, bkz with minimum_lee_wt_codeword (may be non-primitive) ---")
        t=time.time()
        B_red3=bkz(B, beta, [], 3)
        T4.append(time.time()-t)
        B_red3_profile = ell_profile(B_red3)
        print(f"profile after BKZ: {B_red3_profile}")
        print(f"{B_red3[0].hamming_weight()}<=B_red3[0].wtL<={B_red3[0].hamming_weight()*floor(F.cardinality()/2)}")
        print(f"B_red3[0].wtH={B_red3[0].hamming_weight()}, B_red3[0].wtL={WtLeeVec(B_red3[0])}")
        wtH4.append(B_red3[0].hamming_weight())
        wtL4.append(WtLeeVec(B_red3[0]))

    print("\n\n***************************************")
    print(f"******** STATISTICS FOR n={n} ********")
    print("***************************************\n")
    print("Modes:")
    print("(0) Input matrix (proper basis)")
    print("(1) Classic BKZ")
    print("(2) Classic BKZ + selecting lower Lee wt vec")
    print("(3) LEEBKZ + selecting primitive codeword of minimal Lee wt")
    print("(4) LEEBKZ + minimal Lee wt codeword without primitivity (InsertVector)")
    print(f"Reduced {samples} matrices of size {k}x{n}")
    print(f"(0): wtH[0]={numpy.mean(wtH0)} wtL[0]={numpy.mean(wtL0)}")
    print(f"(1): wtH[0]={numpy.mean(wtH1)} wtL[0]={numpy.mean(wtL1)} time={numpy.mean(T1)}")
    print(f"(2): wtH[0]={numpy.mean(wtH2)} wtL[0]={numpy.mean(wtL2)} time={numpy.mean(T2)}")
    print(f"(3): wtH[0]={numpy.mean(wtH3)} wtL[0]={numpy.mean(wtL3)} time={numpy.mean(T3)}")
    print(f"(4): wtH[0]={numpy.mean(wtH4)} wtL[0]={numpy.mean(wtL4)} time={numpy.mean(T4)}")
    #allstat.append([wtH0, wtH1, wtH2, wtH3, wtH4, wtL0, wtL1, wtL2, wtL3, wtL4])
