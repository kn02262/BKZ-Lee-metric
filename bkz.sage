# coding: utf-8

import numpy as np

import sepi
import LEE

def lprint(verbose, *args, **kwargs):
    if verbose:
        print(*args, **kwargs)

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

def or_vector_np(a,b):
    a = np.array(a, dtype=np.uint8)
    b = np.array(b, dtype=np.uint8)
    np.bitwise_or(a, b, out=a)
    return vector(GF(2), a)

def proj_orthogonal_vector(a,b):
    return vector([(a[i]+1)*b[i] for i in range(len(a))])

def proj_orthogonal_vector_np(a,b):
    a = np.array(a, dtype=np.uint8)
    b = np.array(b, dtype=np.uint8)
    np.bitwise_not(a, out=a)
    np.bitwise_and(a, b, out=a)
    return vector(GF(2), a)

def proj_orthogonal_vector_np_pure(a,b):
    a = np.logical_not(a)
    np.logical_and(a, b, out=a)
    return a

def proj_orthogonal_vector_np_bitwise(a,b):
    a = np.bitwise_not(a)
    np.bitwise_and(a, b, out=a)
    return a

def proj_orthogonal_vector_q(a,b):
     # Projection of b orthogonal to a
     res = []
     for i in range(len(a)):
         if a[i] != 0:
             res.append(0)
         else:
             res.append(b[i])
     return vector(res)

def orthogonal_projections(B, i, j, S=None):
    k = B.nrows()
    if B.base_ring() == GF(2):
        #proj_orth = proj_orthogonal_vector
        proj_orth = proj_orthogonal_vector_np
    else:
        proj_orth = proj_orthogonal_vector_q
    assert i < k and i >= 0 and i <= j
    assert j >= 0
    if i == 0:
        return matrix(B.base_ring(), B[i:j+1,:])
    if S == None:
        S = sepi.S_aux(B)
    S.update(B, i)
    prj = []
    for i0 in range(i, min(j+1,k)):
        pr = proj_orth(S[i], B[i0])
        prj.append(pr)
    B = matrix(B.base_ring(), prj)
    return B

def orthogonal_projections_vec(B, i, vec, S=None):
    k = B.nrows()
    if B.base_ring() == GF(2):
        #proj_orth = proj_orthogonal_vector
        proj_orth = proj_orthogonal_vector_np
    else:
        proj_orth = proj_orthogonal_vector_q
    assert i < k and i >= 0
    if i == 0:
        return vec
    if S == None:
        S = sepi.S_aux(B)
    S.update(B, i)
    return proj_orth(S[i], vec)

def insert_primitive(B, p):
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
    B_proj = orthogonal_projections(B, 1, k - 1)
    B_proj_ech = B_proj.echelon_form()

    T = B_proj.solve_left(B_proj_ech)
    assert T.nrows() == k-1
    assert T.ncols() == k-1
    T1 = block_diagonal_matrix(identity_matrix(1), T)
    B = T1 * B
    A = B0.solve_left(B)
    assert (A*B0).rank() == k
    assert (A*B0)[0] == p
    #assert (A*B0)[0].hamming_weight() == p.hamming_weight(), f"{(A*B0)[0].hamming_weight()} != {p.hamming_weight()}"
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

def make_primitive(B, c): # Definition: Claim 4.7 [More basis reduction for linear codes.]
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
    P=Permutation(B1_pivots).to_matrix().change_ring(GF(2))
    B1=(B*P).echelon_form()
    cp = matrix(B1[m])*P^(-1)
    return vector(cp[0])

def MW_Codeword_Hamming_fast(B, lee=0):
    # lee=0: Finds minimum hamming weight codeword
    # lee=1: Finds minimum hamming weight and also minimum lee weight codeword
    zero_pos = []
    Bn=B.ncols()
    for j in range(Bn):
        col = B.column(j)
        if col.is_zero():
            zero_pos.append(j)
    B=B.delete_columns(zero_pos)
    C=LinearCode(B)
    if lee == false:
        p=C[0]
        dw=C.length()
        for c in C:
            if c.hamming_weight() != 0 and c.hamming_weight()<dw:
                dw=c.hamming_weight()
                p=c
    else:
        p=C[0]
        dw=C.length()
        dl=C.length()*floor(B.base_ring().cardinality()/2)
        for c in C:
            if c.hamming_weight() != 0 and c.hamming_weight()<dw:
                dw=c.hamming_weight()
                dl=LEE.WtLeeVec(c)
                p=c
            elif c.hamming_weight != 0 and c.hamming_weight()==dw:
                if LEE.WtLeeVec(c) < dl:
                    dl=LEE.WtLeeVec(c)
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


def MW_Codeword_LEE(B, primitive=False):
    if B.rank() == 0:
        F=B.base_ring()
        p = zero_vector(F, B.ncols())
    else:
        C=LEE.LinearLeeMetricCode(B)
        p = C.minimum_lee_weight_codeword_fast(primitive=primitive)
    return p


def bkz(B, beta, lee=0, verbose=True, STAT=[0,0]):

    # STAT[0] = Number of minimum_weight_codeword() calls;
    # STAT[1] = Number of fails due to potential function growth in LLL

    k = B.nrows()
    assert 2 <= beta and beta <= k
    S_epi = sepi.S_aux(B)
    i = 0
    while i < k - 1:
        lprint(verbose, f"i = {i}", flush=True)
        j = min(i+beta-1, k-1)
        B_proj = orthogonal_projections(B, i, j, S_epi)

        if B_proj[0].hamming_weight() == 1:
            i = detect_tail(B, i, verbose=verbose)
            continue

        p = MW_Codeword_Hamming_fast(B_proj, lee)
        STAT[0] = STAT[0]+1
        if B_proj[0].hamming_weight() == p.hamming_weight():
            i = i + 1
        else:
            lprint(verbose, f"i = {i}, new weight: {p.hamming_weight()}", flush=True)
            if i == 0:
                lprint(verbose, f"codeword = {p}", flush=True)
            A = insert_primitive(B_proj, p)
            T = block_diagonal_matrix(identity_matrix(i), A, identity_matrix(k-j-1))
            B = T * B
            S_epi.set_pos(i)
            if i == 0:
                i += 1 # prevent double calculation for first block
            else:
                i = max(0, i-beta+1)
    return B

def LEELLL(B, beta, verbose=True, STAT=[0,0]):

    # STAT[0] = Number of minimum_weight_codeword() calls;
    # STAT[1] = Number of fails due to potential function growth in LLL

    k = B.nrows()
    assert 2 <= beta and beta <= k
    S_epi = sepi.S_aux(B)
    i = 0
    while i < k - 1:
        lprint(verbose, f"i = {i}", flush=True)
        j = min(i+beta-1, k-1)
        B_proj = orthogonal_projections(B, i, j, S_epi)

        if LEE.WtLeeVec(B_proj[0]) == 1:
            i = detect_tail_LEE(B, i, verbose=verbose)
            continue

        p = MW_Codeword_LEE(B_proj, primitive=False)
        STAT[0] = STAT[0]+1
        B_epi = epipodal_matrix(B, partial=i+2)
        li = LEE.WtLeeVec(B_epi[i])
        li1 = LEE.WtLeeVec(B_epi[i+1])
        if (LEE.WtLeeVec(B_proj[0]) > LEE.WtLeeVec(p) or LEE.WtLeeVec(B_proj[0]) == 0) and (LEE.WtLeeVec(p) > 0):
            lprint(verbose, f"i = {i}, new Lee weight: {LEE.WtLeeVec(p)}", flush=True)
            if i == 0:
                lprint(verbose, f"codeword = {p}", flush=True)
            A = insert_vector(B_proj, p)
            T = block_diagonal_matrix(identity_matrix(i), A, identity_matrix(k-j-1))
            E = T * B
            E_epi = epipodal_matrix(E)
            lip1 = LEE.WtLeeVec(E_epi[i+1])
            if (li-LEE.WtLeeVec(E_epi[i])+li1-lip1 < 1): # lip1 = LEE.WtLeeVec(epipodal_vector(E, i+1))
                lprint(verbose, "LLL Potential growth!")
                STAT[1] = STAT[1] + 1
                i=i+1
                continue
            B = E
            S_epi.set_pos(i)
            if i == 0:
                i += 1 # prevent the first block to be double reduced
            else:
                i = max(0, i-beta+1)
        else:
            i += 1
    return B


def LEEbkz(B, beta, primitive=False, verbose=True, STAT=[0,0]):

    # STAT[0] = Number of minimum_weight_codeword() calls;
    # STAT[1] = Number of fails due to potential function growth in LLL

    k = B.nrows()
    assert 2 <= beta and beta <= k
    S_epi = sepi.S_aux(B)
    i = 0
    while i < k - 1:
        lprint(verbose, f"i = {i}", flush=True)
        j = min(i+beta-1, k-1)
        B_proj = orthogonal_projections(B, i, j, S_epi)

        if LEE.WtLeeVec(B_proj[0]) == 1:
            i = detect_tail_LEE(B, i, verbose=verbose)
            continue

        p = MW_Codeword_LEE(B_proj, primitive=primitive)
        STAT[0] = STAT[0]+1
        if (LEE.WtLeeVec(B_proj[0]) > LEE.WtLeeVec(p)) and (LEE.WtLeeVec(p) > 0):
            lprint(verbose, f"i = {i}, new Lee weight: {LEE.WtLeeVec(p)}", flush=True)
            if i == 0:
                lprint(verbose, f"codeword = {p}", flush=True)
            if primitive==True:
                A = insert_primitive(B_proj, p)
            else:
                A = insert_vector(B_proj, p)
            T = block_diagonal_matrix(identity_matrix(i), A, identity_matrix(k-j-1))
            B = T * B
            S_epi.set_pos(i)
            B_proj = orthogonal_projections(B, i+1, k-1, S_epi) # Matrix of size (k-1-i, n)
            s=0
            while((LEE.WtLeeVec(B_proj[s]) == 0) and (s<k-2-i)):
                s=s+1
            if (s>0):
                B.swap_rows(i+1, i+1+s)
            if i == 0:
                i += 1 # prevent the first block to be double reduced
            else:
                i = max(0, i-beta+1)
        else:
            i += 1
    return B

def detect_tail(B, i, verbose=False):
    k = B.nrows()
    t = walltime()
    B_epi = epipodal_matrix(B)
    lprint(verbose, f"[detect_tail:epipodal_matrix:{walltime()-t} sec.]")
    while i <= k - 1 and B_epi[i].hamming_weight() == 1:
        i += 1
    return i

def detect_tail_LEE(B, i, verbose=False):
    k = B.nrows()
    t = walltime()
    B_epi = epipodal_matrix(B)
    lprint(verbose, f"[detect_tail:epipodal_matrix:{walltime()-t} sec.]")
    while i <= k - 1 and LEE.WtLeeVec(B_epi[i]) <= 1:
        i += 1
    return i

def lll(B):
    return bkz(B, 2)

def lll_v2(B):
    return bkz_v2(B, 2)

def epipodal_vector(B, i):
    # Build a proj system
    p = vector([0 for x in range(B.ncols())])
    for i0 in range(i):
        p = or_vector(p, B[i0])
    return proj_orthogonal_vector_q(p,B[i])

def epipodal_matrix(B, partial=None):
    if B.base_ring() == GF(2):
        return epipodal_matrix_np_bitwise(B, partial)
    Bp = []
    p = vector([0 for x in range(B.ncols())])
    if partial == None:
        s = B.nrows()
    else:
        s = min(B.nrows(), partial)
    for i in range(s):
        Bp.append(proj_orthogonal_vector_q(p, B[i]))
        p = or_vector(p, B[i])
    return matrix(Bp)

def epipodal_matrix_np(B):
    Bp = []
    p = np.zeros(B.ncols(), dtype=bool)
    for i in range(B.nrows()):
        bi_np = np.array(B[i], dtype=bool)
        Bp.append(proj_orthogonal_vector_np_pure(p, bi_np))
        np.logical_or(p, bi_np, out=p)
    return matrix(GF(2), Bp)

def epipodal_matrix_np_bitwise(B, partial=None):
    Bp = []
    p = np.zeros(B.ncols(), dtype=bool)
    p = np.packbits(p)
    if partial == None:
        n0 = B.nrows()
    else:
        n0 = partial
    for i in range(n0):
        bi_np = np.array(B[i], dtype=bool)
        bi_np = np.packbits(bi_np)
        Bp.append(proj_orthogonal_vector_np_bitwise(p, bi_np))
        np.bitwise_or(p, bi_np, out=p)
    return matrix(GF(2), [np.unpackbits(Bp[i]) for i in range(len(Bp))])

def ell(B, i):
    return epipodal_vector(B, i).hamming_weight()

def ell_profile(B, lee=False):
    if B.base_ring() == GF(2):
        B_epi = epipodal_matrix_np_bitwise(B)
        #t = walltime()
        #B_epi = epipodal_matrix(B)
        #print(f"[ell_profile:epipodal_matrix:{walltime() - t}]")
        #t = walltime()
        #B_epi = epipodal_matrix_np(B)
        #print(f"[ell_profile:epipodal_matrix_np:{walltime() - t}]")
    else:
        B_epi = epipodal_matrix(B)
    if lee==True:
        return [LEE.WtLeeVec(B_epi[i]) for i in range(B.nrows())]
    else:
        return [B_epi[i].hamming_weight() for i in range(B.nrows())]
    #return [epipodal_vector(B, i).hamming_weight() for i in range(B.nrows())]

def get_k1(ells):
    return len([x for x in ells if x != 1])

def epi_sort_bin(B):
    k = B.nrows()
    n = B.ncols()
    proj_orth = proj_orthogonal_vector_np_bitwise
    p = np.zeros(n, dtype=bool)
    p = np.packbits(p)
    B_np = []
    for i in range(k):
        bi_np = np.array(B[i], dtype=bool)
        bi_np = np.packbits(bi_np)
        B_np.append(bi_np)
    for i in range(k):
        best_j = i
        best_w = np.bitwise_count(proj_orth(p, B_np[i])).sum()
        for j in range(i+1,k):
            w = np.bitwise_count(proj_orth(p, B_np[j])).sum()
            if w < best_w:
                best_w = w
                best_j = j
        if i != best_j:
            # swap
            B_np[i], B_np[best_j] = B_np[best_j], B_np[i]
        np.bitwise_or(p, B_np[i], out=p)
    return matrix(GF(2), [np.unpackbits(B_np[i]) for i in range(k)])

def epi_sort(B):
    k = B.nrows()
    n = B.ncols()
    if B.base_ring() == GF(2):
        B[:] = epi_sort_bin(B)
        #proj_orth = proj_orthogonal_vector
        #proj_orth = proj_orthogonal_vector_np
        return
    else:
        proj_orth = proj_orthogonal_vector_q
    p = vector([0 for x in range(n)])
    for i in range(k):
        best_j = i
        best_w = proj_orth(p, B[i]).hamming_weight()
        for j in range(i+1,k):
            w = proj_orth(p, B[j]).hamming_weight()
            if w < best_w:
                best_w = w
                best_j = j
        if i != best_j:
            B.swap_rows(i, best_j)
        if B.base_ring() == GF(2):
            p = or_vector_np(p, B[i])
        else:
            p = or_vector(p, B[i])
    return

def subfields_basis(B):
    F = B.base_ring()
    sf = sorted(F.subfields(), key=lambda x: x[0].degree())
    C = LinearCode(B)
    S = matrix(B.base_ring(), 0, B.ncols())
    found = False
    for i in range(len(sf)):
        s = sf[i][0]
        C_s = codes.SubfieldSubcode(C, s, sf[i][1])
        G_s = C_s.generator_matrix()
        for j in range(G_s.nrows()):
            S0 = S.stack(G_s[j])
            if S0.rank() == S.rank() + 1:
                S = S0
            if S.rank() == B.nrows():
                found = True
                break
        if found:
            break
    assert S.nrows() == B.nrows()
    A = B.solve_left(S)
    assert S == A*B
    assert not A.is_singular()
    B = S
    return B

def randomize(B):
    while True:
        G = random_matrix(B.base_ring(), B.nrows(), B.nrows())
        if not G.is_singular():
            break
    return G,G*B

def preprocess(B, algs=["systemize", "episort"]):
    for i in range(len(algs)):
        if algs[i] == 'systemize':
            B = proper_basis(B)
        elif algs[i] == 'randomize':
            _,B = randomize(B)
        elif algs[i] == 'episort':
            B = copy(B)
            epi_sort(B)
        elif algs[i] == 'subfields':
            B = subfields_basis(B)
    return B
