from sage.rings.infinity import Infinity

from sage.coding.linear_code_no_metric import AbstractLinearCodeNoMetric
from sage.coding.linear_code import LinearCodeGeneratorMatrixEncoder
from sage.coding.decoder import Decoder

class AbstractLinearLeeMetricCode(AbstractLinearCodeNoMetric):
    def __init__(self, base_field, length, default_encoder_name, default_decoder_name, basis=None):
        self._generic_constructor = LinearLeeMetricCode
        super().__init__(base_field, length, default_encoder_name, default_decoder_name, "lee")
        
        
class LinearLeeMetricCode(AbstractLinearLeeMetricCode):
    def __init__(self, generator, basis=None):
        base_field = generator.base_ring()
        if not base_field.is_field():
            raise ValueError("'generator' must be defined on a field (not a ring)")

        try:
            gen_basis = None
            if hasattr(generator, "nrows"):  # generator matrix case
                if generator.rank() < generator.nrows():
                    gen_basis = generator.row_space().basis()
            else:
                gen_basis = generator.basis()  # vector space etc. case
            if gen_basis is not None:
                from sage.matrix.constructor import matrix
                generator = matrix(base_field, gen_basis)
                if generator.nrows() == 0:
                    raise ValueError("this linear code contains no nonzero vector")
        except AttributeError:
            # Assume input is an AbstractLinearLeeMetricCode, extract its generator matrix
            generator = generator.generator_matrix()

        self._generator_matrix = generator
        self._length = generator.ncols()
        super().__init__(base_field, self._length, "GeneratorMatrix", "NearestNeighbor", basis)


    def lee_weight(self, word):

        r"""
        Return the Lee weight of the word.

        INPUT:

        - ``word`` -- a vector over the ``base_field`` of ``self``

        """
        q = self.base_field().cardinality()
        n = self.length()
        return sum( int(min(word[i], q-word[i])) for i in range(n) )
        
    def lee_distance(self, w1, w2):

        r"""
        Return the Lee distance between two word w1, w2.

        INPUT:

        - ``w1``, ``w2`` -- vectors over the ``base_field`` of ``self``

        """
        q = self.base_field().cardinality()
        n = self.length()
        return self.lee_weight( w1-w2 )
        
        
    def minimum_lee_distance(self):

        r"""
        Return the minimum distance of ``self``.

        This algorithm simply iterates over all the elements of the code ``self`` and
        returns the minimum weight.
            
        """
        d = Infinity
        for c in self:
            if c == self.zero():
                continue
            d = min(self.lee_weight(c), d)
        return d


    def IsPrimitiveCodeword(self, c): # Definition: Claim 4.7 [GD'24 "More basis reduction for linear codes: backward reduction, BKZ, slide reduction and more"]

        r"""
        Checks, whether the codeword ``c`` is minimal in ``self``,
        i.e. there is no nonzero codeword e, Supp(e) \subset Supp(c)

        """

        B=self.generator_matrix()
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
        return (c == vector(cp)) or (not set(vector(cp).support()).issubset(S))


    def minimum_lee_weight_codeword(self, minHw=false, primitive=false):
        r"""
        Return the codeword of minimum Lee weight of ``self``.
        If minHw==True, additionally selects the codeword with possibly less Hamming weight
        If primitive=True, returns only primitive codeword

        This algorithm simply iterates over all the elements of the code and
        returns the codeword of minimum weight.
            
        """
        d = Infinity
        w = Infinity
        lc = self.zero()
        for c in self:
            if c == self.zero():
                continue
            if self.lee_weight(c) < d:
                if primitive and not self.IsPrimitiveCodeword(c):
                    continue
                lc = c
                d = self.lee_weight(c)
                w = c.hamming_weight()
            elif self.lee_weight(c) == d and c.hamming_weight() < w and minHw:
                if primitive and not self.IsPrimitiveCodeword(c):
                    continue
                w = c.hamming_weight()
                lc = c
        return lc
        
        
    def minimum_lee_weight_codeword_fast(self, minHw=false, primitive=false):
        r"""
        Return the codeword of minimum Lee weight of ``self``.
        If minHw==True, additionally selects the codeword with possibly less Hamming weight
        If primitive=True, returns only primitive codeword

        This algorithm simply iterates over all the elements of the code and
        returns the codeword of minimum weight.
            
        """
        B=self.generator_matrix()
        zero_pos = []
        Bn=B.ncols()
        for j in range(Bn):
            col = B.column(j)
            if col.is_zero():
                zero_pos.append(j)
        B=B.delete_columns(zero_pos)
        C = LinearLeeMetricCode(B)
        d = Infinity
        w = Infinity
        lc = C.zero()
        for c in C:
            if c == C.zero():
                continue
            if C.lee_weight(c) < d:
                if primitive and not C.IsPrimitiveCodeword(c):
                    continue
                lc = c
                d = C.lee_weight(c)
                w = c.hamming_weight()
            elif C.lee_weight(c) == d and c.hamming_weight() < w and minHw:
                if primitive and not C.IsPrimitiveCodeword(c):
                    continue
                w = c.hamming_weight()
                lc = c
        res = []
        ptr = 0
        for i in range(Bn):
            if i in zero_pos:
                res.append(0)
            else:
                res.append(lc[ptr])
                ptr = ptr + 1
        return vector(res)
        
        
    def _repr_(self):
        R = self.base_field()
        if R in Fields():
            return "[%s, %s] linear Lee metric code over GF(%s)" % (self.length(), self.dimension(), R.cardinality())
        else:
            return "[%s, %s] linear Lee metric code over %s" % (self.length(), self.dimension(), R)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.
        """
        return "[%s, %s]\\textnormal{ Linear Lee metric code over }%s"\
                % (self.length(), self.dimension(), self.base_field()._latex_())

    def generator_matrix(self, encoder_name=None, **kwargs):
        r"""
        Return a generator matrix of ``self``.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          used to compute the generator matrix. ``self._generator_matrix``
          will be returned if default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used

        """
        if encoder_name is None or encoder_name == 'GeneratorMatrix':
            g = self._generator_matrix
        else:
            g = super().generator_matrix(encoder_name, **kwargs)
        g.set_immutable()
        return g




####################### ToDo: Implement Decoder ###############################


class LinearLeeMetricCodeNearestNeighborDecoder(Decoder):

    def __init__(self, code):
        super().__init__(code, code.ambient_space(), code._default_encoder_name)

    def __eq__(self, other):
        r"""
        Test equality between LinearLeeMetricCodeNearestNeighborDecoder objects.

        """
        return isinstance(other, LinearLeeMetricCodeNearestNeighborDecoder)\
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        """
        return "Nearest neighbor decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        """
        return "\\textnormal{Nearest neighbor decoder for }%s" % self.code()._latex_()

    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``word`` and returns a codeword.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT: a vector of ``self``'s message space

        """
        C = self.code()
        c_min = C.zero()

        return c_min

    def decoding_radius(self):
        r"""
        Return maximal number of errors ``self`` can decode.

        """
        return (self.code().minimum_distance() - 1) // 2
    
####################### registration ###############################

LinearLeeMetricCode._registered_encoders["GeneratorMatrix"] = LinearCodeGeneratorMatrixEncoder
LinearLeeMetricCode._registered_decoders["NearestNeighbor"] = LinearLeeMetricCodeNearestNeighborDecoder

def WtLee(x): # Coordinate Lee weight
    F=x.base_ring()
    if(x == 0):
        return 0
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
    F=B.base_ring()
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

def bkz(B, beta, preproc=[], lee=0, ignore_potential_growth=false, STAT=[0,0]):
    # lee=0: hamming bkz
    # lee=1: hamming bkz + lower WtLee codeword selection
    # lee=2: lee bkz (inserting primitive codeword)
    # lee=3: lee bkz (inserting codeword may be non-primitive)

    # proj_vector=proj_orthogonal_vector_q: Projection with non-intersecting supports
    # proj_vector=proj_orthogonal_lee_vector_q: Projection with orthogonal elements
    proj_vector=proj_orthogonal_vector_q

    # STAT[0] = Number of minimum_weight_codeword() calls;
    # STAT[1] = Number of fails due to potential function growth in LLL

    k = B.nrows()
    assert 2 <= beta and beta <= k
    i = 0
    while i < k - 1:
        j = min(i+beta-1, k-1)
        B_proj = orthogonal_projections(B, i, j, proj_vector)
        STAT[0] = STAT[0] + 1
        if WtLeeVec(B_proj[0]) == 1:
            i=i+1
            continue
        p = MW_Codeword(B_proj, lee)
        li = WtLeeVec(epipodal_vector(B, i))
        li1 = WtLeeVec(epipodal_vector(B, i+1))
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

            # Try to minimize l_{i+1}' - l_{i+1} (for LLL) <=> minimize l_{i+1}'
            if (lee == 3) and (beta == 2) and (not ignore_potential_growth):
                lip1 = WtLeeVec(epipodal_vector(E, i+1))
                if (li-WtLeeVec(epipodal_vector(E, i))+li1-lip1 < 1): # lip1 = WtLeeVec(epipodal_vector(E, i+1))
                    print("LLL Potential growth!")
                    STAT[1] = STAT[1] + 1
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
