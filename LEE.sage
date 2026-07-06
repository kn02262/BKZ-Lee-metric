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
