import ctypes
import os
import sys
import numpy as np

# Bind and find the path configurations for compiling/loading the .so/.dll shared library
_lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'mwc_lee_avx512.so'))
if not os.path.exists(_lib_path):
    _lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'mwc_lee_avx2.so'))

#_lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'mwc_lee_avx2.so'))
#_lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'mwc_lee_avx512.so'))

# Load ctypes handle mappings
_lib = ctypes.CDLL(_lib_path)

# Define native function parameters types matching signature definitions
_lib.find_min_weight_codeword.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.int16, ndim=2, flags='C_CONTIGUOUS'), # G matrix pointer
    ctypes.c_size_t,                                                    # k (rows)
    ctypes.c_size_t,                                                    # n (columns)
    ctypes.c_int16,                                                     # q
    np.ctypeslib.ndpointer(dtype=np.int16, ndim=1, flags='C_CONTIGUOUS') # output codeword pointer
]
_lib.find_min_weight_codeword.restype = ctypes.c_uint32                  # returns min weight value

def find_min_weight(G, q):
    """
    Finds the minimum weight non-zero codeword using AVX2 acceleration loops.
    
    :param G: numpy array of shape (k, n) and dtype int16
    :param q: Integer modulus boundary conditions
    :return: A tuple containing (min_weight, best_codeword_array)
    """
    # Force array structure formatting constraints safely
    G = np.asarray(G, dtype=np.int16, order='C')
    if G.ndim != 2:
        raise ValueError("Generating matrix G must be a 2D array matrix.")
        
    k, n = G.shape
    
    # Initialize empty destination buffer array matching type alignment bounds
    best_codeword = np.zeros(n, dtype=np.int16, order='C')
    
    # Dispatch extraction query step straight to C pointers
    min_weight = _lib.find_min_weight_codeword(G, k, n, int(q), best_codeword)
    
    return min_weight, best_codeword


if __name__ == '__main__':
    # Test Verification Profile
    # Generator matrix G for k=2, n=4 over Z_5
    G_matrix = np.array([ [1, 0, 3, 1],
        [0, 1, 4, 2]
    ], dtype=np.int16)
    
    q_modulus = 5
    
    weight, codeword = find_min_weight(G_matrix, q_modulus)
    
    print("--- Execution Results ---")
    print(f"Minimum Lee Weight: {weight}")
    print(f"Minimal Codeword:   {codeword}")
