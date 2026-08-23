#include <immintrin.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define ALIGN 32

// Optimized AVX2 Lee Weight evaluation
uint32_t compute_lee_weight_simd(const int16_t *codeword, size_t n, int16_t q) {
    __m256i vq = _mm256_set1_epi16(q);
    __m256i vzero = _mm256_setzero_si256();
    __m256i vsum = _mm256_setzero_si256();
    
    size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m256i x = _mm256_load_si256((const __m256i *)(codeword + i));
        __m256i q_minus_x = _mm256_sub_epi16(vq, x);
        __m256i lee_val = _mm256_min_epi16(x, q_minus_x);
        
        __m256i lo = _mm256_unpacklo_epi16(lee_val, vzero);
        __m256i hi = _mm256_unpackhi_epi16(lee_val, vzero);
        vsum = _mm256_add_epi32(vsum, _mm256_add_epi32(lo, hi));
    }
    
    uint32_t total_weight = 0;
    alignas(ALIGN) int32_t buffer[8];
    _mm256_store_si256((__m256i *)buffer, vsum);
    for (int j = 0; j < 8; j++) {
        total_weight += buffer[j];
    }
    
    for (; i < n; i++) {
        int16_t x = codeword[i];
        int16_t diff = q - x;
        total_weight += (x < diff) ? x : diff;
    }
    
    return total_weight;
}

void vec_add_mod_simd(int16_t *res, const int16_t *a, const int16_t *b, size_t n, int16_t q) {
    __m256i vq = _mm256_set1_epi16(q);
    size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m256i va = _mm256_load_si256((const __m256i *)(a + i));
        __m256i vb = _mm256_load_si256((const __m256i *)(b + i));
        __m256i vsum = _mm256_add_epi16(va, vb);
        
        __m256i mask = _mm256_cmpgt_epi16(vsum, _mm256_sub_epi16(vq, _mm256_set1_epi16(1)));
        __m256i vsub = _mm256_and_si256(mask, vq);
        
        _mm256_store_si256((__m256i *)(res + i), _mm256_sub_epi16(vsum, vsub));
    }
    for (; i < n; i++) {
        int32_t sum = a[i] + b[i];
        res[i] = (sum >= q) ? (sum - q) : sum;
    }
}

void vec_sub_mod_simd(int16_t *res, const int16_t *a, const int16_t *b, size_t n, int16_t q) {
    __m256i vq = _mm256_set1_epi16(q);
    __m256i vzero = _mm256_setzero_si256();
    size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m256i va = _mm256_load_si256((const __m256i *)(a + i));
        __m256i vb = _mm256_load_si256((const __m256i *)(b + i));
        __m256i vdiff = _mm256_sub_epi16(va, vb);
        
        // Correct fix for missing _mm256_cmplt_epi16 instruction
        __m256i mask = _mm256_cmpgt_epi16(vzero, vdiff);
        __m256i vadd = _mm256_and_si256(mask, vq);
        
        _mm256_store_si256((__m256i *)(res + i), _mm256_add_epi16(vdiff, vadd));
    }
    for (; i < n; i++) {
        int32_t diff = a[i] - b[i];
        res[i] = (diff < 0) ? (diff + q) : diff;
    }
}

// C Interface directly mapping to NumPy inputs
uint32_t find_min_weight_codeword(const int16_t *G, size_t k, size_t n, int16_t q, int16_t *best_codeword) {
    size_t padded_n = ((n + 15) / 16) * 16;
    
    int16_t *current_codeword = (int16_t *)_mm_malloc(padded_n * sizeof(int16_t), ALIGN);
    int16_t *aligned_G = (int16_t *)_mm_malloc(k * padded_n * sizeof(int16_t), ALIGN);
    
    memset(current_codeword, 0, padded_n * sizeof(int16_t));
    memset(aligned_G, 0, k * padded_n * sizeof(int16_t));
    
    for (size_t r = 0; r < k; r++) {
        memcpy(aligned_G + r * padded_n, G + r * n, n * sizeof(int16_t));
    }

    uint32_t min_weight = UINT32_MAX;
    
    uint64_t total_codewords = 1;
    for (size_t i = 0; i < k; i++) total_codewords *= q;

    int *u = (int *)calloc(k, sizeof(int));
    int *direction = (int *)malloc(k * sizeof(int));
    for (size_t i = 0; i < k; i++) direction[i] = 1;

    for (uint64_t step = 1; step < total_codewords; step++) {
        size_t change_idx = 0;
        uint64_t temp = step;
        while (temp % q == 0) {
            change_idx++;
            temp /= q;
        }

        const int16_t *row = aligned_G + change_idx * padded_n;
        
        if (direction[change_idx] == 1) {
            vec_add_mod_simd(current_codeword, current_codeword, row, padded_n, q);
            u[change_idx]++;
            if (u[change_idx] == q - 1) direction[change_idx] = -1;
        } else {
            vec_sub_mod_simd(current_codeword, current_codeword, row, padded_n, q);
            u[change_idx]--;
            if (u[change_idx] == 0) direction[change_idx] = 1;
        }

        uint32_t weight = compute_lee_weight_simd(current_codeword, padded_n, q);
        if (weight < min_weight) {
            min_weight = weight;
            memcpy(best_codeword, current_codeword, n * sizeof(int16_t));
        }
    }

    _mm_free(current_codeword);
    _mm_free(aligned_G);
    free(u);
    free(direction);
    
    return min_weight;
}
