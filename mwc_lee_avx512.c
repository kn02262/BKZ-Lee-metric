#include <immintrin.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// Для AVX-512 выравнивание должно быть 64 байта
#define ALIGN 64

// Оптимизированный подсчет веса Ли с использованием AVX-512
uint32_t compute_lee_weight_simd(const int16_t *codeword, size_t n, int16_t q) {
    __m512i vq = _mm512_set1_epi16(q);
    __m512i vsum = _mm512_setzero_si512(); // ИСПРАВЛЕНО c si256 на si512
    
    size_t i = 0;
    for (; i + 32 <= n; i += 32) {
        __m512i x = _mm512_load_si512((const __m512i *)(codeword + i));
        __m512i q_minus_x = _mm512_sub_epi16(vq, x);
        __m512i lee_val = _mm512_min_epi16(x, q_minus_x);
        
        __m512i lo = _mm512_cvtepi16_epi32(_mm512_castsi512_si256(lee_val));
        __m512i hi = _mm512_cvtepi16_epi32(_mm512_extracti64x4_epi64(lee_val, 1));
        
        vsum = _mm512_add_epi32(vsum, _mm512_add_epi32(lo, hi));
    }
    
    uint32_t total_weight = _mm512_reduce_add_epi32(vsum);
    
    for (; i < n; i++) {
        int16_t x = codeword[i];
        int16_t diff = q - x;
        total_weight += (x < diff) ? x : diff;
    }
    
    return total_weight;
}

void vec_add_mod_simd(int16_t *res, const int16_t *a, const int16_t *b, size_t n, int16_t q) {
    __m512i vq = _mm512_set1_epi16(q);
    __m512i vthresh = _mm512_sub_epi16(vq, _mm512_set1_epi16(1));
    
    size_t i = 0;
    for (; i + 32 <= n; i += 32) {
        __m512i va = _mm512_load_si512((const __m512i *)(a + i));
        __m512i vb = _mm512_load_si512((const __m512i *)(b + i));
        __m512i vsum = _mm512_add_epi16(va, vb);
        
        // В AVX-512 сравнение возвращает маску инструкций __mmask32 вместо регистра-маски
        __mmask32 mask = _mm512_cmpgt_epi16_mask(vsum, vthresh);
        
        // Вычитаем q только там, где взведены биты в маске (используем маскированное вычитание)
        __m512i vres = _mm512_mask_sub_epi16(vsum, mask, vsum, vq);
        
        _mm512_store_si512((__m512i *)(res + i), vres);
    }
    for (; i < n; i++) {
        int32_t sum = a[i] + b[i];
        res[i] = (sum >= q) ? (sum - q) : sum;
    }
}

void vec_sub_mod_simd(int16_t *res, const int16_t *a, const int16_t *b, size_t n, int16_t q) {
    __m512i vq = _mm512_set1_epi16(q);
    __m512i vzero = _mm512_setzero_si512(); // ИСПРАВЛЕНО c si256 на si512
    
    size_t i = 0;
    for (; i + 32 <= n; i += 32) {
        __m512i va = _mm512_load_si512((const __m512i *)(a + i));
        __m512i vb = _mm512_load_si512((const __m512i *)(b + i));
        __m512i vdiff = _mm512_sub_epi16(va, vb);
        
        __mmask32 mask = _mm512_cmplt_epi16_mask(vdiff, vzero);
        __m512i vres = _mm512_mask_add_epi16(vdiff, mask, vdiff, vq);
        
        _mm512_store_si512((__m512i *)(res + i), vres);
    }
    for (; i < n; i++) {
        int32_t diff = a[i] - b[i];
        res[i] = (diff < 0) ? (diff + q) : diff;
    }
}

// C-интерфейс для вызова из NumPy
uint32_t find_min_weight_codeword(const int16_t *G, size_t k, size_t n, int16_t q, int16_t *best_codeword) {
    // Выравнивание под 32 элемента (64 байта) для AVX-512
    size_t padded_n = ((n + 31) / 32) * 32;
    
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
