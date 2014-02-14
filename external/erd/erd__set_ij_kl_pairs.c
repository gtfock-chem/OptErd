#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "boys.h"
#include "erd.h"
#define ERD_TABLE_FREE_BOYS_FUNCTIONS

#define TOL 1e-14
#ifdef __x86_64__
#include <x86intrin.h>
#elif defined(__MIC__)
#include <immintrin.h>
#endif

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

static YEP_INLINE double pow3o4(double x) {
    return __builtin_sqrt(x * __builtin_sqrt(x));
}

static YEP_INLINE double square(double x) {
    return x * x;
}

static YEP_INLINE double pow4(double x) {
    return square(square(x));
}

static YEP_INLINE double vector_min(const double *YEP_RESTRICT vector, size_t length) {
    double result = vector[0];
    for (size_t i = 1; i < length; i++) {
        const double element = vector[i];
        result = (element < result) ? element : result;
    }
    return result;
}

static const double c0 = 0x1.0B1A240FD5AF4p-8;
static const double c1 = 0x1.352866F31ED93p+0;
static const double c2 = -0x1.567450B98A180p-1;
static const double c3 = 0x1.9721DD0ADF393p-3;
static const double c4 = -0x1.EF4155EFB6D81p-6;
static const double c5 = 0x1.DFF4D0D064A26p-10;
static const double x0 = 0x1.628C5E7D820BFp-5;

#ifdef __AVX__
    YEP_ALIGN(16) static const uint8_t xmm_pack_table[16*16] = {
        0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           4,    5,    6,    7, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    4,    5,    6,    7, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           8,    9,   10,   11, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    8,    9,   10,   11, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           4,    5,    6,    7,    8,    9,   10,   11, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    4,    5,    6,    7,    8,    9,   10,   11, 0x80, 0x80, 0x80, 0x80,
          12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           4,    5,    6,    7,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    4,    5,    6,    7,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80,
           8,    9,   10,   11,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    8,    9,   10,   11,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80,
           4,    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,   15, 0x80, 0x80, 0x80, 0x80,
           0,    1,    2,    3,    4,    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,   15
    };

    static const __m256d ymm_c0 = {  0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8,   0x1.0B1A240FD5AF4p-8  };
    static const __m256d ymm_c1 = {  0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0,   0x1.352866F31ED93p+0  };
    static const __m256d ymm_c2 = { -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1,  -0x1.567450B98A180p-1  };
    static const __m256d ymm_c3 = {  0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3,   0x1.9721DD0ADF393p-3  };
    static const __m256d ymm_c4 = { -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6,  -0x1.EF4155EFB6D81p-6  };
    static const __m256d ymm_c5 = {  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10,  0x1.DFF4D0D064A26p-10 };
    static const __m256d ymm_x0 = {  0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5,   0x1.628C5E7D820BFp-5  };
    static const __m256d ymm_one = { 0x1.0000000000000p+0,   0x1.0000000000000p+0,   0x1.0000000000000p+0,   0x1.0000000000000p+0  };
    static const __m128i xmm_increment_j = { 0x0000000400000004ull, 0x0000000400000004ull };

    static const __m256d ymm_merge_mask_table[4] = {
        { +0.0, +0.0, +0.0, +0.0 },
        { -0.0, -0.0, -0.0, +0.0 },
        { -0.0, -0.0, +0.0, +0.0 },
        { -0.0, +0.0, +0.0, +0.0 }
    };
    
#endif

static YEP_NOINLINE void set_pairs(
    uint32_t npgtoa, uint32_t npgtob, double rnabsq,
    double *YEP_RESTRICT alphaa, double *YEP_RESTRICT alphab,
    uint32_t *YEP_RESTRICT nij_ptr,
    int *YEP_RESTRICT prima, int *YEP_RESTRICT primb,
    double *YEP_RESTRICT rho,
    double qmin, double smaxcd, double rminsq)
{
    const double csmaxcd = smaxcd * (0x1.C5BF891B4EF6Bp-1 / TOL);
    uint32_t nij = 0;
    if (npgtoa > npgtob) {
        double *YEP_RESTRICT alphaa_copy = alphaa;
        alphaa = alphab;
        alphab = alphaa_copy;
        int *YEP_RESTRICT prima_copy = prima;
        prima = primb;
        primb = prima_copy;
        uint32_t npgtoa_copy = npgtoa;
        npgtoa = npgtob;
        npgtob = npgtoa_copy;
    }
#if defined(__AVX__)
    if (npgtob < 4) {
        const double minus_rnabsq = -rnabsq;
        const double rminsq_qmin = rminsq * qmin;
        for (uint32_t i = 0; i < npgtoa; i += 1) {
            const double a = alphaa[i];
            for (uint32_t j = 0; j < npgtob; j += 1) {
                const double b = alphab[j];
                const double p = a + b;
                const double ab = a * b;
                const double pqp = p + qmin;
                const double pqpinv = 1.0 / pqp;
                const double rhoab = __builtin_exp(ab * minus_rnabsq / p);
                const double t = rminsq_qmin * (p / pqp);

                /* pi/4 * f0(x) == x */
                const double x = (t == 0.0) ? x0 : t;
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const double f0 = c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * c5))));
                const double f = f0 < 1.0 ? f0 : 1.0;
                if (ab*square(square(rhoab*csmaxcd) * (ab*f)) >= square((x*pqp)*square(p))) {
                    rho[nij] = rhoab;
                    prima[nij] = i;
                    primb[nij] = j;
                    nij += 1;
                }
            }
        }
    } else {
        const __m256d ymm_minus_rnabsq = _mm256_set1_pd(-rnabsq);
        const __m256d ymm_csmaxcd = _mm256_set1_pd(csmaxcd);
        const __m256d ymm_qmin = _mm256_set1_pd(qmin);
        const __m256d ymm_rminsq_qmin = _mm256_mul_pd(ymm_qmin, _mm256_set1_pd(rminsq));
        const __m128i xmm_npgtob = _mm_set1_epi32(npgtob);
        uint32_t i = 0, j = 0;
        __m128i xmm_j = _mm_set_epi32(3, 2, 1, 0);
        do {
            __m256d ymm_a = _mm256_broadcast_sd(&alphaa[i]);
            __m128i xmm_i = _mm_set1_epi32(i);
            /* Process in blocks of 4 elements */
            for (const uint32_t j_end = j + ((npgtob - j) & -4); j != j_end; j += 4) {
                const __m256d ymm_b = _mm256_loadu_pd(&alphab[j]);
                const __m256d ymm_p = _mm256_add_pd(ymm_a, ymm_b);
                const __m256d ymm_ab = _mm256_mul_pd(ymm_a, ymm_b);
                const __m256d ymm_pqp = _mm256_add_pd(ymm_p, ymm_qmin);
                const __m256d ymm_t = _mm256_mul_pd(ymm_rminsq_qmin, _mm256_div_pd(ymm_p, ymm_pqp));
                const __m256d ymm_rhoab = _mm256_exp_pd(_mm256_div_pd(_mm256_mul_pd(ymm_ab, ymm_minus_rnabsq), ymm_p));

                /* pi/4 * f0(x) == x */
                const __m256d ymm_x = _mm256_blendv_pd(ymm_t, ymm_x0, _mm256_cmp_pd(ymm_t, _mm256_setzero_pd(), _CMP_EQ_OS));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m256d ymm_f0 = 
                    _mm256_add_pd(ymm_c0, _mm256_mul_pd(ymm_x, 
                        _mm256_add_pd(ymm_c1, _mm256_mul_pd(ymm_x,
                            _mm256_add_pd(ymm_c2, _mm256_mul_pd(ymm_x,
                                _mm256_add_pd(ymm_c3, _mm256_mul_pd(ymm_x,
                                    _mm256_add_pd(ymm_c4, _mm256_mul_pd(ymm_x, 
                                        ymm_c5)
                                ))
                            ))
                        ))
                    ))
                );
                const __m256d ymm_f = _mm256_min_pd(ymm_f0, ymm_one);
                const __m256d ymm_rhoab_csmaxcd = _mm256_mul_pd(ymm_rhoab, ymm_csmaxcd);
                const __m256d ymm_rhoab2_smaxcd2 = _mm256_mul_pd(ymm_rhoab_csmaxcd, ymm_rhoab_csmaxcd);
                const __m256d ymm_ab_f = _mm256_mul_pd(ymm_ab, ymm_f);
                const __m256d ymm_ab_f_rhoab2_smaxcd2 = _mm256_mul_pd(ymm_ab_f, ymm_rhoab2_smaxcd2);
                const __m256d ymm_ab2_f2_rhoab4_smaxcd4 = _mm256_mul_pd(ymm_ab_f_rhoab2_smaxcd2, ymm_ab_f_rhoab2_smaxcd2);
                const __m256d ymm_ab3_f2_rhoab4_smaxcd4 = _mm256_mul_pd(ymm_ab, ymm_ab2_f2_rhoab4_smaxcd4);
                const __m256d ymm_p2 = _mm256_mul_pd(ymm_p, ymm_p);
                const __m256d ymm_x_pqp = _mm256_mul_pd(ymm_x, ymm_pqp);
                const __m256d ymm_x_pqp_p2 = _mm256_mul_pd(ymm_x_pqp, ymm_p2);
                const __m256d ymm_x2_pqp2_p4 = _mm256_mul_pd(ymm_x_pqp_p2, ymm_x_pqp_p2);
                const __m256d ymm_cmp_mask = _mm256_cmp_pd(ymm_ab3_f2_rhoab4_smaxcd4, ymm_x2_pqp2_p4, _CMP_GE_OS);
                uint32_t mask = _mm256_movemask_pd(ymm_cmp_mask);
                const __m128i xmm_pack = _mm_load_si128((const __m128i*)&xmm_pack_table[mask*16]);
                _mm_storeu_si128((__m128i*)&prima[nij], xmm_i);
                _mm_storeu_si128((__m128i*)&primb[nij], _mm_shuffle_epi8(xmm_j, xmm_pack));
                const __m256d ymm_rhoab_out = _mm256_blendv_pd(_mm256_unpackhi_pd(ymm_rhoab, ymm_rhoab), ymm_rhoab, ymm_cmp_mask);
                const __m128d xmm_rhoab_out_lo = _mm256_castpd256_pd128(ymm_rhoab_out);
                const __m128d xmm_rhoab_out_hi = _mm256_extractf128_pd(ymm_rhoab_out, 1);
                _mm_storeu_pd(&rho[nij], xmm_rhoab_out_lo);
                _mm_storeu_pd(&rho[nij] + __builtin_popcount(mask & 0x3), xmm_rhoab_out_hi);
                nij += __builtin_popcount(mask);
                xmm_j = _mm_add_epi32(xmm_j, xmm_increment_j);
            }
            i += 1;
            /* Process the remainder of this j loop + start of the next j loop */
            if (j != npgtob) {
                __m256d ymm_merge_mask = ymm_merge_mask_table[npgtob - j];

                ymm_a = _mm256_blendv_pd(ymm_a, _mm256_broadcast_sd(&alphaa[i]), ymm_merge_mask);
                const __m128i xmm_sub_j = _mm_sub_epi32(xmm_j, xmm_npgtob);
                const __m128i xmm_wrap_j = _mm_min_epu32(xmm_j, xmm_sub_j);
                const __m128i xmm_wrap_i = _mm_sub_epi32(xmm_i, _mm_cmpeq_epi32(xmm_wrap_j, xmm_sub_j));

                __m256d ymm_b = _mm256_loadu_pd(&alphab[j]);
                if (i != npgtoa) {
                    const __m256d ymm_wrap_b = _mm256_maskload_pd(&alphab[(size_t)j - (size_t)npgtob], _mm256_castpd_si256(ymm_merge_mask));
                    ymm_b = _mm256_blendv_pd(ymm_b, ymm_wrap_b, ymm_merge_mask);
                }
                const __m256d ymm_p = _mm256_add_pd(ymm_a, ymm_b);
                const __m256d ymm_ab = _mm256_mul_pd(ymm_a, ymm_b);
                const __m256d ymm_pqp = _mm256_add_pd(ymm_p, ymm_qmin);
                const __m256d ymm_t = _mm256_mul_pd(ymm_rminsq_qmin, _mm256_div_pd(ymm_p, ymm_pqp));
                const __m256d ymm_rhoab = _mm256_exp_pd(_mm256_div_pd(_mm256_mul_pd(ymm_ab, ymm_minus_rnabsq), ymm_p));

                /* pi/4 * f0(x) == x */
                const __m256d ymm_x = _mm256_blendv_pd(ymm_t, ymm_x0, _mm256_cmp_pd(ymm_t, _mm256_setzero_pd(), _CMP_EQ_OS));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m256d ymm_f0 = 
                    _mm256_add_pd(ymm_c0, _mm256_mul_pd(ymm_x, 
                        _mm256_add_pd(ymm_c1, _mm256_mul_pd(ymm_x,
                            _mm256_add_pd(ymm_c2, _mm256_mul_pd(ymm_x,
                                _mm256_add_pd(ymm_c3, _mm256_mul_pd(ymm_x,
                                    _mm256_add_pd(ymm_c4, _mm256_mul_pd(ymm_x, 
                                        ymm_c5)
                                ))
                            ))
                        ))
                    ))
                );
                const __m256d ymm_f = _mm256_min_pd(ymm_f0, ymm_one);
                const __m256d ymm_rhoab_csmaxcd = _mm256_mul_pd(ymm_rhoab, ymm_csmaxcd);
                const __m256d ymm_rhoab2_smaxcd2 = _mm256_mul_pd(ymm_rhoab_csmaxcd, ymm_rhoab_csmaxcd);
                const __m256d ymm_ab_f = _mm256_mul_pd(ymm_ab, ymm_f);
                const __m256d ymm_ab_f_rhoab2_smaxcd2 = _mm256_mul_pd(ymm_ab_f, ymm_rhoab2_smaxcd2);
                const __m256d ymm_ab2_f2_rhoab4_smaxcd4 = _mm256_mul_pd(ymm_ab_f_rhoab2_smaxcd2, ymm_ab_f_rhoab2_smaxcd2);
                const __m256d ymm_ab3_f2_rhoab4_smaxcd4 = _mm256_mul_pd(ymm_ab, ymm_ab2_f2_rhoab4_smaxcd4);
                const __m256d ymm_p2 = _mm256_mul_pd(ymm_p, ymm_p);
                const __m256d ymm_x_pqp = _mm256_mul_pd(ymm_x, ymm_pqp);
                const __m256d ymm_x_pqp_p2 = _mm256_mul_pd(ymm_x_pqp, ymm_p2);
                const __m256d ymm_x2_pqp2_p4 = _mm256_mul_pd(ymm_x_pqp_p2, ymm_x_pqp_p2);
                __m256d ymm_cmp_mask = _mm256_cmp_pd(ymm_ab3_f2_rhoab4_smaxcd4, ymm_x2_pqp2_p4, _CMP_GE_OS);
                if (i == npgtoa)
                    ymm_cmp_mask = _mm256_andnot_pd(ymm_merge_mask, ymm_cmp_mask);
                uint32_t mask = _mm256_movemask_pd(ymm_cmp_mask);
                const __m128i xmm_pack = _mm_load_si128((const __m128i*)&xmm_pack_table[mask*16]);
                _mm_storeu_si128((__m128i*)&prima[nij], _mm_shuffle_epi8(xmm_wrap_i, xmm_pack));
                _mm_storeu_si128((__m128i*)&primb[nij], _mm_shuffle_epi8(xmm_wrap_j, xmm_pack));
                const __m256d ymm_rhoab_out = _mm256_blendv_pd(_mm256_unpackhi_pd(ymm_rhoab, ymm_rhoab), ymm_rhoab, ymm_cmp_mask);
                const __m128d xmm_rhoab_out_lo = _mm256_castpd256_pd128(ymm_rhoab_out);
                const __m128d xmm_rhoab_out_hi = _mm256_extractf128_pd(ymm_rhoab_out, 1);
                _mm_storeu_pd(&rho[nij], xmm_rhoab_out_lo);
                _mm_storeu_pd(&rho[nij] + __builtin_popcount(mask & 0x3), xmm_rhoab_out_hi);
                nij += __builtin_popcount(mask);

                j += 4;
                xmm_j = _mm_add_epi32(xmm_j, xmm_increment_j);
            }
            j -= npgtob;
            xmm_j = _mm_sub_epi32(xmm_j, xmm_npgtob);
        } while (i != npgtoa);
    }
#elif defined(__SSE4_1__)
    if ((npgtob % 2) == 0) {
        const __m128d xmm_minus_rnabsq = _mm_set1_pd(-rnabsq);
        const __m128d xmm_csmaxcd = _mm_set1_pd(csmaxcd);
        const __m128d xmm_qmin = _mm_set1_pd(qmin);
        const __m128d xmm_rminsq_qmin = _mm_mul_pd(xmm_qmin, _mm_set1_pd(rminsq));
        for (uint32_t i = 0; i < npgtoa; i += 1) {
            const __m128d xmm_a = _mm_loaddup_pd(&alphaa[i]);
            for (uint32_t j = 0; j < npgtob; j += 2) {
                const __m128d xmm_b = _mm_load_pd(&alphab[j]);
                const __m128d xmm_p = _mm_add_pd(xmm_a, xmm_b);
                const __m128d xmm_ab = _mm_mul_pd(xmm_a, xmm_b);
                const __m128d xmm_pqp = _mm_add_pd(xmm_p, xmm_qmin);
                const __m128d xmm_t = _mm_mul_pd(xmm_rminsq_qmin, _mm_div_pd(xmm_p, xmm_pqp));
                const __m128d xmm_rhoab = _mm_exp_pd(_mm_div_pd(_mm_mul_pd(xmm_ab, xmm_minus_rnabsq), xmm_p));

                /* pi/4 * f0(x) == x */
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(x0), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(c0),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(c1),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(c2),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(c3),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(c4),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(c5))
                                ))
                            ))
                        ))
                    ))
                );
                const __m128d xmm_f = _mm_min_pd(xmm_f0, _mm_set1_pd(1.0));
                const __m128d xmm_rhoab_csmaxcd = _mm_mul_pd(xmm_rhoab, xmm_csmaxcd);
                const __m128d xmm_rhoab2_smaxcd2 = _mm_mul_pd(xmm_rhoab_csmaxcd, xmm_rhoab_csmaxcd);
                const __m128d xmm_ab_f = _mm_mul_pd(xmm_ab, xmm_f);
                const __m128d xmm_ab_f_rhoab2_smaxcd2 = _mm_mul_pd(xmm_ab_f, xmm_rhoab2_smaxcd2);
                const __m128d xmm_ab2_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab_f_rhoab2_smaxcd2, xmm_ab_f_rhoab2_smaxcd2);
                const __m128d xmm_ab3_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab, xmm_ab2_f2_rhoab4_smaxcd4);
                const __m128d xmm_p2 = _mm_mul_pd(xmm_p, xmm_p);
                const __m128d xmm_x_pqp = _mm_mul_pd(xmm_x, xmm_pqp);
                const __m128d xmm_x_pqp_p2 = _mm_mul_pd(xmm_x_pqp, xmm_p2);
                const __m128d xmm_x2_pqp2_p4 = _mm_mul_pd(xmm_x_pqp_p2, xmm_x_pqp_p2);
                uint32_t mask = _mm_movemask_pd(_mm_cmpge_pd(xmm_ab3_f2_rhoab4_smaxcd4, xmm_x2_pqp2_p4));
                _mm_storel_pd(&rho[nij], xmm_rhoab);
                prima[nij] = i;
                primb[nij] = j;
                nij += mask & 1;
                mask >>= 1;
                _mm_storeh_pd(&rho[nij], xmm_rhoab);
                prima[nij] = i;
                primb[nij] = j + 1;
                nij += mask & 1;
            }
        }
    } else {
        const __m128d xmm_minus_rnabsq = _mm_set1_pd(-rnabsq);
        const __m128d xmm_csmaxcd = _mm_set1_pd(csmaxcd);
        const __m128d xmm_qmin = _mm_set1_pd(qmin);
        const __m128d xmm_rminsq_qmin = _mm_mul_pd(xmm_qmin, _mm_set1_pd(rminsq));
        for (uint32_t i = 0; i + 1 < npgtoa; i += 2) {
            __m128d xmm_a = _mm_loaddup_pd(&alphaa[i]);
            for (uint32_t j = 0; j < npgtob - 1; j += 2) {
                const __m128d xmm_b = _mm_load_pd(&alphab[j]);
                const __m128d xmm_p = _mm_add_pd(xmm_a, xmm_b);
                const __m128d xmm_ab = _mm_mul_pd(xmm_a, xmm_b);
                const __m128d xmm_pqp = _mm_add_pd(xmm_p, xmm_qmin);
                const __m128d xmm_t = _mm_mul_pd(xmm_rminsq_qmin, _mm_div_pd(xmm_p, xmm_pqp));
                const __m128d xmm_rhoab = _mm_exp_pd(_mm_div_pd(_mm_mul_pd(xmm_ab, xmm_minus_rnabsq), xmm_p));

                /* pi/4 * f0(x) == x */
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(x0), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(c0),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(c1),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(c2),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(c3),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(c4),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(c5))
                                ))
                            ))
                        ))
                    ))
                );
                const __m128d xmm_f = _mm_min_pd(xmm_f0, _mm_set1_pd(1.0));
                const __m128d xmm_rhoab_csmaxcd = _mm_mul_pd(xmm_rhoab, xmm_csmaxcd);
                const __m128d xmm_rhoab2_smaxcd2 = _mm_mul_pd(xmm_rhoab_csmaxcd, xmm_rhoab_csmaxcd);
                const __m128d xmm_ab_f = _mm_mul_pd(xmm_ab, xmm_f);
                const __m128d xmm_ab_f_rhoab2_smaxcd2 = _mm_mul_pd(xmm_ab_f, xmm_rhoab2_smaxcd2);
                const __m128d xmm_ab2_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab_f_rhoab2_smaxcd2, xmm_ab_f_rhoab2_smaxcd2);
                const __m128d xmm_ab3_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab, xmm_ab2_f2_rhoab4_smaxcd4);
                const __m128d xmm_p2 = _mm_mul_pd(xmm_p, xmm_p);
                const __m128d xmm_x_pqp = _mm_mul_pd(xmm_x, xmm_pqp);
                const __m128d xmm_x_pqp_p2 = _mm_mul_pd(xmm_x_pqp, xmm_p2);
                const __m128d xmm_x2_pqp2_p4 = _mm_mul_pd(xmm_x_pqp_p2, xmm_x_pqp_p2);
                uint32_t mask = _mm_movemask_pd(_mm_cmpge_pd(xmm_ab3_f2_rhoab4_smaxcd4, xmm_x2_pqp2_p4));
                _mm_storel_pd(&rho[nij], xmm_rhoab);
                prima[nij] = i;
                primb[nij] = j;
                nij += mask & 1;
                mask >>= 1;
                _mm_storeh_pd(&rho[nij], xmm_rhoab);
                prima[nij] = i;
                primb[nij] = j + 1;
                nij += mask & 1;
            }
            xmm_a = _mm_loadh_pd(xmm_a, &alphaa[i+1]);
            {
                const __m128d xmm_b = _mm_loadh_pd(_mm_load_sd(&alphab[npgtob-1]), &alphab[0]);
                const __m128d xmm_p = _mm_add_pd(xmm_a, xmm_b);
                const __m128d xmm_ab = _mm_mul_pd(xmm_a, xmm_b);
                const __m128d xmm_pqp = _mm_add_pd(xmm_p, xmm_qmin);
                const __m128d xmm_t = _mm_mul_pd(xmm_rminsq_qmin, _mm_div_pd(xmm_p, xmm_pqp));
                const __m128d xmm_rhoab = _mm_exp_pd(_mm_div_pd(_mm_mul_pd(xmm_ab, xmm_minus_rnabsq), xmm_p));

                /* pi/4 * f0(x) == x */
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(x0), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(c0),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(c1),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(c2),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(c3),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(c4),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(c5))
                                ))
                            ))
                        ))
                    ))
                );
                const __m128d xmm_f = _mm_min_pd(xmm_f0, _mm_set1_pd(1.0));
                const __m128d xmm_rhoab_csmaxcd = _mm_mul_pd(xmm_rhoab, xmm_csmaxcd);
                const __m128d xmm_rhoab2_smaxcd2 = _mm_mul_pd(xmm_rhoab_csmaxcd, xmm_rhoab_csmaxcd);
                const __m128d xmm_ab_f = _mm_mul_pd(xmm_ab, xmm_f);
                const __m128d xmm_ab_f_rhoab2_smaxcd2 = _mm_mul_pd(xmm_ab_f, xmm_rhoab2_smaxcd2);
                const __m128d xmm_ab2_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab_f_rhoab2_smaxcd2, xmm_ab_f_rhoab2_smaxcd2);
                const __m128d xmm_ab3_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab, xmm_ab2_f2_rhoab4_smaxcd4);
                const __m128d xmm_p2 = _mm_mul_pd(xmm_p, xmm_p);
                const __m128d xmm_x_pqp = _mm_mul_pd(xmm_x, xmm_pqp);
                const __m128d xmm_x_pqp_p2 = _mm_mul_pd(xmm_x_pqp, xmm_p2);
                const __m128d xmm_x2_pqp2_p4 = _mm_mul_pd(xmm_x_pqp_p2, xmm_x_pqp_p2);
                uint32_t mask = _mm_movemask_pd(_mm_cmpge_pd(xmm_ab3_f2_rhoab4_smaxcd4, xmm_x2_pqp2_p4));
                _mm_storel_pd(&rho[nij], xmm_rhoab);
                prima[nij] = i;
                primb[nij] = npgtob-1;
                nij += mask & 1;
                mask >>= 1;
                _mm_storeh_pd(&rho[nij], xmm_rhoab);
                prima[nij] = i + 1;
                primb[nij] = 0;
                nij += mask & 1;
            }
            xmm_a = _mm_unpackhi_pd(xmm_a, xmm_a);
            for (uint32_t j = 1; j < npgtob; j += 2) {
                const __m128d xmm_b = _mm_loadu_pd(&alphab[j]);
                const __m128d xmm_p = _mm_add_pd(xmm_a, xmm_b);
                const __m128d xmm_ab = _mm_mul_pd(xmm_a, xmm_b);
                const __m128d xmm_pqp = _mm_add_pd(xmm_p, xmm_qmin);
                const __m128d xmm_t = _mm_mul_pd(xmm_rminsq_qmin, _mm_div_pd(xmm_p, xmm_pqp));
                const __m128d xmm_rhoab = _mm_exp_pd(_mm_div_pd(_mm_mul_pd(xmm_ab, xmm_minus_rnabsq), xmm_p));

                /* pi/4 * f0(x) == x */
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(x0), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(c0),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(c1),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(c2),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(c3),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(c4),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(c5))
                                ))
                            ))
                        ))
                    ))
                );
                const __m128d xmm_f = _mm_min_pd(xmm_f0, _mm_set1_pd(1.0));
                const __m128d xmm_rhoab_csmaxcd = _mm_mul_pd(xmm_rhoab, xmm_csmaxcd);
                const __m128d xmm_rhoab2_smaxcd2 = _mm_mul_pd(xmm_rhoab_csmaxcd, xmm_rhoab_csmaxcd);
                const __m128d xmm_ab_f = _mm_mul_pd(xmm_ab, xmm_f);
                const __m128d xmm_ab_f_rhoab2_smaxcd2 = _mm_mul_pd(xmm_ab_f, xmm_rhoab2_smaxcd2);
                const __m128d xmm_ab2_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab_f_rhoab2_smaxcd2, xmm_ab_f_rhoab2_smaxcd2);
                const __m128d xmm_ab3_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab, xmm_ab2_f2_rhoab4_smaxcd4);
                const __m128d xmm_p2 = _mm_mul_pd(xmm_p, xmm_p);
                const __m128d xmm_x_pqp = _mm_mul_pd(xmm_x, xmm_pqp);
                const __m128d xmm_x_pqp_p2 = _mm_mul_pd(xmm_x_pqp, xmm_p2);
                const __m128d xmm_x2_pqp2_p4 = _mm_mul_pd(xmm_x_pqp_p2, xmm_x_pqp_p2);
                uint32_t mask = _mm_movemask_pd(_mm_cmpge_pd(xmm_ab3_f2_rhoab4_smaxcd4, xmm_x2_pqp2_p4));
                _mm_storel_pd(&rho[nij], xmm_rhoab);
                prima[nij] = i + 1;
                primb[nij] = j;
                nij += mask & 1;
                mask >>= 1;
                _mm_storeh_pd(&rho[nij], xmm_rhoab);
                prima[nij] = i + 1;
                primb[nij] = j + 1;
                nij += mask & 1;
            }
        }
        if ((npgtoa % 2) == 1) {
            __m128d xmm_a = _mm_loaddup_pd(&alphaa[npgtoa-1]);
            for (uint32_t j = 0; j < npgtob - 1; j += 2) {
                const __m128d xmm_b = _mm_load_pd(&alphab[j]);
                const __m128d xmm_p = _mm_add_pd(xmm_a, xmm_b);
                const __m128d xmm_ab = _mm_mul_pd(xmm_a, xmm_b);
                const __m128d xmm_pqp = _mm_add_pd(xmm_p, xmm_qmin);
                const __m128d xmm_t = _mm_mul_pd(xmm_rminsq_qmin, _mm_div_pd(xmm_p, xmm_pqp));
                const __m128d xmm_rhoab = _mm_exp_pd(_mm_div_pd(_mm_mul_pd(xmm_ab, xmm_minus_rnabsq), xmm_p));

                /* pi/4 * f0(x) == x */
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(x0), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(c0),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(c1),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(c2),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(c3),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(c4),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(c5))
                                ))
                            ))
                        ))
                    ))
                );
                const __m128d xmm_f = _mm_min_pd(xmm_f0, _mm_set1_pd(1.0));
                const __m128d xmm_rhoab_csmaxcd = _mm_mul_pd(xmm_rhoab, xmm_csmaxcd);
                const __m128d xmm_rhoab2_smaxcd2 = _mm_mul_pd(xmm_rhoab_csmaxcd, xmm_rhoab_csmaxcd);
                const __m128d xmm_ab_f = _mm_mul_pd(xmm_ab, xmm_f);
                const __m128d xmm_ab_f_rhoab2_smaxcd2 = _mm_mul_pd(xmm_ab_f, xmm_rhoab2_smaxcd2);
                const __m128d xmm_ab2_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab_f_rhoab2_smaxcd2, xmm_ab_f_rhoab2_smaxcd2);
                const __m128d xmm_ab3_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab, xmm_ab2_f2_rhoab4_smaxcd4);
                const __m128d xmm_p2 = _mm_mul_pd(xmm_p, xmm_p);
                const __m128d xmm_x_pqp = _mm_mul_pd(xmm_x, xmm_pqp);
                const __m128d xmm_x_pqp_p2 = _mm_mul_pd(xmm_x_pqp, xmm_p2);
                const __m128d xmm_x2_pqp2_p4 = _mm_mul_pd(xmm_x_pqp_p2, xmm_x_pqp_p2);
                uint32_t mask = _mm_movemask_pd(_mm_cmpge_pd(xmm_ab3_f2_rhoab4_smaxcd4, xmm_x2_pqp2_p4));
                _mm_storel_pd(&rho[nij], xmm_rhoab);
                prima[nij] = npgtoa-1;
                primb[nij] = j;
                nij += mask & 1;
                mask >>= 1;
                _mm_storeh_pd(&rho[nij], xmm_rhoab);
                prima[nij] = npgtoa-1;
                primb[nij] = j + 1;
                nij += mask & 1;
            }
            {
                const __m128d xmm_b = _mm_load_sd(&alphab[npgtob-1]);
                const __m128d xmm_p = _mm_add_pd(xmm_a, xmm_b);
                const __m128d xmm_ab = _mm_mul_pd(xmm_a, xmm_b);
                const __m128d xmm_pqp = _mm_add_pd(xmm_p, xmm_qmin);
                const __m128d xmm_t = _mm_mul_pd(xmm_rminsq_qmin, _mm_div_pd(xmm_p, xmm_pqp));
                const __m128d xmm_rhoab = _mm_exp_pd(_mm_div_pd(_mm_mul_pd(xmm_ab, xmm_minus_rnabsq), xmm_p));

                /* pi/4 * f0(x) == x */
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(x0), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(c0),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(c1),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(c2),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(c3),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(c4),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(c5))
                                ))
                            ))
                        ))
                    ))
                );
                const __m128d xmm_f = _mm_min_pd(xmm_f0, _mm_set1_pd(1.0));
                const __m128d xmm_rhoab_csmaxcd = _mm_mul_pd(xmm_rhoab, xmm_csmaxcd);
                const __m128d xmm_rhoab2_smaxcd2 = _mm_mul_pd(xmm_rhoab_csmaxcd, xmm_rhoab_csmaxcd);
                const __m128d xmm_ab_f = _mm_mul_pd(xmm_ab, xmm_f);
                const __m128d xmm_ab_f_rhoab2_smaxcd2 = _mm_mul_pd(xmm_ab_f, xmm_rhoab2_smaxcd2);
                const __m128d xmm_ab2_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab_f_rhoab2_smaxcd2, xmm_ab_f_rhoab2_smaxcd2);
                const __m128d xmm_ab3_f2_rhoab4_smaxcd4 = _mm_mul_pd(xmm_ab, xmm_ab2_f2_rhoab4_smaxcd4);
                const __m128d xmm_p2 = _mm_mul_pd(xmm_p, xmm_p);
                const __m128d xmm_x_pqp = _mm_mul_pd(xmm_x, xmm_pqp);
                const __m128d xmm_x_pqp_p2 = _mm_mul_pd(xmm_x_pqp, xmm_p2);
                const __m128d xmm_x2_pqp2_p4 = _mm_mul_pd(xmm_x_pqp_p2, xmm_x_pqp_p2);
                uint32_t mask = _mm_movemask_pd(_mm_cmpge_pd(xmm_ab3_f2_rhoab4_smaxcd4, xmm_x2_pqp2_p4));
                _mm_storel_pd(&rho[nij], xmm_rhoab);
                prima[nij] = npgtoa-1;
                primb[nij] = npgtob-1;
                nij += mask & 1;
            }
        }
    }
#else
    const double minus_rnabsq = -rnabsq;
    const double rminsq_qmin = rminsq * qmin;
    for (uint32_t i = 0; i < npgtoa; i += 1) {
        const double a = alphaa[i];
        for (uint32_t j = 0; j < npgtob; j += 1) {
            const double b = alphab[j];
            const double p = a + b;
            const double ab = a * b;
            const double pqp = p + qmin;
            const double pqpinv = 1.0 / pqp;
            const double rhoab = __builtin_exp(ab * minus_rnabsq / p);
            const double t = rminsq_qmin * (p / pqp);

            /* pi/4 * f0(x) == x */
            const double x = (t == 0.0) ? x0 : t;
            /* approximates square(erf(sqrt(x))) on [0, 5] */
            const double f0 = c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * c5))));
            const double f = f0 < 1.0 ? f0 : 1.0;
            if (ab*square(square(rhoab*csmaxcd) * (ab*f)) >= square((x*pqp)*square(p))) {
                rho[nij] = rhoab;
                prima[nij] = i;
                primb[nij] = j;
                nij += 1;
            }
        }
    }
#endif
    *nij_ptr = nij;
}

int erd__set_ij_kl_pairs(
    int npgtoa, int npgtob, int npgtoc, int npgtod,
    double xa, double ya, double za,
    double xb, double yb, double zb,
    double xc, double yc, double zc,
    double xd, double yd, double zd,
    double rnabsq, double rncdsq, double prefact,
    double *YEP_RESTRICT alphaa,
    double *YEP_RESTRICT alphab,
    double *YEP_RESTRICT alphac,
    double *YEP_RESTRICT alphad, int screen,
    int *YEP_RESTRICT empty, int *YEP_RESTRICT nij_ptr,
    int *YEP_RESTRICT nkl_ptr, int *YEP_RESTRICT prima,
    int *YEP_RESTRICT primb, int *YEP_RESTRICT primc,
    int *YEP_RESTRICT primd, double *YEP_RESTRICT rho)
{
    *empty = 0;

    uint32_t nij = 0;
    uint32_t nkl = 0;
    uint32_t padnij = 0;
    
    // if not screening
    if (!(screen)) {
        for (int i = 0; i < npgtoa; i += 1) {
            const double a = alphaa[i];
            for (int j = 0; j < npgtob; j += 1) {
                const double b = alphab[j];
                rho[nij] = __builtin_exp(-a * b * rnabsq / (a + b));
                prima[nij] = i;
                primb[nij] = j;
                nij += 1;
            }
        }
        padnij = PAD_LEN (nij);
        
        for (int k = 0; k < npgtoc; k += 1) {
            const double c = alphac[k];
            for (int l = 0; l < npgtod; l += 1) {
                const double d = alphad[l];
                rho[padnij + nkl] = __builtin_exp(-c * d * rncdsq / (c + d));
                primc[nkl] = k;
                primd[nkl] = l;
                nkl += 1;
            }
        }
        return 0;
    }

    // compute min
    const double rminsq = erd__dsqmin_line_segments(xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd);

    const double a = vector_min(alphaa, npgtoa);
    const double b = vector_min(alphab, npgtob);
    const double c = vector_min(alphac, npgtoc);
    const double d = vector_min(alphad, npgtod);

    const double pmin = a + b;
    const double qmin = c + d;
    const double abmin = a * b;
    const double cdmin = c * d;
    const double pinv = 1.0 / pmin;
    const double qinv = 1.0 / qmin;
    const double smaxab = prefact * pow3o4(abmin) * __builtin_exp(-abmin * rnabsq * pinv) * pinv;
    const double smaxcd = prefact * pow3o4(cdmin) * __builtin_exp(-cdmin * rncdsq * qinv) * qinv;

    /* ...perform K2 primitive screening on A,B part. */
    set_pairs(npgtoa, npgtob, rnabsq, alphaa, alphab, &nij, prima, primb, rho, qmin, smaxcd, rminsq);
    if (nij == 0) {
        *empty = 1;

        *nij_ptr = nij;
        *nkl_ptr = nkl;
        return 0;
    }
    padnij = PAD_LEN (nij);
    
    set_pairs(npgtoc, npgtod, rncdsq, alphac, alphad, &nkl, primc, primd, &rho[padnij], pmin, smaxab, rminsq);
    if (nkl == 0) {
        *empty = 1;

        *nij_ptr = nij;
        *nkl_ptr = nkl;
        return 0;
    }

    *nij_ptr = nij;
    *nkl_ptr = nkl;
    return 0;
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
