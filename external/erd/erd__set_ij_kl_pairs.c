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

static YEP_NOINLINE void set_pairs(
    uint32_t npgtoa, uint32_t npgtob, double rnabsq,
    double *YEP_RESTRICT alphaa, double *YEP_RESTRICT alphab,
    uint32_t *YEP_RESTRICT nij_ptr,
    int *YEP_RESTRICT prima, int *YEP_RESTRICT primb,
    double *YEP_RESTRICT rho,
    double qmin, double smaxcd, double rminsq)
{
    const double csmaxcd = smaxcd * (0x1.C5BF891B4EF6Bp-1 / TOL);
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
    uint32_t nij = 0;
#if defined(__SSE4_1__)
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
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(0.043279823828983313427), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(0.00407565479675641252609948418213623998912072458154598242200929737220744709843048711660458),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(1.20764773784216992679627679198009535561546045160994177718758068495792765807646134349002377),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(-0.668856165551048279761732866991215362420284709052392489348530890378707598253294818725157054),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(0.198795058149610610004637933292818184140451276354350294406659821950295705343602350343494287),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(-0.0302279795858857468323941559145865292550787118904121714353356215879652816000677315137445060),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(0.00183088802814224943984687247959824735309951766250626784494885964622113988142120217488811367))
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
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(0.043279823828983313427), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(0.00407565479675641252609948418213623998912072458154598242200929737220744709843048711660458),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(1.20764773784216992679627679198009535561546045160994177718758068495792765807646134349002377),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(-0.668856165551048279761732866991215362420284709052392489348530890378707598253294818725157054),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(0.198795058149610610004637933292818184140451276354350294406659821950295705343602350343494287),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(-0.0302279795858857468323941559145865292550787118904121714353356215879652816000677315137445060),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(0.00183088802814224943984687247959824735309951766250626784494885964622113988142120217488811367))
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
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(0.043279823828983313427), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(0.00407565479675641252609948418213623998912072458154598242200929737220744709843048711660458),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(1.20764773784216992679627679198009535561546045160994177718758068495792765807646134349002377),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(-0.668856165551048279761732866991215362420284709052392489348530890378707598253294818725157054),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(0.198795058149610610004637933292818184140451276354350294406659821950295705343602350343494287),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(-0.0302279795858857468323941559145865292550787118904121714353356215879652816000677315137445060),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(0.00183088802814224943984687247959824735309951766250626784494885964622113988142120217488811367))
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
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(0.043279823828983313427), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(0.00407565479675641252609948418213623998912072458154598242200929737220744709843048711660458),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(1.20764773784216992679627679198009535561546045160994177718758068495792765807646134349002377),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(-0.668856165551048279761732866991215362420284709052392489348530890378707598253294818725157054),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(0.198795058149610610004637933292818184140451276354350294406659821950295705343602350343494287),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(-0.0302279795858857468323941559145865292550787118904121714353356215879652816000677315137445060),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(0.00183088802814224943984687247959824735309951766250626784494885964622113988142120217488811367))
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
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(0.043279823828983313427), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(0.00407565479675641252609948418213623998912072458154598242200929737220744709843048711660458),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(1.20764773784216992679627679198009535561546045160994177718758068495792765807646134349002377),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(-0.668856165551048279761732866991215362420284709052392489348530890378707598253294818725157054),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(0.198795058149610610004637933292818184140451276354350294406659821950295705343602350343494287),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(-0.0302279795858857468323941559145865292550787118904121714353356215879652816000677315137445060),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(0.00183088802814224943984687247959824735309951766250626784494885964622113988142120217488811367))
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
                const __m128d xmm_x = _mm_blendv_pd(xmm_t, _mm_set1_pd(0.043279823828983313427), _mm_cmpeq_pd(xmm_t, _mm_setzero_pd()));
                /* approximates square(erf(sqrt(x))) on [0, 5] */
                const __m128d xmm_f0 = _mm_add_pd(
                    _mm_set1_pd(0.00407565479675641252609948418213623998912072458154598242200929737220744709843048711660458),
                    _mm_mul_pd(xmm_x, _mm_add_pd(
                        _mm_set1_pd(1.20764773784216992679627679198009535561546045160994177718758068495792765807646134349002377),
                        _mm_mul_pd(xmm_x, _mm_add_pd(
                            _mm_set1_pd(-0.668856165551048279761732866991215362420284709052392489348530890378707598253294818725157054),
                            _mm_mul_pd(xmm_x, _mm_add_pd(
                                _mm_set1_pd(0.198795058149610610004637933292818184140451276354350294406659821950295705343602350343494287),
                                _mm_mul_pd(xmm_x, _mm_add_pd(
                                    _mm_set1_pd(-0.0302279795858857468323941559145865292550787118904121714353356215879652816000677315137445060),
                                        _mm_mul_pd(xmm_x, 
                                            _mm_set1_pd(0.00183088802814224943984687247959824735309951766250626784494885964622113988142120217488811367))
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
            const double x = (t == 0.0) ? 0.043279823828983313427 : t;
            /* approximates square(erf(sqrt(x))) on [0, 5] */
            const double f0 = 0.00407565479675641252609948418213623998912072458154598242200929737220744709843048711660458
                + x * (1.20764773784216992679627679198009535561546045160994177718758068495792765807646134349002377
                    + x * (-0.668856165551048279761732866991215362420284709052392489348530890378707598253294818725157054
                        + x * (0.198795058149610610004637933292818184140451276354350294406659821950295705343602350343494287
                            + x * (-0.0302279795858857468323941559145865292550787118904121714353356215879652816000677315137445060
                                + x * 0.00183088802814224943984687247959824735309951766250626784494885964622113988142120217488811367
                            )
                        )
                    )
                );
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
