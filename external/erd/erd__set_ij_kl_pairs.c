#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "boys.h"
#include "erd.h"
#define ERD_TABLE_FREE_BOYS_FUNCTIONS

#define TOL 1e-14

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
    int npgtoa, int npgtob, double rnabsq,
    double *YEP_RESTRICT alphaa, double *YEP_RESTRICT alphab,
    uint32_t *YEP_RESTRICT nij_ptr,
    int *YEP_RESTRICT prima, int *YEP_RESTRICT primb,
    double *YEP_RESTRICT rho,
    double qmin, double smaxcd, double rminsq)
{
    const double smaxcd_scaled = smaxcd * (0x1.C5BF891B4EF6Bp-1 / TOL);
    uint32_t nij = 0;
    for (int i = 0; i < npgtoa; i += 1) {
        const double a = alphaa[i];
        for (int j = 0; j < npgtob; j += 1) {
           // __assume_aligned (alphab, 64);
            const double b = alphab[j];
            const double p = a + b;
            const double ab = a * b;
            const double pinv = 1.0 / p;
            const double pqpinv = 1.0 / (p + qmin);
            const double rhoab = __builtin_exp(-ab * rnabsq * pinv);
            const double t = rminsq * p * qmin * pqpinv;

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
            const double f0_bounded = f0 < 1.0 ? f0 : 1.0;
            if (ab*ab*ab*pow4(rhoab*smaxcd_scaled*pinv) * square(f0_bounded * pqpinv) >= square(x)) {
                rho[nij] = rhoab;
                prima[nij] = i;
                primb[nij] = j;
                nij += 1;
            }
        }
    }
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

        for (int k = 0; k < npgtoc; k += 1) {
            const double c = alphac[k];
            for (int l = 0; l < npgtod; l += 1) {
                const double d = alphad[l];
                rho[nij + nkl] = __builtin_exp(-c * d * rncdsq / (c + d));
                primc[nkl] = k;
                primd[nkl] = l;
                nkl += 1;
            }
        }
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

    set_pairs(npgtoc, npgtod, rncdsq, alphac, alphad, &nkl, primc, primd, &rho[nij], pmin, smaxab, rminsq);
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
