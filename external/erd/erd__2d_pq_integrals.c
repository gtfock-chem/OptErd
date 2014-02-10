#include <stdio.h>
#include <immintrin.h>

#include "erd.h"

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif


/*  Copyright (c) 2003-2010 University of Florida */

/*  This program is free software; you can redistribute it and/or modify */
/*  it under the terms of the GNU General Public License as published by */
/*  the Free Software Foundation; either version 2 of the License, or */
/*  (at your option) any later version. */
/*  This program is distributed in the hope that it will be useful, */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*  GNU General Public License for more details. */
/*  The GNU General Public License is included in this distribution */
/*  in the file COPYRIGHT. */
/* Subroutine */ int
erd__2d_pq_integrals (int shellp, int shellq,
                        int ngqexq, double * wts, double * b00,
                        double * b01, double * b10, double * c00x,
                        double * c00y, double * c00z,
                        double * d00x, double * d00y,
                        double * d00z, int case2d,
                        double * int2dx, double * int2dy,
                        double * int2dz)
{
    /* System generated locals */
    int int2d_dim1, int2d_dim2;

    /* Local variables */
    int i, k, n, n1;
    double b0, b1;
    double weight;

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__2D_PQ_INTEGRALS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation calculates a full table of 2D PQ X,Y,Z */
/*                integrals using the Rys vertical recurrence scheme */
/*                VRR explained below. */

/*                The Rys weight is multiplied to the 2DX PQ integral */
/*                to reduce overall FLOP count. Note, that the Rys weight */
/*                factor needs to be introduced only three times for the */
/*                starting 2DX PQ integrals for the recurrence scheme, */
/*                namely to the (0,0), (1,0) and (0,1) elements. The */
/*                weight factor is then automatically propagated */
/*                through the vertical transfer equations (see below). */
/*                The recurrence scheme VRR is due to Rys, Dupuis and */
/*                King, J. Comp. Chem. 4, p.154-157 (1983). */


/*                   INT2D (0,0) = 1.D0    (* WEIGHT for the 2DX case) */
/*                   INT2D (1,0) = C00     (* WEIGHT for the 2DX case) */
/*                   INT2D (0,1) = D00     (* WEIGHT for the 2DX case) */

/*                   For I = 1,...,SHELLP-1 */
/*                       INT2D (I+1,0) = I * B10 * INT2D (I-1,0) */
/*                                         + C00 * INT2D (I,0) */
/*                   For K = 1,...,SHELLQ-1 */
/*                       INT2D (0,K+1) = K * B01 * INT2D (0,K-1) */
/*                                         + D00 * INT2D (0,K) */
/*                   For I = 1,...,SHELLP */
/*                       INT2D (I,1)   = I * B00 * INT2D (I-1,0) */
/*                                         + D00 * INT2D (I,0) */
/*                   For K = 2,...,SHELLQ */
/*                       INT2D (1,K)   = K * B00 * INT2D (0,K-1) */
/*                                         + C00 * INT2D (0,K) */
/*                   For K = 2,...,SHELLQ */
/*                   For I = 2,...,SHELLP */
/*                       INT2D (I,K)   = (I-1) * B10 * INT2D (I-2,K) */
/*                                         + K * B00 * INT2D (I-1,K-1) */
/*                                             + C00 * INT2D (I-1,K) */


/*                The 2D PQ integrals are calculated for all roots (info */
/*                already present in transmitted VRR coefficients!) and */
/*                for all exponent quadruples simultaneously and placed */
/*                into a 3-dimensional array. */


/*                  Input: */

/*                    SHELLx      =  maximum shell type for electrons */
/*                                   1 and 2 (x = P,Q) */
/*                    NGQEXQ      =  product of # of gaussian quadrature */
/*                                   points times exponent quadruplets */
/*                    WTS         =  all quadrature weights */
/*                    B00,B01,B10 =  VRR expansion coefficients */
/*                                   (cartesian coordinate independent) */
/*                    C00x,D00x   =  cartesian coordinate dependent */
/*                                   VRR expansion coefficients */
/*                                   (x = X,Y,Z) */
/*                    CASE2D      =  logical flag for simplifications */
/*                                   in 2D integral evaluation for */
/*                                   low quantum numbers */


/*                  Output: */

/*                    INT2Dx      =  all 2D PQ integrals for each */
/*                                   cartesian component (x = X,Y,Z) */


/*  AUTHOR      : Norbert Flocke */
/* ------------------------------------------------------------------------ */



/*             ...include files and declare variables. */




/* ------------------------------------------------------------------------ */


/*             ...jump according to the 4 different cases that can arise: */

/*                  P-shell = s- or higher angular momentum */
/*                  Q-shell = s- or higher angular momentum */

/*                each leading to simplifications in the VRR formulas. */
/*                The case present has been evaluated outside this */
/*                routine and is transmitted via argument. */


    /* Parameter adjustments */
    int2d_dim1 = ngqexq;
    int2d_dim2 = shellp + 1;

#ifdef __AVX__
  __m256d one_256 = _mm256_set1_pd(1.0);
#endif

    //int2dx_dim3 = (shellp + 1) * (shellq + 1);
    //printf("%d, %d, ngqexq = %d\n", shellp, shellq, ngqexq);
    /* Function Body */
    switch (case2d)
    {
    case 1:
        goto L1;
    case 2:
        goto L3;
    case 3:
        goto L3;
    case 4:
        goto L2;
    case 5:
        goto L4;
    case 6:
        goto L4;
    case 7:
        goto L2;
    case 8:
        goto L4;
    case 9:
        goto L4;
    }


/*             ...the case P = s-shell and Q = s-shell. */

  L1:
    for (n = 0; n < ngqexq; n+=SIMDW)
    {
#ifdef __AVX__
        __m256d wts_256 = _mm256_load_pd(&wts[n]);
        _mm256_store_pd(&int2dx[n], wts_256);
        _mm256_store_pd(&int2dy[n], one_256);
        _mm256_store_pd(&int2dz[n], one_256);
#else
#pragma vector aligned
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            int2dx[n + n1] = wts[n + n1];
            int2dy[n + n1] = 1.;
            int2dz[n + n1] = 1.;
        }
#endif
    }
    return 0;


/*             ...the cases P = s-shell and Q >= p-shell. */
/*                Evaluate I=0 and K=0,1. */


  L2:
    for (n = 0; n < ngqexq; n+=SIMDW)
    {

#ifdef __AVX__
        __m256d int2dx_0_256, int2dx_1_256, int2dx_2_256;
        __m256d int2dy_0_256, int2dy_1_256, int2dy_2_256;
        __m256d int2dz_0_256, int2dz_1_256, int2dz_2_256;
        __m256d d00x_256, d00y_256, d00z_256;

        int2dx_2_256 = _mm256_load_pd(&wts[n]);
        _mm256_store_pd(&int2dx[n], int2dx_2_256);
        d00x_256 = _mm256_load_pd(&d00x[n]);
        int2dx_1_256 = _mm256_mul_pd(d00x_256, int2dx_2_256);
        _mm256_store_pd(&int2dx[n + int2d_dim2 * int2d_dim1], int2dx_1_256);

        int2dy_2_256 = one_256;
        _mm256_store_pd(&int2dy[n], int2dy_2_256);
        d00y_256 = _mm256_load_pd(&d00y[n]);
        int2dy_1_256 = d00y_256;
        _mm256_store_pd(&int2dy[n + int2d_dim2 * int2d_dim1], int2dy_1_256);

        int2dz_2_256 = one_256;
        _mm256_store_pd(&int2dz[n], int2dz_2_256);
        d00z_256 = _mm256_load_pd(&d00z[n]);
        int2dz_1_256 = d00z_256;
        _mm256_store_pd(&int2dz[n + int2d_dim2 * int2d_dim1], int2dz_1_256);

#else
        double int2dx_0[SIMDW], int2dx_1[SIMDW], int2dx_2[SIMDW];
        double int2dy_0[SIMDW], int2dy_1[SIMDW], int2dy_2[SIMDW];
        double int2dz_0[SIMDW], int2dz_1[SIMDW], int2dz_2[SIMDW];
#pragma vector aligned
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            weight = wts[n + n1];
            int2dx[n + n1] = int2dx_2[n1] = weight;
            int2dy[n + n1] = int2dy_2[n1] = 1.;
            int2dz[n + n1] = int2dz_2[n1] = 1.;
        }
#pragma vector aligned
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            weight = wts[n + n1];
            int2dx[n + n1 + int2d_dim2 * int2d_dim1] = int2dx_1[n1] = d00x[n + n1] * weight;
            int2dy[n + n1 + int2d_dim2 * int2d_dim1] = int2dy_1[n1] = d00y[n + n1];
            int2dz[n + n1 + int2d_dim2 * int2d_dim1] = int2dz_1[n1] = d00z[n + n1];
        }
#endif

/*             ...evaluate I=0 and K=2,SHELLQ (if any). */
        for (k = 2; k <= shellq; ++k)
        {
            double k1 = k - 1;

#ifdef __AVX__
            __m256d k1_256 = _mm256_broadcast_sd(&k1);
            __m256d b01_256 = _mm256_load_pd(&b01[n]);
            __m256d b1_256 = _mm256_mul_pd(k1_256, b01_256);

            int2dx_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dx_2_256),
                                         _mm256_mul_pd(d00x_256, int2dx_1_256));
            int2dx_2_256 = int2dx_1_256;
            int2dx_1_256 = int2dx_0_256;
            _mm256_store_pd(&int2dx[n + k * int2d_dim2 * int2d_dim1], int2dx_0_256);

            int2dy_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dy_2_256),
                                         _mm256_mul_pd(d00y_256, int2dy_1_256));
            int2dy_2_256 = int2dy_1_256;
            int2dy_1_256 = int2dy_0_256;
            _mm256_store_pd(&int2dy[n + k * int2d_dim2 * int2d_dim1], int2dy_0_256);

            int2dz_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dz_2_256),
                                         _mm256_mul_pd(d00z_256, int2dz_1_256));
            int2dz_2_256 = int2dz_1_256;
            int2dz_1_256 = int2dz_0_256;
            _mm256_store_pd(&int2dz[n + k * int2d_dim2 * int2d_dim1], int2dz_0_256);

#else
#pragma vector aligned
            for(n1 = 0; n1 < SIMDW; n1++)
            {
                b1 = k1 * b01[n + n1];
                int2dx_0[n1] = b1 * int2dx_2[n1] + d00x[n + n1] * int2dx_1[n1];
                int2dx_2[n1] = int2dx_1[n1];
                int2dx_1[n1] = int2dx_0[n1];
                int2dx[n + n1 + k * int2d_dim2 * int2d_dim1] = int2dx_0[n1];

                int2dy_0[n1] = b1 * int2dy_2[n1] + d00y[n + n1] * int2dy_1[n1];
                int2dy_2[n1] = int2dy_1[n1];
                int2dy_1[n1] = int2dy_0[n1];
                int2dy[n + n1 + k * int2d_dim2 * int2d_dim1] = int2dy_0[n1];

                int2dz_0[n1] = b1 * int2dz_2[n1] + d00z[n + n1] * int2dz_1[n1];
                int2dz_2[n1] = int2dz_1[n1];
                int2dz_1[n1] = int2dz_0[n1];
                int2dz[n + n1 + k * int2d_dim2 * int2d_dim1] = int2dz_0[n1];
            }
#endif
        }
    }
    return 0;


/*             ...the cases P >= p-shell and Q = s-shell. */
/*                Evaluate I=0,1 and K=0. */


  L3:
    for (n = 0; n < ngqexq; n+=SIMDW)
    {
#ifdef __AVX__
        __m256d int2dx_0_256, int2dx_1_256, int2dx_2_256;
        __m256d int2dy_0_256, int2dy_1_256, int2dy_2_256;
        __m256d int2dz_0_256, int2dz_1_256, int2dz_2_256;
        __m256d c00x_256, c00y_256, c00z_256;

        int2dx_2_256 = _mm256_load_pd(&wts[n]);
        _mm256_store_pd(&int2dx[n], int2dx_2_256);
        c00x_256 = _mm256_load_pd(&c00x[n]);
        int2dx_1_256 = _mm256_mul_pd(c00x_256, int2dx_2_256);
        _mm256_store_pd(&int2dx[n + int2d_dim1], int2dx_1_256);

        int2dy_2_256 = one_256;
        _mm256_store_pd(&int2dy[n], int2dy_2_256);
        c00y_256 = _mm256_load_pd(&c00y[n]);
        int2dy_1_256 = c00y_256;
        _mm256_store_pd(&int2dy[n + int2d_dim1], int2dy_1_256);

        int2dz_2_256 = one_256;
        _mm256_store_pd(&int2dz[n], int2dz_2_256);
        c00z_256 = _mm256_load_pd(&c00z[n]);
        int2dz_1_256 = c00z_256;
        _mm256_store_pd(&int2dz[n + int2d_dim1], int2dz_1_256);

#else
        double int2dx_0[SIMDW], int2dx_1[SIMDW], int2dx_2[SIMDW];
        double int2dy_0[SIMDW], int2dy_1[SIMDW], int2dy_2[SIMDW];
        double int2dz_0[SIMDW], int2dz_1[SIMDW], int2dz_2[SIMDW];

#pragma vector aligned
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            weight = wts[n + n1];
            int2dx[n + n1] = int2dx_2[n1] = weight;
            int2dy[n + n1] = int2dy_2[n1] = 1.;
            int2dz[n + n1] = int2dz_2[n1] = 1.;
        }
#pragma vector aligned
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            weight = wts[n + n1];
            int2dx[n + n1 + int2d_dim1] = int2dx_1[n1] = c00x[n + n1] * weight;
            int2dy[n + n1 + int2d_dim1] = int2dy_1[n1] = c00y[n + n1];
            int2dz[n + n1 + int2d_dim1] = int2dz_1[n1] = c00z[n + n1];
        }
#endif
/*             ...evaluate I=2,SHELLP (if any) and K=0. */

        for (i = 2; i <= shellp; ++i)
        {
            double i1 = i - 1;

#ifdef __AVX__

            __m256d i1_256 = _mm256_broadcast_sd(&i1);
            __m256d b10_256 = _mm256_load_pd(&b10[n]);
            __m256d b1_256 = _mm256_mul_pd(i1_256, b10_256);

            int2dx_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dx_2_256),
                                         _mm256_mul_pd(c00x_256, int2dx_1_256));
            int2dx_2_256 = int2dx_1_256;
            int2dx_1_256 = int2dx_0_256;
            _mm256_store_pd(&int2dx[n + i * int2d_dim1], int2dx_0_256);

            int2dy_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dy_2_256),
                                         _mm256_mul_pd(c00y_256, int2dy_1_256));
            int2dy_2_256 = int2dy_1_256;
            int2dy_1_256 = int2dy_0_256;
            _mm256_store_pd(&int2dy[n + i * int2d_dim1], int2dy_0_256);

            int2dz_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dz_2_256),
                                         _mm256_mul_pd(c00z_256, int2dz_1_256));
            int2dz_2_256 = int2dz_1_256;
            int2dz_1_256 = int2dz_0_256;
            _mm256_store_pd(&int2dz[n + i * int2d_dim1], int2dz_0_256);

#else
#pragma vector aligned
            for(n1 = 0; n1 < SIMDW; n1++)
            {
                b1 = i1 * b10[n + n1];
                int2dx_0[n1] = b1 * int2dx_2[n1] + c00x[n + n1] * int2dx_1[n1];
                int2dx_2[n1] = int2dx_1[n1];
                int2dx_1[n1] = int2dx_0[n1];
                int2dx[n + n1 + i * int2d_dim1] = int2dx_0[n1];

                int2dy_0[n1] = b1 * int2dy_2[n1] + c00y[n + n1] * int2dy_1[n1];
                int2dy_2[n1] = int2dy_1[n1];
                int2dy_1[n1] = int2dy_0[n1];
                int2dy[n + n1 + i * int2d_dim1] = int2dy_0[n1];

                int2dz_0[n1] = b1 * int2dz_2[n1] + c00z[n + n1] * int2dz_1[n1];
                int2dz_2[n1] = int2dz_1[n1];
                int2dz_1[n1] = int2dz_0[n1];
                int2dz[n + n1 + i * int2d_dim1] = int2dz_0[n1];
            }
#endif
        }
    }
    return 0;


/*             ...the cases P >= p-shell and Q >= p-shell. */
/*                Evaluate I=0,SHELLP       I=0 */
/*                         K=0        and   K=0,SHELLQ */


  L4:
    for (n = 0; n < ngqexq; n+=SIMDW)
    {
#ifdef __AVX__
        __m256d int2dx_0_256, int2dx_i1_256, int2dx_k1_256, int2dx_2_256;
        __m256d int2dy_0_256, int2dy_i1_256, int2dy_k1_256, int2dy_2_256;
        __m256d int2dz_0_256, int2dz_i1_256, int2dz_k1_256, int2dz_2_256;
        __m256d c00x_256, c00y_256, c00z_256;

        int2dx_2_256 = _mm256_load_pd(&wts[n]);
        _mm256_store_pd(&int2dx[n], int2dx_2_256);
        c00x_256 = _mm256_load_pd(&c00x[n]);
        int2dx_i1_256 = _mm256_mul_pd(c00x_256, int2dx_2_256);
        _mm256_store_pd(&int2dx[n + int2d_dim1], int2dx_i1_256);

        int2dy_2_256 = one_256;
        _mm256_store_pd(&int2dy[n], int2dy_2_256);
        c00y_256 = _mm256_load_pd(&c00y[n]);
        int2dy_i1_256 = c00y_256;
        _mm256_store_pd(&int2dy[n + int2d_dim1], int2dy_i1_256);

        int2dz_2_256 = one_256;
        _mm256_store_pd(&int2dz[n], int2dz_2_256);
        c00z_256 = _mm256_load_pd(&c00z[n]);
        int2dz_i1_256 = c00z_256;
        _mm256_store_pd(&int2dz[n + int2d_dim1], int2dz_i1_256);

#else
        double int2dx_0[SIMDW], int2dx_i1[SIMDW], int2dx_k1[SIMDW], int2dx_2[SIMDW];
        double int2dy_0[SIMDW], int2dy_i1[SIMDW], int2dy_k1[SIMDW], int2dy_2[SIMDW];
        double int2dz_0[SIMDW], int2dz_i1[SIMDW], int2dz_k1[SIMDW], int2dz_2[SIMDW];

#pragma vector aligned
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            weight = wts[n + n1];
            int2dx[n + n1] = int2dx_2[n1] = weight;
            int2dy[n + n1] = int2dy_2[n1] = 1.;
            int2dz[n + n1] = int2dz_2[n1] = 1.;
        }
#pragma vector aligned
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            weight = wts[n + n1];
            int2dx[n + n1 + int2d_dim1] = int2dx_i1[n1] = c00x[n + n1] * weight;
            int2dy[n + n1 + int2d_dim1] = int2dy_i1[n1] = c00y[n + n1];
            int2dz[n + n1 + int2d_dim1] = int2dz_i1[n1] = c00z[n + n1];
        }
#endif

        for (i = 2; i <= shellp; ++i)
        {
            double i1 = i - 1;

#ifdef __AVX__

            __m256d i1_256 = _mm256_broadcast_sd(&i1);
            __m256d b10_256 = _mm256_load_pd(&b10[n]);
            __m256d b1_256 = _mm256_mul_pd(i1_256, b10_256);

            int2dx_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dx_2_256),
                    _mm256_mul_pd(c00x_256, int2dx_i1_256));
            int2dx_2_256 = int2dx_i1_256;
            int2dx_i1_256 = int2dx_0_256;
            _mm256_store_pd(&int2dx[n + i * int2d_dim1], int2dx_0_256);

            int2dy_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dy_2_256),
                    _mm256_mul_pd(c00y_256, int2dy_i1_256));
            int2dy_2_256 = int2dy_i1_256;
            int2dy_i1_256 = int2dy_0_256;
            _mm256_store_pd(&int2dy[n + i * int2d_dim1], int2dy_0_256);

            int2dz_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dz_2_256),
                    _mm256_mul_pd(c00z_256, int2dz_i1_256));
            int2dz_2_256 = int2dz_i1_256;
            int2dz_i1_256 = int2dz_0_256;
            _mm256_store_pd(&int2dz[n + i * int2d_dim1], int2dz_0_256);

#else
#pragma vector aligned
            for(n1 = 0; n1 < SIMDW; n1++)
            {
                b1 = i1 * b10[n + n1];
                int2dx_0[n1] = b1 * int2dx_2[n1] + c00x[n + n1] * int2dx_i1[n1];
                int2dx_2[n1] = int2dx_i1[n1];
                int2dx_i1[n1] = int2dx_0[n1];
                int2dx[n + n1 + i * int2d_dim1] = int2dx_0[n1];

                int2dy_0[n1] = b1 * int2dy_2[n1] + c00y[n + n1] * int2dy_i1[n1];
                int2dy_2[n1] = int2dy_i1[n1];
                int2dy_i1[n1] = int2dy_0[n1];
                int2dy[n + n1 + i * int2d_dim1] = int2dy_0[n1];

                int2dz_0[n1] = b1 * int2dz_2[n1] + c00z[n + n1] * int2dz_i1[n1];
                int2dz_2[n1] = int2dz_i1[n1];
                int2dz_i1[n1] = int2dz_0[n1];
                int2dz[n + n1 + i * int2d_dim1] = int2dz_0[n1];
            }
#endif
        }

#ifdef __AVX__

        int2dx_2_256 = _mm256_load_pd(&wts[n]);
        int2dy_2_256 = one_256;
        int2dz_2_256 = one_256;

        __m256d d00x_256 = _mm256_load_pd(&d00x[n]);
        int2dx_k1_256 = _mm256_mul_pd(d00x_256, int2dx_2_256);
        _mm256_store_pd(&int2dx[n + int2d_dim2 * int2d_dim1], int2dx_k1_256);

        __m256d d00y_256 = _mm256_load_pd(&d00y[n]);
        int2dy_k1_256 = d00y_256;
        _mm256_store_pd(&int2dy[n + int2d_dim2 * int2d_dim1], int2dy_k1_256);

        __m256d d00z_256 = _mm256_load_pd(&d00z[n]);
        int2dz_k1_256 = d00z_256;
        _mm256_store_pd(&int2dz[n + int2d_dim2 * int2d_dim1], int2dz_k1_256);

#else
#pragma vector aligned
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            weight = wts[n + n1];
            int2dx_2[n1] = weight;
            int2dy_2[n1] = 1.;
            int2dz_2[n1] = 1.;
            int2dx[n + n1 + int2d_dim2 * int2d_dim1] = int2dx_k1[n1] = d00x[n + n1] * weight;
            int2dy[n + n1 + int2d_dim2 * int2d_dim1] = int2dy_k1[n1] = d00y[n + n1];
            int2dz[n + n1 + int2d_dim2 * int2d_dim1] = int2dz_k1[n1] = d00z[n + n1];
        }
#endif

        for (k = 2; k <= shellq; ++k)
        {
            double k1 = k - 1;

#ifdef __AVX__
            __m256d k1_256 = _mm256_broadcast_sd(&k1);
            __m256d b01_256 = _mm256_load_pd(&b01[n]);
            __m256d b1_256 = _mm256_mul_pd(k1_256, b01_256);

            int2dx_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dx_2_256),
                    _mm256_mul_pd(d00x_256, int2dx_k1_256));
            int2dx_2_256 = int2dx_k1_256;
            int2dx_k1_256 = int2dx_0_256;
            _mm256_store_pd(&int2dx[n + k * int2d_dim2 * int2d_dim1], int2dx_0_256);

            int2dy_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dy_2_256),
                    _mm256_mul_pd(d00y_256, int2dy_k1_256));
            int2dy_2_256 = int2dy_k1_256;
            int2dy_k1_256 = int2dy_0_256;
            _mm256_store_pd(&int2dy[n + k * int2d_dim2 * int2d_dim1], int2dy_0_256);

            int2dz_0_256 = _mm256_add_pd(_mm256_mul_pd(b1_256, int2dz_2_256),
                    _mm256_mul_pd(d00z_256, int2dz_k1_256));
            int2dz_2_256 = int2dz_k1_256;
            int2dz_k1_256 = int2dz_0_256;
            _mm256_store_pd(&int2dz[n + k * int2d_dim2 * int2d_dim1], int2dz_0_256);

#else
#pragma vector aligned
            for(n1 = 0; n1 < SIMDW; n1++)
            {
                b1 = k1 * b01[n + n1];
                int2dx_0[n1] = b1 * int2dx_2[n1] + d00x[n + n1] * int2dx_k1[n1];
                int2dx_2[n1] = int2dx_k1[n1];
                int2dx_k1[n1] = int2dx_0[n1];
                int2dx[n + n1 + k * int2d_dim2 * int2d_dim1] = int2dx_0[n1];

                int2dy_0[n1] = b1 * int2dy_2[n1] + d00y[n + n1] * int2dy_k1[n1];
                int2dy_2[n1] = int2dy_k1[n1];
                int2dy_k1[n1] = int2dy_0[n1];
                int2dy[n + n1 + k * int2d_dim2 * int2d_dim1] = int2dy_0[n1];

                int2dz_0[n1] = b1 * int2dz_2[n1] + d00z[n + n1] * int2dz_k1[n1];
                int2dz_2[n1] = int2dz_k1[n1];
                int2dz_k1[n1] = int2dz_0[n1];
                int2dz[n + n1 + k * int2d_dim2 * int2d_dim1] = int2dz_0[n1];
            }
#endif
        }
    }


/*             ...evaluate I=1,SHELLP and K=1,SHELLQ (if any) */
/*                in most economical way. */


    if (shellq <= shellp)
    {
        for (n = 0; n < ngqexq; n+=SIMDW)
        {
#ifdef __AVX__
            __m256d int2dx_00_256, int2dx_10_256, int2dx_20_256, int2dx_11_256;
            __m256d int2dy_00_256, int2dy_10_256, int2dy_20_256, int2dy_11_256;
            __m256d int2dz_00_256, int2dz_10_256, int2dz_20_256, int2dz_11_256;

#else
            double int2dx_00[SIMDW], int2dx_10[SIMDW], int2dx_20[SIMDW], int2dx_11[SIMDW];
            double int2dy_00[SIMDW], int2dy_10[SIMDW], int2dy_20[SIMDW], int2dy_11[SIMDW];
            double int2dz_00[SIMDW], int2dz_10[SIMDW], int2dz_20[SIMDW], int2dz_11[SIMDW];
#endif

            for (k = 1; k <= shellq; ++k)
            {
                int k1 = k - 1;

#ifdef __AVX__
                double k_double = k;
                __m256d k_256 = _mm256_broadcast_sd(&k_double);
                __m256d b00_256 = _mm256_load_pd(&b00[n]);
                __m256d b0_256 = _mm256_mul_pd(k_256, b00_256);

                __m256d int2dx_k1_0_256 = _mm256_load_pd(&int2dx[n + k1 * int2d_dim2 * int2d_dim1]);
                int2dx_20_256 = _mm256_load_pd(&int2dx[n + k * int2d_dim2 * int2d_dim1]);
                __m256d c00x_256 = _mm256_load_pd(&c00x[n]);
                int2dx_10_256 = _mm256_add_pd(_mm256_mul_pd(b0_256, int2dx_k1_0_256),
                                              _mm256_mul_pd(c00x_256, int2dx_20_256));
                _mm256_store_pd(&int2dx[n + (k * int2d_dim2 + 1) * int2d_dim1], int2dx_10_256);

                __m256d int2dy_k1_0_256 = _mm256_load_pd(&int2dy[n + k1 * int2d_dim2 * int2d_dim1]);
                int2dy_20_256 = _mm256_load_pd(&int2dy[n + k * int2d_dim2 * int2d_dim1]);
                __m256d c00y_256 = _mm256_load_pd(&c00y[n]);
                int2dy_10_256 = _mm256_add_pd(_mm256_mul_pd(b0_256, int2dy_k1_0_256),
                                              _mm256_mul_pd(c00y_256, int2dy_20_256));
                _mm256_store_pd(&int2dy[n + (k * int2d_dim2 + 1) * int2d_dim1], int2dy_10_256);

                __m256d int2dz_k1_0_256 = _mm256_load_pd(&int2dz[n + k1 * int2d_dim2 * int2d_dim1]);
                int2dz_20_256 = _mm256_load_pd(&int2dz[n + k * int2d_dim2 * int2d_dim1]);
                __m256d c00z_256 = _mm256_load_pd(&c00z[n]);
                int2dz_10_256 = _mm256_add_pd(_mm256_mul_pd(b0_256, int2dz_k1_0_256),
                                              _mm256_mul_pd(c00z_256, int2dz_20_256));
                _mm256_store_pd(&int2dz[n + (k * int2d_dim2 + 1) * int2d_dim1], int2dz_10_256);

#else
#pragma vector aligned
                for(n1 = 0; n1 < SIMDW; n1++)
                {
                    b0 = k * b00[n + n1];

                    int2dx_10[n1] = b0 * int2dx[n + n1 + k1 * int2d_dim2 * int2d_dim1] +
                        c00x[n + n1] * int2dx[n + n1 + k * int2d_dim2 * int2d_dim1];

                    int2dy_10[n1] = b0 * int2dy[n + n1 + k1 * int2d_dim2 * int2d_dim1] +
                        c00y[n + n1] * int2dy[n + n1 + k * int2d_dim2 * int2d_dim1];

                    int2dz_10[n1] = b0 * int2dz[n + n1 + k1 * int2d_dim2 * int2d_dim1] +
                        c00z[n + n1] * int2dz[n + n1 + k * int2d_dim2 * int2d_dim1];

                    int2dx_20[n1] = int2dx[n + n1 + (k * int2d_dim2) * int2d_dim1];
                    int2dy_20[n1] = int2dy[n + n1 + (k * int2d_dim2) * int2d_dim1];
                    int2dz_20[n1] = int2dz[n + n1 + (k * int2d_dim2) * int2d_dim1];
                }

#pragma vector aligned
                for(n1 = 0; n1 < SIMDW; n1++)
                {
                    int2dx[n + n1 + (k * int2d_dim2 + 1) * int2d_dim1] = int2dx_10[n1];
                    int2dy[n + n1 + (k * int2d_dim2 + 1) * int2d_dim1] = int2dy_10[n1];
                    int2dz[n + n1 + (k * int2d_dim2 + 1) * int2d_dim1] = int2dz_10[n1];
                }
#endif
                for (i = 2; i <= shellp; ++i)
                {
                    int i1 = i - 1;

#ifdef __AVX__

                    double i1_double = i1;
                    __m256d i1_256 = _mm256_broadcast_sd(&i1_double);
                    __m256d b10_256 = _mm256_load_pd(&b10[n]);
                    __m256d b1_256 = _mm256_mul_pd(i1_256, b10_256);

                    int2dx_11_256 = _mm256_load_pd(&int2dx[n + (i1 + k1 * int2d_dim2) * int2d_dim1]);
                    int2dx_00_256 = _mm256_add_pd(_mm256_add_pd(
                                                  _mm256_mul_pd(b0_256, int2dx_11_256),
                                                  _mm256_mul_pd(b1_256, int2dx_20_256)),
                                                  _mm256_mul_pd(c00x_256, int2dx_10_256));
                    int2dx_20_256 = int2dx_10_256;
                    int2dx_10_256 = int2dx_00_256;
                    _mm256_store_pd(&int2dx[n + (i + k * int2d_dim2) * int2d_dim1], int2dx_00_256);

                    int2dy_11_256 = _mm256_load_pd(&int2dy[n + (i1 + k1 * int2d_dim2) * int2d_dim1]);
                    int2dy_00_256 = _mm256_add_pd(_mm256_add_pd(
                                                  _mm256_mul_pd(b0_256, int2dy_11_256),
                                                  _mm256_mul_pd(b1_256, int2dy_20_256)),
                                                  _mm256_mul_pd(c00y_256, int2dy_10_256));
                    int2dy_20_256 = int2dy_10_256;
                    int2dy_10_256 = int2dy_00_256;
                    _mm256_store_pd(&int2dy[n + (i + k * int2d_dim2) * int2d_dim1], int2dy_00_256);

                    int2dz_11_256 = _mm256_load_pd(&int2dz[n + (i1 + k1 * int2d_dim2) * int2d_dim1]);
                    int2dz_00_256 = _mm256_add_pd(_mm256_add_pd(
                                                  _mm256_mul_pd(b0_256, int2dz_11_256),
                                                  _mm256_mul_pd(b1_256, int2dz_20_256)),
                                                  _mm256_mul_pd(c00z_256, int2dz_10_256));
                    int2dz_20_256 = int2dz_10_256;
                    int2dz_10_256 = int2dz_00_256;
                    _mm256_store_pd(&int2dz[n + (i + k * int2d_dim2) * int2d_dim1], int2dz_00_256);

#else
#pragma vector aligned
                    for(n1 = 0; n1 < SIMDW; n1++)
                    {
                        b0 = k * b00[n + n1];
                        b1 = i1 * b10[n + n1];
                        int2dx_11[n1] = int2dx[n + n1 + (i1 + k1 * int2d_dim2) * int2d_dim1];
                        int2dx_00[n1] = b0 * int2dx_11[n1] + b1 * int2dx_20[n1] + c00x[n + n1] * int2dx_10[n1];
                        int2dx_20[n1] = int2dx_10[n1];
                        int2dx_10[n1] = int2dx_00[n1];

                        int2dy_11[n1] = int2dy[n + n1 + (i1 + k1 * int2d_dim2) * int2d_dim1];
                        int2dy_00[n1] = b0 * int2dy_11[n1] + b1 * int2dy_20[n1] + c00y[n + n1] * int2dy_10[n1];
                        int2dy_20[n1] = int2dy_10[n1];
                        int2dy_10[n1] = int2dy_00[n1];

                        int2dz_11[n1] = int2dz[n + n1 + (i1 + k1 * int2d_dim2) * int2d_dim1];
                        int2dz_00[n1] = b0 * int2dz_11[n1] + b1 * int2dz_20[n1] + c00z[n + n1] * int2dz_10[n1];
                        int2dz_20[n1] = int2dz_10[n1];
                        int2dz_10[n1] = int2dz_00[n1];
                    }
#pragma vector aligned
                    for(n1 = 0; n1 < SIMDW; n1++)
                    {
                        int2dx[n + n1 + (i + k * int2d_dim2) * int2d_dim1] = int2dx_00[n1];
                        int2dy[n + n1 + (i + k * int2d_dim2) * int2d_dim1] = int2dy_00[n1];
                        int2dz[n + n1 + (i + k * int2d_dim2) * int2d_dim1] = int2dz_00[n1];
                    }
#endif
                }
            }
        }
    }
    else
    {
        for (n = 0; n < ngqexq; n+=SIMDW)
        {
#ifdef __AVX__
            __m256d int2dx_00_256, int2dx_01_256, int2dx_02_256, int2dx_11_256;
            __m256d int2dy_00_256, int2dy_01_256, int2dy_02_256, int2dy_11_256;
            __m256d int2dz_00_256, int2dz_01_256, int2dz_02_256, int2dz_11_256;

#else
            double int2dx_00[SIMDW], int2dx_01[SIMDW], int2dx_02[SIMDW], int2dx_11[SIMDW];
            double int2dy_00[SIMDW], int2dy_01[SIMDW], int2dy_02[SIMDW], int2dy_11[SIMDW];
            double int2dz_00[SIMDW], int2dz_01[SIMDW], int2dz_02[SIMDW], int2dz_11[SIMDW];
#endif

            for (i = 1; i <= shellp; ++i)
            {
                int i1 = i - 1;

#ifdef __AVX__

                double i_double = i;
                __m256d i_256 = _mm256_broadcast_sd(&i_double);
                __m256d b00_256 = _mm256_load_pd(&b00[n]);
                __m256d b0_256 = _mm256_mul_pd(i_256, b00_256);

                __m256d int2dx_0_i1_256 = _mm256_load_pd(&int2dx[n + i1 * int2d_dim1]);
                int2dx_02_256 = _mm256_load_pd(&int2dx[n + i * int2d_dim1]);
                __m256d d00x_256 = _mm256_load_pd(&d00x[n]);
                int2dx_01_256 = _mm256_add_pd(_mm256_mul_pd(b0_256, int2dx_0_i1_256),
                                              _mm256_mul_pd(d00x_256, int2dx_02_256));
                _mm256_store_pd(&int2dx[n + (i + int2d_dim2) * int2d_dim1], int2dx_01_256);

                __m256d int2dy_0_i1_256 = _mm256_load_pd(&int2dy[n + i1 * int2d_dim1]);
                int2dy_02_256 = _mm256_load_pd(&int2dy[n + i * int2d_dim1]);
                __m256d d00y_256 = _mm256_load_pd(&d00y[n]);
                int2dy_01_256 = _mm256_add_pd(_mm256_mul_pd(b0_256, int2dy_0_i1_256),
                                              _mm256_mul_pd(d00y_256, int2dy_02_256));
                _mm256_store_pd(&int2dy[n + (i + int2d_dim2) * int2d_dim1], int2dy_01_256);

                __m256d int2dz_0_i1_256 = _mm256_load_pd(&int2dz[n + i1 * int2d_dim1]);
                int2dz_02_256 = _mm256_load_pd(&int2dz[n + i * int2d_dim1]);
                __m256d d00z_256 = _mm256_load_pd(&d00z[n]);
                int2dz_01_256 = _mm256_add_pd(_mm256_mul_pd(b0_256, int2dz_0_i1_256),
                                              _mm256_mul_pd(d00z_256, int2dz_02_256));
                _mm256_store_pd(&int2dz[n + (i + int2d_dim2) * int2d_dim1], int2dz_01_256);

#else
#pragma vector aligned
                for(n1 = 0; n1 < SIMDW; n1++)
                {
                    b0 = i * b00[n + n1];
                    int2dx_01[n1] = b0 * int2dx[n + n1 + i1 * int2d_dim1]
                        + d00x[n + n1] * int2dx[n + n1 + i * int2d_dim1];

                    int2dy_01[n1] = b0 * int2dy[n + n1 + i1 * int2d_dim1]
                        + d00y[n + n1] * int2dy[n + n1 + i * int2d_dim1];

                    int2dz_01[n1] = b0 * int2dz[n + n1 + i1 * int2d_dim1]
                        + d00z[n + n1] * int2dz[n + n1 + i * int2d_dim1];

                    int2dx_02[n1] = int2dx[n + n1 + (i) * int2d_dim1];
                    int2dy_02[n1] = int2dy[n + n1 + (i) * int2d_dim1];
                    int2dz_02[n1] = int2dz[n + n1 + (i) * int2d_dim1];
                }
#pragma vector aligned
                for(n1 = 0; n1 < SIMDW; n1++)
                {
                    int2dx[n + n1 + (i + int2d_dim2) * int2d_dim1] = int2dx_01[n1];
                    int2dy[n + n1 + (i + int2d_dim2) * int2d_dim1] = int2dy_01[n1];
                    int2dz[n + n1 + (i + int2d_dim2) * int2d_dim1] = int2dz_01[n1];
                }
#endif

                for (k = 2; k <= shellq; ++k)
                {
                    int k1 = k - 1;

#ifdef __AVX__

                    double k1_double = k1;
                    __m256d k1_256 = _mm256_broadcast_sd(&k1_double);
                    __m256d b01_256 = _mm256_load_pd(&b01[n]);
                    __m256d b1_256 = _mm256_mul_pd(k1_256, b01_256);

                    int2dx_11_256 = _mm256_load_pd(&int2dx[n + (i1 + k1 * int2d_dim2) * int2d_dim1]);
                    int2dx_00_256 = _mm256_add_pd(_mm256_add_pd(
                                                  _mm256_mul_pd(b0_256, int2dx_11_256),
                                                  _mm256_mul_pd(b1_256, int2dx_02_256)),
                                                  _mm256_mul_pd(d00x_256, int2dx_01_256));
                    int2dx_02_256 = int2dx_01_256;
                    int2dx_01_256 = int2dx_00_256;
                    _mm256_store_pd(&int2dx[n + (i + k * int2d_dim2) * int2d_dim1], int2dx_00_256);

                    int2dy_11_256 = _mm256_load_pd(&int2dy[n + (i1 + k1 * int2d_dim2) * int2d_dim1]);
                    int2dy_00_256 = _mm256_add_pd(_mm256_add_pd(
                                                  _mm256_mul_pd(b0_256, int2dy_11_256),
                                                  _mm256_mul_pd(b1_256, int2dy_02_256)),
                                                  _mm256_mul_pd(d00y_256, int2dy_01_256));
                    int2dy_02_256 = int2dy_01_256;
                    int2dy_01_256 = int2dy_00_256;
                    _mm256_store_pd(&int2dy[n + (i + k * int2d_dim2) * int2d_dim1], int2dy_00_256);

                    int2dz_11_256 = _mm256_load_pd(&int2dz[n + (i1 + k1 * int2d_dim2) * int2d_dim1]);
                    int2dz_00_256 = _mm256_add_pd(_mm256_add_pd(
                                                  _mm256_mul_pd(b0_256, int2dz_11_256),
                                                  _mm256_mul_pd(b1_256, int2dz_02_256)),
                                                  _mm256_mul_pd(d00z_256, int2dz_01_256));
                    int2dz_02_256 = int2dz_01_256;
                    int2dz_01_256 = int2dz_00_256;
                    _mm256_store_pd(&int2dz[n + (i + k * int2d_dim2) * int2d_dim1], int2dz_00_256);

#else
#pragma vector aligned
                    for(n1 = 0; n1 < SIMDW; n1++)
                    {
                        b0 = i * b00[n + n1];
                        b1 = k1 * b01[n + n1];

                        int2dx_11[n1] = int2dx[n + n1 + (i1 + k1 * int2d_dim2) * int2d_dim1];
                        int2dx_00[n1] = b0 * int2dx_11[n1] + b1 * int2dx_02[n1] + d00x[n + n1] * int2dx_01[n1];
                        int2dx_02[n1] = int2dx_01[n1];
                        int2dx_01[n1] = int2dx_00[n1];

                        int2dy_11[n1] = int2dy[n + n1 + (i1 + k1 * int2d_dim2) * int2d_dim1];
                        int2dy_00[n1] = b0 * int2dy_11[n1] + b1 * int2dy_02[n1] + d00y[n + n1] * int2dy_01[n1];
                        int2dy_02[n1] = int2dy_01[n1];
                        int2dy_01[n1] = int2dy_00[n1];

                        int2dz_11[n1] = int2dz[n + n1 + (i1 + k1 * int2d_dim2) * int2d_dim1];
                        int2dz_00[n1] = b0 * int2dz_11[n1] + b1 * int2dz_02[n1] + d00z[n + n1] * int2dz_01[n1];
                        int2dz_02[n1] = int2dz_01[n1];
                        int2dz_01[n1] = int2dz_00[n1];
                    }
#pragma vector aligned
                    for(n1 = 0; n1 < SIMDW; n1++)
                    {
                        int2dx[n + n1 + (i + k * int2d_dim2) * int2d_dim1] = int2dx_00[n1];
                        int2dy[n + n1 + (i + k * int2d_dim2) * int2d_dim1] = int2dy_00[n1];
                        int2dz[n + n1 + (i + k * int2d_dim2) * int2d_dim1] = int2dz_00[n1];
                    }
#endif
                }
            }
        }
    }


/*             ...ready! */


    return 0;
}                               /* erd__2d_pq_integrals__ */


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
