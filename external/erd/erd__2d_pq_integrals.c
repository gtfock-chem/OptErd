#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


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
/* ------------------------------------------------------------------------ */
int erd__2d_pq_integrals (int shellp, int shellq, int ngqexq,
                          double *wts, double *b00, double *b01, double *b10,
                          double *c00x, double *c00y, double *c00z,
                          double *d00x, double *d00y, double *d00z,
                          int case2d, double *int2dx,
                          double *int2dy, double *int2dz)
{
    double f;
    int i, k, n;
    double b0, b1, f1, f2;
    int i1, i2, k1, k2;
    double weight;

/*             ...jump according to the 4 different cases that can arise: */
/*                  P-shell = s- or higher angular momentum */
/*                  Q-shell = s- or higher angular momentum */
/*                each leading to simplifications in the VRR formulas. */
/*                The case present has been evaluated outside this */
/*                routine and is transmitted via argument. */
    goto L4;
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
    for (n = 0; n < ngqexq; ++n)
    {
        int2dx[n] = wts[n];
        int2dy[n] = 1.;
        int2dz[n] = 1.;
    }
    return 0;


/*             ...the cases P = s-shell and Q >= p-shell. */
/*                Evaluate I=0 and K=0,1. */
  L2:
    for (n = 0; n < ngqexq; ++n)
    {
        weight = wts[n];
        int2dx[n] = weight;
        int2dx[n + (shellp + 1) * ngqexq] = d00x[n] * weight;
        int2dy[n] = 1.;
        int2dy[n + (shellp + 1) * ngqexq] = d00y[n];
        int2dz[n] = 1.;
        int2dz[n + (shellp + 1) * ngqexq] = d00z[n];
    }


/*             ...evaluate I=0 and K=2,SHELLQ (if any). */
    f = 1.;
    for (k = 2; k <= shellq; ++k)
    {
        k1 = k - 1;
        k2 = k - 2;
        for (n = 0; n < ngqexq; ++n)
        {
            b1 = f * b01[n];
            int2dx[n + k * (shellp + 1) * ngqexq] =
                b1 * int2dx[n + k2 * (shellp + 1) * ngqexq]
                + d00x[n] * int2dx[n + k1 * (shellp + 1) * ngqexq];
            int2dy[n + k * (shellp + 1) * ngqexq] =
                b1 * int2dy[n + k2 * (shellp + 1) * ngqexq] +
                d00y[n] * int2dy[n + k1 * (shellp + 1) * ngqexq];
            int2dz[n + k * (shellp + 1) * ngqexq] =
                b1 * int2dz[n + k2 * (shellp + 1) * ngqexq] +
                d00z[n] * int2dz[n + k1 * (shellp + 1) * ngqexq];
        }
        f += 1.;
    }
    return 0;


/*             ...the cases P >= p-shell and Q = s-shell. */
/*                Evaluate I=0,1 and K=0. */
  L3:
    for (n = 0; n < ngqexq; ++n)
    {
        weight = wts[n];
        int2dx[n] = weight;
        int2dx[n + ngqexq] = c00x[n] * weight;
        int2dy[n] = 1.;
        int2dy[n + ngqexq] = c00y[n];
        int2dz[n] = 1.;
        int2dz[n + ngqexq] = c00z[n];
    }


/*             ...evaluate I=2,SHELLP (if any) and K=0. */
    f = 1.;
    for (i = 2; i <= shellp; ++i)
    {
        i1 = i - 1;
        i2 = i - 2;
        for (n = 0; n < ngqexq; ++n)
        {
            b1 = f * b10[n];
            int2dx[n + i * ngqexq] = b1 * int2dx[n + i2 * ngqexq]
                + c00x[n] * int2dx[n + i1 * ngqexq];
            int2dy[n + i * ngqexq] = b1 * int2dy[n + i2 * ngqexq]
                + c00y[n] * int2dy[n + i1 * ngqexq];
            int2dz[n + i * ngqexq] = b1 * int2dz[n + i2 * ngqexq]
                + c00z[n] * int2dz[n + i1 * ngqexq];
        }
        f += 1.;
    }
    return 0;


/*             ...the cases P >= p-shell and Q >= p-shell. */
/*                Evaluate I=0,SHELLP       I=0 */
/*                         K=0        and   K=0,SHELLQ */
  L4:
    for (n = 0; n < ngqexq; ++n)
    {
        weight = wts[n];
        int2dx[n] = weight;
        int2dx[n + ngqexq] = c00x[n] * weight;
        int2dx[n + (shellp + 1) * ngqexq] = d00x[n] * weight;
        int2dy[n] = 1.;
        int2dy[n + ngqexq] = c00y[n];
        int2dy[n + (shellp + 1) * ngqexq] = d00y[n];
        int2dz[n] = 1.;
        int2dz[n + ngqexq] = c00z[n];
        int2dz[n + (shellp + 1) * ngqexq] = d00z[n];
    }
    f = 1.;
    for (i = 2; i <= shellp; ++i)
    {
        i1 = i - 1;
        i2 = i - 2;
        for (n = 0; n < ngqexq; ++n)
        {
            b1 = f * b10[n];
            int2dx[n + i * ngqexq] = b1 * int2dx[n + i2 * ngqexq]
                + c00x[n] * int2dx[n + i1 * ngqexq];
            int2dy[n + i * ngqexq] = b1 * int2dy[n + i2 * ngqexq]
                + c00y[n] * int2dy[n + i1 * ngqexq];
            int2dz[n + i * ngqexq] = b1 * int2dz[n + i2 * ngqexq]
                + c00z[n] * int2dz[n + i1 * ngqexq];
        }
        f += 1.;
    }
    f = 1.;
    for (k = 2; k <= shellq; ++k)
    {
        k1 = k - 1;
        k2 = k - 2;
        for (n = 0; n < ngqexq; ++n)
        {
            b1 = f * b01[n];
            int2dx[n + k * (shellp + 1) * ngqexq] = b1 *
                int2dx[n + k2 * (shellp + 1) * ngqexq]
                + d00x[n] * int2dx[n + k1 * (shellp + 1) * ngqexq];
            int2dy[n + k * (shellp + 1) * ngqexq] =
                b1 * int2dy[n + k2 * (shellp + 1) * ngqexq] +
                d00y[n] * int2dy[n + k1 * (shellp + 1) * ngqexq];
            int2dz[n + k * (shellp + 1) * ngqexq] =
                b1 * int2dz[n + k2 * (shellp + 1) * ngqexq] +
                d00z[n] * int2dz[n + k1 * (shellp + 1) * ngqexq];
        }
        f += 1.;
    }


/*             ...evaluate I=1,SHELLP and K=1,SHELLQ (if any) */
/*                in most economical way. */
    if (shellq <= shellp)
    {
        f1 = 1.;
        for (k = 1; k <= shellq; ++k)
        {
            k1 = k - 1;
            for (n = 0; n < ngqexq; ++n)
            {
                b0 = f1 * b00[n];
                int2dx[n + (k * (shellp + 1) + 1) * ngqexq] =
                    b0 * int2dx[n + k1 * (shellp + 1) * ngqexq] +
                    c00x[n] * int2dx[n + k * (shellp + 1) * ngqexq];
                int2dy[n + (k * (shellp + 1) + 1) * ngqexq] =
                    b0 * int2dy[n + k1 * (shellp + 1) * ngqexq] +
                    c00y[n] * int2dy[n + k * (shellp + 1) * ngqexq];
                int2dz[n + (k * (shellp + 1) + 1) * ngqexq] =
                    b0 * int2dz[n + k1 * (shellp + 1) * ngqexq] +
                    c00z[n] * int2dz[n + k * (shellp + 1) * ngqexq];
            }
            f2 = 1.;
            for (i = 2; i <= shellp; ++i)
            {
                i1 = i - 1;
                i2 = i - 2;
                for (n = 0; n < ngqexq; ++n)
                {
                    b0 = f1 * b00[n];
                    b1 = f2 * b10[n];
                    int2dx[n + (i + k * (shellp + 1)) * ngqexq] =
                        b0 * int2dx[n + (i1 + k1 * (shellp + 1)) * ngqexq] +
                        b1 * int2dx[n + (i2 + k * (shellp + 1)) * ngqexq] +
                        c00x[n] * int2dx[n + (i1 + k * (shellp + 1)) * ngqexq];
                    int2dy[n + (i + k * (shellp + 1)) * ngqexq] =
                        b0 * int2dy[n + (i1 + k1 * (shellp + 1)) * ngqexq] +
                        b1 * int2dy[n + (i2 + k * (shellp + 1)) * ngqexq] +
                        c00y[n] * int2dy[n + (i1 + k * (shellp + 1)) * ngqexq];
                    int2dz[n + (i + k * (shellp + 1)) * ngqexq] =
                        b0 * int2dz[n + (i1 + k1 * (shellp + 1)) * ngqexq] +
                        b1 * int2dz[n + (i2 + k * (shellp + 1)) * ngqexq] +
                        c00z[n] * int2dz[n + (i1 + k * (shellp + 1)) * ngqexq];
                }
                f2 += 1.;
            }
            f1 += 1.;
        }
    }
    else
    {
        f1 = 1.;
        for (i = 1; i <= shellp; ++i)
        {
            i1 = i - 1;
            for (n = 0; n < ngqexq; ++n)
            {
                b0 = f1 * b00[n];
                int2dx[n + (i + (shellp + 1)) * ngqexq] =
                    b0 * int2dx[n + i1 * ngqexq]
                    + d00x[n] * int2dx[n + i * ngqexq];
                int2dy[n + (i + (shellp + 1)) * ngqexq] =
                    b0 * int2dy[n + i1 * ngqexq]
                    + d00y[n] * int2dy[n + i * ngqexq];
                int2dz[n + (i + (shellp + 1)) * ngqexq] =
                    b0 * int2dz[n + i1 * ngqexq]
                    + d00z[n] * int2dz[n + i * ngqexq];
            }
            f2 = 1.;
            for (k = 2; k <= shellq; ++k)
            {
                k1 = k - 1;
                k2 = k - 2;
                for (n = 0; n < ngqexq; ++n)
                {
                    b0 = f1 * b00[n];
                    b1 = f2 * b01[n];
                    int2dx[n + (i + k * (shellp + 1)) * ngqexq] =
                        b0 * int2dx[n + (i1 + k1 * (shellp + 1)) * ngqexq] +
                        b1 * int2dx[n + (i + k2 * (shellp + 1)) * ngqexq] +
                        d00x[n] * int2dx[n + (i + k1 * (shellp + 1)) * ngqexq];
                    int2dy[n + (i + k * (shellp + 1)) * ngqexq] =
                        b0 * int2dy[n + (i1 + k1 * (shellp + 1)) * ngqexq] +
                        b1 * int2dy[n + (i + k2 * (shellp + 1)) * ngqexq] +
                        d00y[n] * int2dy[n + (i + k1 * (shellp + 1)) * ngqexq];
                    int2dz[n + (i + k * (shellp + 1)) * ngqexq] =
                        b0 * int2dz[n + (i1 + k1 * (shellp + 1)) * ngqexq] +
                        b1 * int2dz[n + (i + k2 * (shellp + 1)) * ngqexq] +
                        d00z[n] * int2dz[n + (i + k1 * (shellp + 1)) * ngqexq];
                }
                f2 += 1.;
            }
            f1 += 1.;
        }
    }


    return 0;
}
