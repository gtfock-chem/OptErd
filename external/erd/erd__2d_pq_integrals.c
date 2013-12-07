/* erd__2d_pq_integrals.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

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
erd__2d_pq_integrals_ (int * shellp, int * shellq,
                        int * ngqexq, double * wts, double * b00,
                        double * b01, double * b10, double * c00x,
                        double * c00y, double * c00z,
                        double * d00x, double * d00y,
                        double * d00z, int * case2d,
                        double * int2dx, double * int2dy,
                        double * int2dz)
{
    /* System generated locals */
    int int2dx_dim1, int2dx_dim2, int2dx_offset, int2dy_dim1, int2dy_dim2,
        int2dy_offset, int2dz_dim1, int2dz_dim2, int2dz_offset, i__1,
        i__2, i__3;

    /* Local variables */
    static double f;
    static int i__, k, n;
    static double b0, b1, f1, f2;
    static int i1, i2, k1, k2;
    static double weight;

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
    int2dz_dim1 = *ngqexq - 1 + 1;
    int2dz_dim2 = *shellp - 0 + 1;
    int2dz_offset = 1 + int2dz_dim1 * (0 + int2dz_dim2 * 0);
    int2dz -= int2dz_offset;
    int2dy_dim1 = *ngqexq - 1 + 1;
    int2dy_dim2 = *shellp - 0 + 1;
    int2dy_offset = 1 + int2dy_dim1 * (0 + int2dy_dim2 * 0);
    int2dy -= int2dy_offset;
    int2dx_dim1 = *ngqexq - 1 + 1;
    int2dx_dim2 = *shellp - 0 + 1;
    int2dx_offset = 1 + int2dx_dim1 * (0 + int2dx_dim2 * 0);
    int2dx -= int2dx_offset;
    --d00z;
    --d00y;
    --d00x;
    --c00z;
    --c00y;
    --c00x;
    --b10;
    --b01;
    --b00;
    --wts;

    /* Function Body */
    switch (*case2d)
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
    i__1 = *ngqexq;
    for (n = 1; n <= i__1; ++n)
    {
        int2dx[n] = wts[n];
        int2dy[n] = 1.;
        int2dz[n] = 1.;
/* L100: */
    }
    return 0;


/*             ...the cases P = s-shell and Q >= p-shell. */
/*                Evaluate I=0 and K=0,1. */


  L2:
    i__1 = *ngqexq;
    for (n = 1; n <= i__1; ++n)
    {
        weight = wts[n];
        int2dx[n] = weight;
        int2dx[n + int2dx_dim2 * int2dx_dim1] = d00x[n] * weight;
        int2dy[n] = 1.;
        int2dy[n + int2dy_dim2 * int2dy_dim1] = d00y[n];
        int2dz[n] = 1.;
        int2dz[n + int2dz_dim2 * int2dz_dim1] = d00z[n];
/* L200: */
    }


/*             ...evaluate I=0 and K=2,SHELLQ (if any). */


    f = 1.;
    i__1 = *shellq;
    for (k = 2; k <= i__1; ++k)
    {
        k1 = k - 1;
        k2 = k - 2;
        i__2 = *ngqexq;
        for (n = 1; n <= i__2; ++n)
        {
            b1 = f * b01[n];
            int2dx[n + k * int2dx_dim2 * int2dx_dim1] = b1 * int2dx[n + k2 *
                                                                    int2dx_dim2
                                                                    *
                                                                    int2dx_dim1]
                + d00x[n] * int2dx[n + k1 * int2dx_dim2 * int2dx_dim1];
            int2dy[n + k * int2dy_dim2 * int2dy_dim1] =
                b1 * int2dy[n + k2 * int2dy_dim2 * int2dy_dim1] +
                d00y[n] * int2dy[n + k1 * int2dy_dim2 * int2dy_dim1];
            int2dz[n + k * int2dz_dim2 * int2dz_dim1] =
                b1 * int2dz[n + k2 * int2dz_dim2 * int2dz_dim1] +
                d00z[n] * int2dz[n + k1 * int2dz_dim2 * int2dz_dim1];
/* L212: */
        }
        f += 1.;
/* L210: */
    }
    return 0;


/*             ...the cases P >= p-shell and Q = s-shell. */
/*                Evaluate I=0,1 and K=0. */


  L3:
    i__1 = *ngqexq;
    for (n = 1; n <= i__1; ++n)
    {
        weight = wts[n];
        int2dx[n] = weight;
        int2dx[n + int2dx_dim1] = c00x[n] * weight;
        int2dy[n] = 1.;
        int2dy[n + int2dy_dim1] = c00y[n];
        int2dz[n] = 1.;
        int2dz[n + int2dz_dim1] = c00z[n];
/* L300: */
    }


/*             ...evaluate I=2,SHELLP (if any) and K=0. */


    f = 1.;
    i__1 = *shellp;
    for (i__ = 2; i__ <= i__1; ++i__)
    {
        i1 = i__ - 1;
        i2 = i__ - 2;
        i__2 = *ngqexq;
        for (n = 1; n <= i__2; ++n)
        {
            b1 = f * b10[n];
            int2dx[n + i__ * int2dx_dim1] = b1 * int2dx[n + i2 * int2dx_dim1]
                + c00x[n] * int2dx[n + i1 * int2dx_dim1];
            int2dy[n + i__ * int2dy_dim1] = b1 * int2dy[n + i2 * int2dy_dim1]
                + c00y[n] * int2dy[n + i1 * int2dy_dim1];
            int2dz[n + i__ * int2dz_dim1] = b1 * int2dz[n + i2 * int2dz_dim1]
                + c00z[n] * int2dz[n + i1 * int2dz_dim1];
/* L312: */
        }
        f += 1.;
/* L310: */
    }
    return 0;


/*             ...the cases P >= p-shell and Q >= p-shell. */
/*                Evaluate I=0,SHELLP       I=0 */
/*                         K=0        and   K=0,SHELLQ */


  L4:
    i__1 = *ngqexq;
    for (n = 1; n <= i__1; ++n)
    {
        weight = wts[n];
        int2dx[n] = weight;
        int2dx[n + int2dx_dim1] = c00x[n] * weight;
        int2dx[n + int2dx_dim2 * int2dx_dim1] = d00x[n] * weight;
        int2dy[n] = 1.;
        int2dy[n + int2dy_dim1] = c00y[n];
        int2dy[n + int2dy_dim2 * int2dy_dim1] = d00y[n];
        int2dz[n] = 1.;
        int2dz[n + int2dz_dim1] = c00z[n];
        int2dz[n + int2dz_dim2 * int2dz_dim1] = d00z[n];
/* L400: */
    }
    f = 1.;
    i__1 = *shellp;
    for (i__ = 2; i__ <= i__1; ++i__)
    {
        i1 = i__ - 1;
        i2 = i__ - 2;
        i__2 = *ngqexq;
        for (n = 1; n <= i__2; ++n)
        {
            b1 = f * b10[n];
            int2dx[n + i__ * int2dx_dim1] = b1 * int2dx[n + i2 * int2dx_dim1]
                + c00x[n] * int2dx[n + i1 * int2dx_dim1];
            int2dy[n + i__ * int2dy_dim1] = b1 * int2dy[n + i2 * int2dy_dim1]
                + c00y[n] * int2dy[n + i1 * int2dy_dim1];
            int2dz[n + i__ * int2dz_dim1] = b1 * int2dz[n + i2 * int2dz_dim1]
                + c00z[n] * int2dz[n + i1 * int2dz_dim1];
/* L412: */
        }
        f += 1.;
/* L410: */
    }
    f = 1.;
    i__1 = *shellq;
    for (k = 2; k <= i__1; ++k)
    {
        k1 = k - 1;
        k2 = k - 2;
        i__2 = *ngqexq;
        for (n = 1; n <= i__2; ++n)
        {
            b1 = f * b01[n];
            int2dx[n + k * int2dx_dim2 * int2dx_dim1] = b1 * int2dx[n + k2 *
                                                                    int2dx_dim2
                                                                    *
                                                                    int2dx_dim1]
                + d00x[n] * int2dx[n + k1 * int2dx_dim2 * int2dx_dim1];
            int2dy[n + k * int2dy_dim2 * int2dy_dim1] =
                b1 * int2dy[n + k2 * int2dy_dim2 * int2dy_dim1] +
                d00y[n] * int2dy[n + k1 * int2dy_dim2 * int2dy_dim1];
            int2dz[n + k * int2dz_dim2 * int2dz_dim1] =
                b1 * int2dz[n + k2 * int2dz_dim2 * int2dz_dim1] +
                d00z[n] * int2dz[n + k1 * int2dz_dim2 * int2dz_dim1];
/* L416: */
        }
        f += 1.;
/* L414: */
    }


/*             ...evaluate I=1,SHELLP and K=1,SHELLQ (if any) */
/*                in most economical way. */


    if (*shellq <= *shellp)
    {
        f1 = 1.;
        i__1 = *shellq;
        for (k = 1; k <= i__1; ++k)
        {
            k1 = k - 1;
            i__2 = *ngqexq;
            for (n = 1; n <= i__2; ++n)
            {
                b0 = f1 * b00[n];
                int2dx[n + (k * int2dx_dim2 + 1) * int2dx_dim1] =
                    b0 * int2dx[n + k1 * int2dx_dim2 * int2dx_dim1] +
                    c00x[n] * int2dx[n + k * int2dx_dim2 * int2dx_dim1];
                int2dy[n + (k * int2dy_dim2 + 1) * int2dy_dim1] =
                    b0 * int2dy[n + k1 * int2dy_dim2 * int2dy_dim1] +
                    c00y[n] * int2dy[n + k * int2dy_dim2 * int2dy_dim1];
                int2dz[n + (k * int2dz_dim2 + 1) * int2dz_dim1] =
                    b0 * int2dz[n + k1 * int2dz_dim2 * int2dz_dim1] +
                    c00z[n] * int2dz[n + k * int2dz_dim2 * int2dz_dim1];
/* L421: */
            }
            f2 = 1.;
            i__2 = *shellp;
            for (i__ = 2; i__ <= i__2; ++i__)
            {
                i1 = i__ - 1;
                i2 = i__ - 2;
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    b0 = f1 * b00[n];
                    b1 = f2 * b10[n];
                    int2dx[n + (i__ + k * int2dx_dim2) * int2dx_dim1] = b0 *
                        int2dx[n + (i1 + k1 * int2dx_dim2) * int2dx_dim1]
                        + b1 * int2dx[n + (i2 + k * int2dx_dim2) *
                                      int2dx_dim1] + c00x[n] * int2dx[n +
                                                                      (i1 +
                                                                       k *
                                                                       int2dx_dim2)
                                                                      *
                                                                      int2dx_dim1];
                    int2dy[n + (i__ + k * int2dy_dim2) * int2dy_dim1] =
                        b0 * int2dy[n +
                                    (i1 + k1 * int2dy_dim2) * int2dy_dim1] +
                        b1 * int2dy[n +
                                    (i2 + k * int2dy_dim2) * int2dy_dim1] +
                        c00y[n] * int2dy[n +
                                         (i1 +
                                          k * int2dy_dim2) * int2dy_dim1];
                    int2dz[n + (i__ + k * int2dz_dim2) * int2dz_dim1] =
                        b0 * int2dz[n +
                                    (i1 + k1 * int2dz_dim2) * int2dz_dim1] +
                        b1 * int2dz[n +
                                    (i2 + k * int2dz_dim2) * int2dz_dim1] +
                        c00z[n] * int2dz[n +
                                         (i1 +
                                          k * int2dz_dim2) * int2dz_dim1];
/* L423: */
                }
                f2 += 1.;
/* L422: */
            }
            f1 += 1.;
/* L420: */
        }
    }
    else
    {
        f1 = 1.;
        i__1 = *shellp;
        for (i__ = 1; i__ <= i__1; ++i__)
        {
            i1 = i__ - 1;
            i__2 = *ngqexq;
            for (n = 1; n <= i__2; ++n)
            {
                b0 = f1 * b00[n];
                int2dx[n + (i__ + int2dx_dim2) * int2dx_dim1] = b0 * int2dx[n
                                                                            +
                                                                            i1
                                                                            *
                                                                            int2dx_dim1]
                    + d00x[n] * int2dx[n + i__ * int2dx_dim1];
                int2dy[n + (i__ + int2dy_dim2) * int2dy_dim1] =
                    b0 * int2dy[n + i1 * int2dy_dim1] + d00y[n] * int2dy[n +
                                                                         i__ *
                                                                         int2dy_dim1];
                int2dz[n + (i__ + int2dz_dim2) * int2dz_dim1] =
                    b0 * int2dz[n + i1 * int2dz_dim1] + d00z[n] * int2dz[n +
                                                                         i__ *
                                                                         int2dz_dim1];
/* L431: */
            }
            f2 = 1.;
            i__2 = *shellq;
            for (k = 2; k <= i__2; ++k)
            {
                k1 = k - 1;
                k2 = k - 2;
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    b0 = f1 * b00[n];
                    b1 = f2 * b01[n];
                    int2dx[n + (i__ + k * int2dx_dim2) * int2dx_dim1] = b0 *
                        int2dx[n + (i1 + k1 * int2dx_dim2) * int2dx_dim1]
                        + b1 * int2dx[n + (i__ + k2 * int2dx_dim2) *
                                      int2dx_dim1] + d00x[n] * int2dx[n +
                                                                      (i__ +
                                                                       k1 *
                                                                       int2dx_dim2)
                                                                      *
                                                                      int2dx_dim1];
                    int2dy[n + (i__ + k * int2dy_dim2) * int2dy_dim1] =
                        b0 * int2dy[n +
                                    (i1 + k1 * int2dy_dim2) * int2dy_dim1] +
                        b1 * int2dy[n +
                                    (i__ + k2 * int2dy_dim2) * int2dy_dim1] +
                        d00y[n] * int2dy[n +
                                         (i__ +
                                          k1 * int2dy_dim2) * int2dy_dim1];
                    int2dz[n + (i__ + k * int2dz_dim2) * int2dz_dim1] =
                        b0 * int2dz[n +
                                    (i1 + k1 * int2dz_dim2) * int2dz_dim1] +
                        b1 * int2dz[n +
                                    (i__ + k2 * int2dz_dim2) * int2dz_dim1] +
                        d00z[n] * int2dz[n +
                                         (i__ +
                                          k1 * int2dz_dim2) * int2dz_dim1];
/* L433: */
                }
                f2 += 1.;
/* L432: */
            }
            f1 += 1.;
/* L430: */
        }
    }


/*             ...ready! */


    return 0;
}                               /* erd__2d_pq_integrals__ */
