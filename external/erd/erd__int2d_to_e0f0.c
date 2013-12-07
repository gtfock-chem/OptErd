/* erd__int2d_to_e0f0.f -- translated by f2c (version 20100827).
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
#include "erd.h"

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
erd__int2d_to_e0f0_ (int * shella, int * shellp,
                      int * shellc, int * shellq, int * ngqp,
                      int * nexq, int * ngqexq, int * nxyzet,
                      int * nxyzft, int * nxyzp, int * nxyzq,
                      double * int2dx, double * int2dy,
                      double * int2dz, double * temp1,
                      double * temp2, double * scale,
                      double * batch)
{
    /* System generated locals */
    int batch_dim1, batch_dim2, batch_offset, int2dx_dim1, int2dx_dim2,
        int2dx_offset, int2dy_dim1, int2dy_dim2, int2dy_offset,
        int2dz_dim1, int2dz_dim2, int2dz_offset, i__1, i__2, i__3, i__4,
        i__5, i__6, i__7, i__8;

    int i__, j, k, m, n, se, sf, xe, ye, ze, xf, yf, zf, xep, xfp;
    double sum;
    int xye, xyf, xyep, xyfp, seend, sfend, yeend, yfend, xemax,
        xfmax, nxyze, nxyzf;

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__INT2D_TO_E0F0 */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine assembles the set of batches of cartesian */
/*                eris [E0|F0] , E = A to P, F = C to Q, adding up all */
/*                the contributions from all the 2D PQ integrals. */

/*                The routine uses the reduced Rys multiplication scheme */
/*                as suggested in R.Lindh, U.Ryu and B.Liu, JCP 95, 5889. */
/*                This scheme reuses intermediate products between */
/*                2DX and 2DY integrals, which can be achieved by having */
/*                the outer loops run over all possible x and y monomial */
/*                parts and the inner loops over all allowed E and F */
/*                shell combinations. The price to pay for such loop */
/*                ordering is the scattered addressing of locations */
/*                within the batch array, which has its rows and columns */
/*                ordered such that the E and F shells are increasing */
/*                and within each E and F shell the monomials are */
/*                ordered such that x>y>z in exponents. */

/*                An example follows: */
/*                ------------------- */

/*                Let E = 0,2 and F = 0,1. Then we have the left and */
/*                right hand of the batch array ordered as follows: */

/*                          left xyz        right xyz */

/*                             000             000 */
/*                             ---             --- */
/*                             100             100 */
/*                             010             010 */
/*                             001             001 */
/*                             --- */
/*                             200 -> -5 */
/*                             110 */
/*                             101 -> -3 */
/*                             020 */
/*                             011 */
/*                             002 -> 0 */

/*                The batch would thus have dimensions 10 x 4. For */
/*                the left side (and analogous for the right side) the */
/*                reduced multiplication scheme would have its outer */
/*                most loop run over x=0,2, followed by the next loop */
/*                y=0,2-x. The innermost loop would then run over the */
/*                allowed shells E=E(max),max(E(min),x+y). In this */
/*                case all x,y-pairs can be reused for all appropriate */
/*                shell combinations. */

/*                To find the address of a specific x,y,z,E combination */
/*                inside the batch array, we first note that the z-part */
/*                is dependent on the x,y-parts and is hence not needed. */
/*                Lets look at the E-part first. The E-part is evaluated */
/*                from its dimension formula (E+1)*(E+2)/2. Organizing */
/*                the inner E-loop to run from E(max) always, the */
/*                dimension for E(max) is passed as the argument NXYZP */
/*                and all lower E dimensions are calculated by the */
/*                formula relating dimensions between E and E+1: */

/*                           dim(E+1) = dim(E) + E + 2 */

/*                In this way multiplications in the E-part can be */
/*                entirely avoided. The x,y-part is defined as the */
/*                part which has to be subtracted from dim(E) to */
/*                reach the xyz monomial position inside the E shell. */
/*                It can be divided into an x-part and a y-part. The */
/*                x-part is given by the fomula: */

/*                          x-part = - x*E + x(x-3)/2 */

/*                and for the example above has been given for E=2 */
/*                and marked with arrows ->. The last term of the x-part */
/*                involves 1 multiplication and division, however it */
/*                can be changed to: */

/*                                               x-1 */
/*                          x-part = - x*E - x + sum i */
/*                                               i=0 */

/*                and clever additions inside the x-loop avoid the use */
/*                of multiplications and divisions. The y-part is trivial */
/*                and is simply equal to -y. The overall conclusion is */
/*                thus that the location of a specific x,y,z,E quadruple */
/*                inside the batch comes at the cost of one x*E(max) */
/*                multiplication in the outermost x-loops, since the */
/*                other x*E ones can again be reached via stepwise */
/*                subtraction of x from x*E(max). */

/*                Due to the very computational intensive steps inside */
/*                the x,y,z,E loops, special sections of identical */
/*                x,y,z,E loop structers have been given for each */
/*                # of roots =< 9, thus saving considerable computing */
/*                time over the general case. */

/*                For comments on how the x,y,z,E loop structures are */
/*                coded please refer to the general root case. */


/*                  Input: */

/*                    SHELLx      =  shell types for individual csh */
/*                                   x=A,C and csh sums P=A+B,Q=C+D */
/*                    NGQP        =  # of gaussian quadrature points */
/*                                   (roots) */
/*                    NEXQ        =  current # of exponent quadruplets */
/*                    NGQEXQ      =  product of # of gaussian quadrature */
/*                                   points times exponent quadruplets */
/*                    NXYZE(F)T   =  sum of # of cartesian monomials */
/*                                   for all shells in the range */
/*                                   E = A,...,P=A+B and in the range */
/*                                   F = C,...,Q=C+D */
/*                    NXYZy       =  # of cartesian monomials for */
/*                                   y = P,Q shells */
/*                    INT2Dx      =  all current 2D PQ integrals for */
/*                                   each cartesian component */
/*                                   (x = X,Y,Z) */
/*                    TEMP1(2)    =  scratch arrays holding intermediate */
/*                                   2D PQ integral products */
/*                    SCALE       =  the NGQEXQ scaling factors */


/*                  Output: */

/*                    BATCH       =  batch of primitive cartesian */
/*                                   [E0|F0] integrals corresponding */
/*                                   to all current exponent quadruplets */


/*  AUTHOR      : Norbert Flocke */
/* ------------------------------------------------------------------------ */


/*             ...include files and declare variables. */




/* ------------------------------------------------------------------------ */


/*             ...jump according to number of roots. */


    /* Parameter adjustments */
    --scale;
    --temp2;
    --temp1;
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
    batch_dim1 = *nexq - 1 + 1;
    batch_dim2 = *nxyzet - 1 + 1;
    batch_offset = 1 + batch_dim1 * (1 + batch_dim2 * 1);
    batch -= batch_offset;

    switch (MIN (*ngqp, 10))
    {
    case 1:
        goto L1;
    case 2:
        goto L2;
    case 3:
        goto L3;
    case 4:
        goto L4;
    case 5:
        goto L5;
    case 6:
        goto L6;
    case 7:
        goto L7;
    case 8:
        goto L8;
    case 9:
        goto L9;
    case 10:
        goto L10;
    }


/*                       ******************** */
/*                       *  # of roots = 1  * */
/*                       ******************** */


  L1:
    xfp = *nxyzft + 3;
    i__1 = *shellq;
    for (xf = 0; xf <= i__1; ++xf)
    {
        xfp = xfp + xf - 2;
        xfmax = xf * *shellq;
        yfend = *shellq - xf;
        xep = *nxyzet + 3;
        i__2 = *shellp;
        for (xe = 0; xe <= i__2; ++xe)
        {
            xep = xep + xe - 2;
            xemax = xe * *shellp;
            yeend = *shellp - xe;
            i__3 = *nexq;
            for (m = 1; m <= i__3; ++m)
            {
                temp1[m] = scale[m] * int2dx[m + (xe + xf * int2dx_dim2) *
                                             int2dx_dim1];
            }
            xyfp = xfp - xfmax;
            i__3 = yfend;
            for (yf = 0; yf <= i__3; ++yf)
            {
                xyf = xf + yf;
                --xyfp;
                sfend = max (*shellc, xyf);
                xyep = xep - xemax;
                i__4 = yeend;
                for (ye = 0; ye <= i__4; ++ye)
                {
                    xye = xe + ye;
                    --xyep;
                    seend = max (*shella, xye);
                    if (ye + yf == 0)
                    {
                        i__5 = *nexq;
                        for (m = 1; m <= i__5; ++m)
                        {
                            temp2[m] = temp1[m];
                        }
                    }
                    else
                    {
                        i__5 = *nexq;
                        for (m = 1; m <= i__5; ++m)
                        {
                            temp2[m] = temp1[m] * int2dy[m + (ye + yf *
                                                              int2dy_dim2) *
                                                         int2dy_dim1];
                        }
                    }
                    j = xyfp;
                    nxyzf = *nxyzq;
                    i__5 = sfend;
                    for (sf = *shellq; sf >= i__5; --sf)
                    {
                        zf = sf - xyf;
                        i__ = xyep;
                        nxyze = *nxyzp;
                        i__6 = seend;
                        for (se = *shellp; se >= i__6; --se)
                        {
                            ze = se - xye;
                            if (ze + zf == 0)
                            {
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[m];
                                }
                            }
                            else
                            {
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[m] * int2dz[m
                                                                          +
                                                                          (ze
                                                                           +
                                                                           zf
                                                                           *
                                                                           int2dz_dim2)
                                                                          *
                                                                          int2dz_dim1];
                                }
                            }
                            i__ = i__ - nxyze + xe;
                            nxyze = nxyze - se - 1;
/* L122: */
                        }
                        j = j - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
/* L120: */
                    }
/* L112: */
                }
/* L110: */
            }
/* L102: */
        }
/* L100: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 2  * */
/*                       ******************** */


  L2:
    xfp = *nxyzft + 3;
    i__1 = *shellq;
    for (xf = 0; xf <= i__1; ++xf)
    {
        xfp = xfp + xf - 2;
        xfmax = xf * *shellq;
        yfend = *shellq - xf;
        xep = *nxyzet + 3;
        i__2 = *shellp;
        for (xe = 0; xe <= i__2; ++xe)
        {
            xep = xep + xe - 2;
            xemax = xe * *shellp;
            yeend = *shellp - xe;
            i__3 = *ngqexq;
            for (m = 1; m <= i__3; ++m)
            {
                temp1[m] = scale[m] * int2dx[m + (xe + xf * int2dx_dim2) *
                                             int2dx_dim1];
            }
            xyfp = xfp - xfmax;
            i__3 = yfend;
            for (yf = 0; yf <= i__3; ++yf)
            {
                xyf = xf + yf;
                --xyfp;
                sfend = max (*shellc, xyf);
                xyep = xep - xemax;
                i__4 = yeend;
                for (ye = 0; ye <= i__4; ++ye)
                {
                    xye = xe + ye;
                    --xyep;
                    seend = max (*shella, xye);
                    if (ye + yf == 0)
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n];
                        }
                    }
                    else
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n] * int2dy[n + (ye + yf *
                                                              int2dy_dim2) *
                                                         int2dy_dim1];
                        }
                    }
                    j = xyfp;
                    nxyzf = *nxyzq;
                    i__5 = sfend;
                    for (sf = *shellq; sf >= i__5; --sf)
                    {
                        zf = sf - xyf;
                        i__ = xyep;
                        nxyze = *nxyzp;
                        i__6 = seend;
                        for (se = *shellp; se >= i__6; --se)
                        {
                            ze = se - xye;
                            if (ze + zf == 0)
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] + temp2[k
                                                                         + 1];
                                    k += 2;
                                }
                            }
                            else
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] * int2dz[k
                                                                          +
                                                                          (ze
                                                                           +
                                                                           zf
                                                                           *
                                                                           int2dz_dim2)
                                                                          *
                                                                          int2dz_dim1]
                                        + temp2[k + 1] * int2dz[k + 1 +
                                                                (ze +
                                                                 zf *
                                                                 int2dz_dim2)
                                                                *
                                                                int2dz_dim1];
                                    k += 2;
                                }
                            }
                            i__ = i__ - nxyze + xe;
                            nxyze = nxyze - se - 1;
/* L222: */
                        }
                        j = j - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
/* L220: */
                    }
/* L212: */
                }
/* L210: */
            }
/* L202: */
        }
/* L200: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 3  * */
/*                       ******************** */


  L3:
    xfp = *nxyzft + 3;
    i__1 = *shellq;
    for (xf = 0; xf <= i__1; ++xf)
    {
        xfp = xfp + xf - 2;
        xfmax = xf * *shellq;
        yfend = *shellq - xf;
        xep = *nxyzet + 3;
        i__2 = *shellp;
        for (xe = 0; xe <= i__2; ++xe)
        {
            xep = xep + xe - 2;
            xemax = xe * *shellp;
            yeend = *shellp - xe;
            i__3 = *ngqexq;
            for (m = 1; m <= i__3; ++m)
            {
                temp1[m] = scale[m] * int2dx[m + (xe + xf * int2dx_dim2) *
                                             int2dx_dim1];
            }
            xyfp = xfp - xfmax;
            i__3 = yfend;
            for (yf = 0; yf <= i__3; ++yf)
            {
                xyf = xf + yf;
                --xyfp;
                sfend = max (*shellc, xyf);
                xyep = xep - xemax;
                i__4 = yeend;
                for (ye = 0; ye <= i__4; ++ye)
                {
                    xye = xe + ye;
                    --xyep;
                    seend = max (*shella, xye);
                    if (ye + yf == 0)
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n];
                        }
                    }
                    else
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n] * int2dy[n + (ye + yf *
                                                              int2dy_dim2) *
                                                         int2dy_dim1];
                        }
                    }
                    j = xyfp;
                    nxyzf = *nxyzq;
                    i__5 = sfend;
                    for (sf = *shellq; sf >= i__5; --sf)
                    {
                        zf = sf - xyf;
                        i__ = xyep;
                        nxyze = *nxyzp;
                        i__6 = seend;
                        for (se = *shellp; se >= i__6; --se)
                        {
                            ze = se - xye;
                            if (ze + zf == 0)
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] + temp2[k
                                                                         +
                                                                         1] +
                                        temp2[k + 2];
                                    k += 3;
                                }
                            }
                            else
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] * int2dz[k
                                                                          +
                                                                          (ze
                                                                           +
                                                                           zf
                                                                           *
                                                                           int2dz_dim2)
                                                                          *
                                                                          int2dz_dim1]
                                        + temp2[k + 1] * int2dz[k + 1 +
                                                                (ze +
                                                                 zf *
                                                                 int2dz_dim2)
                                                                *
                                                                int2dz_dim1] +
                                        temp2[k + 2] * int2dz[k + 2 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1];
                                    k += 3;
                                }
                            }
                            i__ = i__ - nxyze + xe;
                            nxyze = nxyze - se - 1;
/* L322: */
                        }
                        j = j - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
/* L320: */
                    }
/* L312: */
                }
/* L310: */
            }
/* L302: */
        }
/* L300: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 4  * */
/*                       ******************** */


  L4:
    xfp = *nxyzft + 3;
    i__1 = *shellq;
    for (xf = 0; xf <= i__1; ++xf)
    {
        xfp = xfp + xf - 2;
        xfmax = xf * *shellq;
        yfend = *shellq - xf;
        xep = *nxyzet + 3;
        i__2 = *shellp;
        for (xe = 0; xe <= i__2; ++xe)
        {
            xep = xep + xe - 2;
            xemax = xe * *shellp;
            yeend = *shellp - xe;
            i__3 = *ngqexq;
            for (m = 1; m <= i__3; ++m)
            {
                temp1[m] = scale[m] * int2dx[m + (xe + xf * int2dx_dim2) *
                                             int2dx_dim1];
            }
            xyfp = xfp - xfmax;
            i__3 = yfend;
            for (yf = 0; yf <= i__3; ++yf)
            {
                xyf = xf + yf;
                --xyfp;
                sfend = max (*shellc, xyf);
                xyep = xep - xemax;
                i__4 = yeend;
                for (ye = 0; ye <= i__4; ++ye)
                {
                    xye = xe + ye;
                    --xyep;
                    seend = max (*shella, xye);
                    if (ye + yf == 0)
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n];
                        }
                    }
                    else
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n] * int2dy[n + (ye + yf *
                                                              int2dy_dim2) *
                                                         int2dy_dim1];
                        }
                    }
                    j = xyfp;
                    nxyzf = *nxyzq;
                    i__5 = sfend;
                    for (sf = *shellq; sf >= i__5; --sf)
                    {
                        zf = sf - xyf;
                        i__ = xyep;
                        nxyze = *nxyzp;
                        i__6 = seend;
                        for (se = *shellp; se >= i__6; --se)
                        {
                            ze = se - xye;
                            if (ze + zf == 0)
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] + temp2[k
                                                                         +
                                                                         1] +
                                        temp2[k + 2] + temp2[k + 3];
                                    k += 4;
                                }
                            }
                            else
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] * int2dz[k
                                                                          +
                                                                          (ze
                                                                           +
                                                                           zf
                                                                           *
                                                                           int2dz_dim2)
                                                                          *
                                                                          int2dz_dim1]
                                        + temp2[k + 1] * int2dz[k + 1 +
                                                                (ze +
                                                                 zf *
                                                                 int2dz_dim2)
                                                                *
                                                                int2dz_dim1] +
                                        temp2[k + 2] * int2dz[k + 2 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 3] * int2dz[k + 3 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1];
                                    k += 4;
                                }
                            }
                            i__ = i__ - nxyze + xe;
                            nxyze = nxyze - se - 1;
/* L422: */
                        }
                        j = j - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
/* L420: */
                    }
/* L412: */
                }
/* L410: */
            }
/* L402: */
        }
/* L400: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 5  * */
/*                       ******************** */


  L5:
    xfp = *nxyzft + 3;
    i__1 = *shellq;
    for (xf = 0; xf <= i__1; ++xf)
    {
        xfp = xfp + xf - 2;
        xfmax = xf * *shellq;
        yfend = *shellq - xf;
        xep = *nxyzet + 3;
        i__2 = *shellp;
        for (xe = 0; xe <= i__2; ++xe)
        {
            xep = xep + xe - 2;
            xemax = xe * *shellp;
            yeend = *shellp - xe;
            i__3 = *ngqexq;
            for (m = 1; m <= i__3; ++m)
            {
                temp1[m] = scale[m] * int2dx[m + (xe + xf * int2dx_dim2) *
                                             int2dx_dim1];
            }
            xyfp = xfp - xfmax;
            i__3 = yfend;
            for (yf = 0; yf <= i__3; ++yf)
            {
                xyf = xf + yf;
                --xyfp;
                sfend = max (*shellc, xyf);
                xyep = xep - xemax;
                i__4 = yeend;
                for (ye = 0; ye <= i__4; ++ye)
                {
                    xye = xe + ye;
                    --xyep;
                    seend = max (*shella, xye);
                    if (ye + yf == 0)
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n];
                        }
                    }
                    else
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n] * int2dy[n + (ye + yf *
                                                              int2dy_dim2) *
                                                         int2dy_dim1];
                        }
                    }
                    j = xyfp;
                    nxyzf = *nxyzq;
                    i__5 = sfend;
                    for (sf = *shellq; sf >= i__5; --sf)
                    {
                        zf = sf - xyf;
                        i__ = xyep;
                        nxyze = *nxyzp;
                        i__6 = seend;
                        for (se = *shellp; se >= i__6; --se)
                        {
                            ze = se - xye;
                            if (ze + zf == 0)
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] + temp2[k
                                                                         +
                                                                         1] +
                                        temp2[k + 2] + temp2[k + 3] +
                                        temp2[k + 4];
                                    k += 5;
                                }
                            }
                            else
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] * int2dz[k
                                                                          +
                                                                          (ze
                                                                           +
                                                                           zf
                                                                           *
                                                                           int2dz_dim2)
                                                                          *
                                                                          int2dz_dim1]
                                        + temp2[k + 1] * int2dz[k + 1 +
                                                                (ze +
                                                                 zf *
                                                                 int2dz_dim2)
                                                                *
                                                                int2dz_dim1] +
                                        temp2[k + 2] * int2dz[k + 2 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 3] * int2dz[k + 3 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 4] * int2dz[k + 4 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1];
                                    k += 5;
                                }
                            }
                            i__ = i__ - nxyze + xe;
                            nxyze = nxyze - se - 1;
/* L522: */
                        }
                        j = j - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
/* L520: */
                    }
/* L512: */
                }
/* L510: */
            }
/* L502: */
        }
/* L500: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 6  * */
/*                       ******************** */


  L6:
    xfp = *nxyzft + 3;
    i__1 = *shellq;
    for (xf = 0; xf <= i__1; ++xf)
    {
        xfp = xfp + xf - 2;
        xfmax = xf * *shellq;
        yfend = *shellq - xf;
        xep = *nxyzet + 3;
        i__2 = *shellp;
        for (xe = 0; xe <= i__2; ++xe)
        {
            xep = xep + xe - 2;
            xemax = xe * *shellp;
            yeend = *shellp - xe;
            i__3 = *ngqexq;
            for (m = 1; m <= i__3; ++m)
            {
                temp1[m] = scale[m] * int2dx[m + (xe + xf * int2dx_dim2) *
                                             int2dx_dim1];
            }
            xyfp = xfp - xfmax;
            i__3 = yfend;
            for (yf = 0; yf <= i__3; ++yf)
            {
                xyf = xf + yf;
                --xyfp;
                sfend = max (*shellc, xyf);
                xyep = xep - xemax;
                i__4 = yeend;
                for (ye = 0; ye <= i__4; ++ye)
                {
                    xye = xe + ye;
                    --xyep;
                    seend = max (*shella, xye);
                    if (ye + yf == 0)
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n];
                        }
                    }
                    else
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n] * int2dy[n + (ye + yf *
                                                              int2dy_dim2) *
                                                         int2dy_dim1];
                        }
                    }
                    j = xyfp;
                    nxyzf = *nxyzq;
                    i__5 = sfend;
                    for (sf = *shellq; sf >= i__5; --sf)
                    {
                        zf = sf - xyf;
                        i__ = xyep;
                        nxyze = *nxyzp;
                        i__6 = seend;
                        for (se = *shellp; se >= i__6; --se)
                        {
                            ze = se - xye;
                            if (ze + zf == 0)
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] + temp2[k
                                                                         +
                                                                         1] +
                                        temp2[k + 2] + temp2[k + 3] +
                                        temp2[k + 4] + temp2[k + 5];
                                    k += 6;
                                }
                            }
                            else
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] * int2dz[k
                                                                          +
                                                                          (ze
                                                                           +
                                                                           zf
                                                                           *
                                                                           int2dz_dim2)
                                                                          *
                                                                          int2dz_dim1]
                                        + temp2[k + 1] * int2dz[k + 1 +
                                                                (ze +
                                                                 zf *
                                                                 int2dz_dim2)
                                                                *
                                                                int2dz_dim1] +
                                        temp2[k + 2] * int2dz[k + 2 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 3] * int2dz[k + 3 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 4] * int2dz[k + 4 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 5] * int2dz[k + 5 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1];
                                    k += 6;
                                }
                            }
                            i__ = i__ - nxyze + xe;
                            nxyze = nxyze - se - 1;
/* L622: */
                        }
                        j = j - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
/* L620: */
                    }
/* L612: */
                }
/* L610: */
            }
/* L602: */
        }
/* L600: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 7  * */
/*                       ******************** */


  L7:
    xfp = *nxyzft + 3;
    i__1 = *shellq;
    for (xf = 0; xf <= i__1; ++xf)
    {
        xfp = xfp + xf - 2;
        xfmax = xf * *shellq;
        yfend = *shellq - xf;
        xep = *nxyzet + 3;
        i__2 = *shellp;
        for (xe = 0; xe <= i__2; ++xe)
        {
            xep = xep + xe - 2;
            xemax = xe * *shellp;
            yeend = *shellp - xe;
            i__3 = *ngqexq;
            for (m = 1; m <= i__3; ++m)
            {
                temp1[m] = scale[m] * int2dx[m + (xe + xf * int2dx_dim2) *
                                             int2dx_dim1];
            }
            xyfp = xfp - xfmax;
            i__3 = yfend;
            for (yf = 0; yf <= i__3; ++yf)
            {
                xyf = xf + yf;
                --xyfp;
                sfend = max (*shellc, xyf);
                xyep = xep - xemax;
                i__4 = yeend;
                for (ye = 0; ye <= i__4; ++ye)
                {
                    xye = xe + ye;
                    --xyep;
                    seend = max (*shella, xye);
                    if (ye + yf == 0)
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n];
                        }
                    }
                    else
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n] * int2dy[n + (ye + yf *
                                                              int2dy_dim2) *
                                                         int2dy_dim1];
                        }
                    }
                    j = xyfp;
                    nxyzf = *nxyzq;
                    i__5 = sfend;
                    for (sf = *shellq; sf >= i__5; --sf)
                    {
                        zf = sf - xyf;
                        i__ = xyep;
                        nxyze = *nxyzp;
                        i__6 = seend;
                        for (se = *shellp; se >= i__6; --se)
                        {
                            ze = se - xye;
                            if (ze + zf == 0)
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] + temp2[k
                                                                         +
                                                                         1] +
                                        temp2[k + 2] + temp2[k + 3] +
                                        temp2[k + 4] + temp2[k + 5] +
                                        temp2[k + 6];
                                    k += 7;
                                }
                            }
                            else
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] * int2dz[k
                                                                          +
                                                                          (ze
                                                                           +
                                                                           zf
                                                                           *
                                                                           int2dz_dim2)
                                                                          *
                                                                          int2dz_dim1]
                                        + temp2[k + 1] * int2dz[k + 1 +
                                                                (ze +
                                                                 zf *
                                                                 int2dz_dim2)
                                                                *
                                                                int2dz_dim1] +
                                        temp2[k + 2] * int2dz[k + 2 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 3] * int2dz[k + 3 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 4] * int2dz[k + 4 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 5] * int2dz[k + 5 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 6] * int2dz[k + 6 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1];
                                    k += 7;
                                }
                            }
                            i__ = i__ - nxyze + xe;
                            nxyze = nxyze - se - 1;
/* L722: */
                        }
                        j = j - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
/* L720: */
                    }
/* L712: */
                }
/* L710: */
            }
/* L702: */
        }
/* L700: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 8  * */
/*                       ******************** */


  L8:
    xfp = *nxyzft + 3;
    i__1 = *shellq;
    for (xf = 0; xf <= i__1; ++xf)
    {
        xfp = xfp + xf - 2;
        xfmax = xf * *shellq;
        yfend = *shellq - xf;
        xep = *nxyzet + 3;
        i__2 = *shellp;
        for (xe = 0; xe <= i__2; ++xe)
        {
            xep = xep + xe - 2;
            xemax = xe * *shellp;
            yeend = *shellp - xe;
            i__3 = *ngqexq;
            for (m = 1; m <= i__3; ++m)
            {
                temp1[m] = scale[m] * int2dx[m + (xe + xf * int2dx_dim2) *
                                             int2dx_dim1];
            }
            xyfp = xfp - xfmax;
            i__3 = yfend;
            for (yf = 0; yf <= i__3; ++yf)
            {
                xyf = xf + yf;
                --xyfp;
                sfend = max (*shellc, xyf);
                xyep = xep - xemax;
                i__4 = yeend;
                for (ye = 0; ye <= i__4; ++ye)
                {
                    xye = xe + ye;
                    --xyep;
                    seend = max (*shella, xye);
                    if (ye + yf == 0)
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n];
                        }
                    }
                    else
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n] * int2dy[n + (ye + yf *
                                                              int2dy_dim2) *
                                                         int2dy_dim1];
                        }
                    }
                    j = xyfp;
                    nxyzf = *nxyzq;
                    i__5 = sfend;
                    for (sf = *shellq; sf >= i__5; --sf)
                    {
                        zf = sf - xyf;
                        i__ = xyep;
                        nxyze = *nxyzp;
                        i__6 = seend;
                        for (se = *shellp; se >= i__6; --se)
                        {
                            ze = se - xye;
                            if (ze + zf == 0)
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] + temp2[k
                                                                         +
                                                                         1] +
                                        temp2[k + 2] + temp2[k + 3] +
                                        temp2[k + 4] + temp2[k + 5] +
                                        temp2[k + 6] + temp2[k + 7];
                                    k += 8;
                                }
                            }
                            else
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] * int2dz[k
                                                                          +
                                                                          (ze
                                                                           +
                                                                           zf
                                                                           *
                                                                           int2dz_dim2)
                                                                          *
                                                                          int2dz_dim1]
                                        + temp2[k + 1] * int2dz[k + 1 +
                                                                (ze +
                                                                 zf *
                                                                 int2dz_dim2)
                                                                *
                                                                int2dz_dim1] +
                                        temp2[k + 2] * int2dz[k + 2 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 3] * int2dz[k + 3 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 4] * int2dz[k + 4 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 5] * int2dz[k + 5 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 6] * int2dz[k + 6 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 7] * int2dz[k + 7 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1];
                                    k += 8;
                                }
                            }
                            i__ = i__ - nxyze + xe;
                            nxyze = nxyze - se - 1;
/* L822: */
                        }
                        j = j - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
/* L820: */
                    }
/* L812: */
                }
/* L810: */
            }
/* L802: */
        }
/* L800: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 9  * */
/*                       ******************** */


  L9:
    xfp = *nxyzft + 3;
    i__1 = *shellq;
    for (xf = 0; xf <= i__1; ++xf)
    {
        xfp = xfp + xf - 2;
        xfmax = xf * *shellq;
        yfend = *shellq - xf;
        xep = *nxyzet + 3;
        i__2 = *shellp;
        for (xe = 0; xe <= i__2; ++xe)
        {
            xep = xep + xe - 2;
            xemax = xe * *shellp;
            yeend = *shellp - xe;
            i__3 = *ngqexq;
            for (m = 1; m <= i__3; ++m)
            {
                temp1[m] = scale[m] * int2dx[m + (xe + xf * int2dx_dim2) *
                                             int2dx_dim1];
            }
            xyfp = xfp - xfmax;
            i__3 = yfend;
            for (yf = 0; yf <= i__3; ++yf)
            {
                xyf = xf + yf;
                --xyfp;
                sfend = max (*shellc, xyf);
                xyep = xep - xemax;
                i__4 = yeend;
                for (ye = 0; ye <= i__4; ++ye)
                {
                    xye = xe + ye;
                    --xyep;
                    seend = max (*shella, xye);
                    if (ye + yf == 0)
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n];
                        }
                    }
                    else
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n] * int2dy[n + (ye + yf *
                                                              int2dy_dim2) *
                                                         int2dy_dim1];
                        }
                    }
                    j = xyfp;
                    nxyzf = *nxyzq;
                    i__5 = sfend;
                    for (sf = *shellq; sf >= i__5; --sf)
                    {
                        zf = sf - xyf;
                        i__ = xyep;
                        nxyze = *nxyzp;
                        i__6 = seend;
                        for (se = *shellp; se >= i__6; --se)
                        {
                            ze = se - xye;
                            if (ze + zf == 0)
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] + temp2[k
                                                                         +
                                                                         1] +
                                        temp2[k + 2] + temp2[k + 3] +
                                        temp2[k + 4] + temp2[k + 5] +
                                        temp2[k + 6] + temp2[k + 7] +
                                        temp2[k + 8];
                                    k += 9;
                                }
                            }
                            else
                            {
                                k = 1;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = temp2[k] * int2dz[k
                                                                          +
                                                                          (ze
                                                                           +
                                                                           zf
                                                                           *
                                                                           int2dz_dim2)
                                                                          *
                                                                          int2dz_dim1]
                                        + temp2[k + 1] * int2dz[k + 1 +
                                                                (ze +
                                                                 zf *
                                                                 int2dz_dim2)
                                                                *
                                                                int2dz_dim1] +
                                        temp2[k + 2] * int2dz[k + 2 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 3] * int2dz[k + 3 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 4] * int2dz[k + 4 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 5] * int2dz[k + 5 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 6] * int2dz[k + 6 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 7] * int2dz[k + 7 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1] +
                                        temp2[k + 8] * int2dz[k + 8 +
                                                              (ze +
                                                               zf *
                                                               int2dz_dim2) *
                                                              int2dz_dim1];
                                    k += 9;
                                }
                            }
                            i__ = i__ - nxyze + xe;
                            nxyze = nxyze - se - 1;
/* L922: */
                        }
                        j = j - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
/* L920: */
                    }
/* L912: */
                }
/* L910: */
            }
/* L902: */
        }
/* L900: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots > 9  * */
/*                       ******************** */

/*             ...outer loops over x,x-pairs. No skipping of the */
/*                x,x-contribution of 0,0-type can be done here, */
/*                since the 2DX integrals carry the Rys weight! */


  L10:
    xfp = *nxyzft + 3;
    i__1 = *shellq;
    for (xf = 0; xf <= i__1; ++xf)
    {
        xfp = xfp + xf - 2;
        xfmax = xf * *shellq;
        yfend = *shellq - xf;
        xep = *nxyzet + 3;
        i__2 = *shellp;
        for (xe = 0; xe <= i__2; ++xe)
        {
            xep = xep + xe - 2;
            xemax = xe * *shellp;
            yeend = *shellp - xe;
            i__3 = *ngqexq;
            for (m = 1; m <= i__3; ++m)
            {
                temp1[m] = scale[m] * int2dx[m + (xe + xf * int2dx_dim2) *
                                             int2dx_dim1];
            }


/*             ...middle loops over y,y-pairs. Skip multiplication */
/*                of y,y-contributions, if we have a 0,0-pair, as */
/*                then the 2DY integral is equal to 1. */


            xyfp = xfp - xfmax;
            i__3 = yfend;
            for (yf = 0; yf <= i__3; ++yf)
            {
                xyf = xf + yf;
                --xyfp;
                sfend = max (*shellc, xyf);
                xyep = xep - xemax;
                i__4 = yeend;
                for (ye = 0; ye <= i__4; ++ye)
                {
                    xye = xe + ye;
                    --xyep;
                    seend = max (*shella, xye);
                    if (ye + yf == 0)
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n];
                        }
                    }
                    else
                    {
                        i__5 = *ngqexq;
                        for (n = 1; n <= i__5; ++n)
                        {
                            temp2[n] = temp1[n] * int2dy[n + (ye + yf *
                                                              int2dy_dim2) *
                                                         int2dy_dim1];
                        }
                    }


/*             ...inner loops over E,F-pairs. Skip multiplication */
/*                of z,z-contributions, if we have a 0,0-pair, as */
/*                then the 2DZ integral is equal to 1. */


                    j = xyfp;
                    nxyzf = *nxyzq;
                    i__5 = sfend;
                    for (sf = *shellq; sf >= i__5; --sf)
                    {
                        zf = sf - xyf;
                        i__ = xyep;
                        nxyze = *nxyzp;
                        i__6 = seend;
                        for (se = *shellp; se >= i__6; --se)
                        {
                            ze = se - xye;


/*             ...all info concerning all three x,x-, y,y- and z,z-pairs */
/*                have been collected for all exponent quadruplets at */
/*                once. Sum up the 2D X,Y,Z integral products to the */
/*                appropriate place of the [E0|F0] batch. */


                            if (ze + zf == 0)
                            {
                                k = 0;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    sum = 0.;
                                    i__8 = *ngqp;
                                    for (n = 1; n <= i__8; ++n)
                                    {
                                        sum += temp2[k + n];
                                    }
                                    k += *ngqp;
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = sum;
                                }
                            }
                            else
                            {
                                k = 0;
                                i__7 = *nexq;
                                for (m = 1; m <= i__7; ++m)
                                {
                                    sum = 0.;
                                    i__8 = *ngqp;
                                    for (n = 1; n <= i__8; ++n)
                                    {
                                        sum +=
                                            temp2[k + n] * int2dz[k + n +
                                                                  (ze +
                                                                   zf *
                                                                   int2dz_dim2)
                                                                  *
                                                                  int2dz_dim1];
                                    }
                                    k += *ngqp;
                                    batch[m + (i__ + j * batch_dim2) *
                                          batch_dim1] = sum;
                                }
                            }


/*             ...next z,z-pair. */


                            i__ = i__ - nxyze + xe;
                            nxyze = nxyze - se - 1;
/* L1022: */
                        }
                        j = j - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
/* L1020: */
                    }


/*             ...next y,y-pair and next x,x-pair. */


/* L1012: */
                }
/* L1010: */
            }
/* L1002: */
        }
/* L1000: */
    }

    return 0;
}
