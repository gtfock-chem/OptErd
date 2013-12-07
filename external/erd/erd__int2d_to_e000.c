/* erd__int2d_to_e000.f -- translated by f2c (version 20100827).
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
erd__int2d_to_e000_ (int * shella, int * shellp,
                      int * ngqp, int * nexq, int * ngqexq,
                      int * nxyzet, int * nxyzp, double * int2dx,
                      double * int2dy, double * int2dz,
                      double * temp1, double * temp2,
                      double * scale, double * batch)
{
    /* System generated locals */
    int batch_dim1, batch_offset, int2dx_dim1, int2dx_offset, int2dy_dim1,
        int2dy_offset, int2dz_dim1, int2dz_offset, i__1, i__2, i__3,
        i__4, i__5;

    /* Local variables */
    static int i__, k, m, n, se, xe, ye, ze, xep;
    static double sum;
    static int xye, xyep, seend, yeend, xemax, nxyze;

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__INT2D_TO_E000 */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine assembles the set of batches of cartesian */
/*                eris: */

/*                        [E0|00] or [00|E0] , E = A to P */

/*                adding up all the contributions from all the */
/*                respective 2D integrals: */

/*                                   P0 or 0P */

/*                Simplified version of the general E0F0 routine to */
/*                reduce loop overheads for those cases where there is */
/*                at least one s-shell on the bra or ket side. For */
/*                comments and details see the general E0F0 routine. */


/*                  Input: */

/*                    SHELLx      =  shell types for individual csh */
/*                                   x=A and csh sum P=A+B */
/*                    NGQP        =  # of gaussian quadrature points */
/*                                   (roots) */
/*                    NEXQ        =  current # of exponent quadruplets */
/*                    NGQEXQ      =  product of # of gaussian quadrature */
/*                                   points times exponent quadruplets */
/*                    NXYZET      =  sum of # of cartesian monomials */
/*                                   for all shells in the range */
/*                                   E = A,...,P=A+B */
/*                    NXYZP       =  # of cartesian monomials for the */
/*                                   P=A+B shell */
/*                    INT2Dx      =  all current 2D P0/0P integrals for */
/*                                   each cartesian component */
/*                                   (x = X,Y,Z) */
/*                    TEMP1(2)    =  scratch arrays holding intermediate */
/*                                   2D P0/0P integral products */
/*                    SCALE       =  the NGQEXQ scaling factors */


/*                  Output: */

/*                    BATCH       =  batch of primitive cartesian */
/*                                   [E0|00] integrals corresponding */
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
    int2dz_offset = 1 + int2dz_dim1 * 0;
    int2dz -= int2dz_offset;
    int2dy_dim1 = *ngqexq - 1 + 1;
    int2dy_offset = 1 + int2dy_dim1 * 0;
    int2dy -= int2dy_offset;
    int2dx_dim1 = *ngqexq - 1 + 1;
    int2dx_offset = 1 + int2dx_dim1 * 0;
    int2dx -= int2dx_offset;
    batch_dim1 = *nexq - 1 + 1;
    batch_offset = 1 + batch_dim1 * 1;
    batch -= batch_offset;

    /* Function Body */
    switch (min (*ngqp, 10))
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
    xep = *nxyzet + 3;
    i__1 = *shellp;
    for (xe = 0; xe <= i__1; ++xe)
    {
        xep = xep + xe - 2;
        xemax = xe * *shellp;
        yeend = *shellp - xe;
        i__2 = *nexq;
        for (m = 1; m <= i__2; ++m)
        {
            temp1[m] = scale[m] * int2dx[m + xe * int2dx_dim1];
        }
        xyep = xep - xemax;
        i__2 = yeend;
        for (ye = 0; ye <= i__2; ++ye)
        {
            xye = xe + ye;
            --xyep;
            seend = max (*shella, xye);
            if (ye == 0)
            {
                i__3 = *nexq;
                for (m = 1; m <= i__3; ++m)
                {
                    temp2[m] = temp1[m];
                }
            }
            else
            {
                i__3 = *nexq;
                for (m = 1; m <= i__3; ++m)
                {
                    temp2[m] = temp1[m] * int2dy[m + ye * int2dy_dim1];
                }
            }
            i__ = xyep;
            nxyze = *nxyzp;
            i__3 = seend;
            for (se = *shellp; se >= i__3; --se)
            {
                ze = se - xye;
                if (ze == 0)
                {
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[m];
                    }
                }
                else
                {
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[m] * int2dz[m +
                                                                        ze *
                                                                        int2dz_dim1];
                    }
                }
                i__ = i__ - nxyze + xe;
                nxyze = nxyze - se - 1;
/* L120: */
            }
/* L110: */
        }
/* L100: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 2  * */
/*                       ******************** */


  L2:
    xep = *nxyzet + 3;
    i__1 = *shellp;
    for (xe = 0; xe <= i__1; ++xe)
    {
        xep = xep + xe - 2;
        xemax = xe * *shellp;
        yeend = *shellp - xe;
        i__2 = *ngqexq;
        for (m = 1; m <= i__2; ++m)
        {
            temp1[m] = scale[m] * int2dx[m + xe * int2dx_dim1];
        }
        xyep = xep - xemax;
        i__2 = yeend;
        for (ye = 0; ye <= i__2; ++ye)
        {
            xye = xe + ye;
            --xyep;
            seend = max (*shella, xye);
            if (ye == 0)
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n];
                }
            }
            else
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n] * int2dy[n + ye * int2dy_dim1];
                }
            }
            i__ = xyep;
            nxyze = *nxyzp;
            i__3 = seend;
            for (se = *shellp; se >= i__3; --se)
            {
                ze = se - xye;
                if (ze == 0)
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] + temp2[k + 1];
                        k += 2;
                    }
                }
                else
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] * int2dz[k +
                                                                        ze *
                                                                        int2dz_dim1]
                            + temp2[k + 1] * int2dz[k + 1 + ze * int2dz_dim1];
                        k += 2;
                    }
                }
                i__ = i__ - nxyze + xe;
                nxyze = nxyze - se - 1;
/* L220: */
            }
/* L210: */
        }
/* L200: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 3  * */
/*                       ******************** */


  L3:
    xep = *nxyzet + 3;
    i__1 = *shellp;
    for (xe = 0; xe <= i__1; ++xe)
    {
        xep = xep + xe - 2;
        xemax = xe * *shellp;
        yeend = *shellp - xe;
        i__2 = *ngqexq;
        for (m = 1; m <= i__2; ++m)
        {
            temp1[m] = scale[m] * int2dx[m + xe * int2dx_dim1];
        }
        xyep = xep - xemax;
        i__2 = yeend;
        for (ye = 0; ye <= i__2; ++ye)
        {
            xye = xe + ye;
            --xyep;
            seend = max (*shella, xye);
            if (ye == 0)
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n];
                }
            }
            else
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n] * int2dy[n + ye * int2dy_dim1];
                }
            }
            i__ = xyep;
            nxyze = *nxyzp;
            i__3 = seend;
            for (se = *shellp; se >= i__3; --se)
            {
                ze = se - xye;
                if (ze == 0)
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] + temp2[k + 1]
                            + temp2[k + 2];
                        k += 3;
                    }
                }
                else
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] * int2dz[k +
                                                                        ze *
                                                                        int2dz_dim1]
                            + temp2[k + 1] * int2dz[k + 1 +
                                                    ze * int2dz_dim1] +
                            temp2[k + 2] * int2dz[k + 2 + ze * int2dz_dim1];
                        k += 3;
                    }
                }
                i__ = i__ - nxyze + xe;
                nxyze = nxyze - se - 1;
/* L320: */
            }
/* L310: */
        }
/* L300: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 4  * */
/*                       ******************** */


  L4:
    xep = *nxyzet + 3;
    i__1 = *shellp;
    for (xe = 0; xe <= i__1; ++xe)
    {
        xep = xep + xe - 2;
        xemax = xe * *shellp;
        yeend = *shellp - xe;
        i__2 = *ngqexq;
        for (m = 1; m <= i__2; ++m)
        {
            temp1[m] = scale[m] * int2dx[m + xe * int2dx_dim1];
        }
        xyep = xep - xemax;
        i__2 = yeend;
        for (ye = 0; ye <= i__2; ++ye)
        {
            xye = xe + ye;
            --xyep;
            seend = max (*shella, xye);
            if (ye == 0)
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n];
                }
            }
            else
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n] * int2dy[n + ye * int2dy_dim1];
                }
            }
            i__ = xyep;
            nxyze = *nxyzp;
            i__3 = seend;
            for (se = *shellp; se >= i__3; --se)
            {
                ze = se - xye;
                if (ze == 0)
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] + temp2[k + 1]
                            + temp2[k + 2] + temp2[k + 3];
                        k += 4;
                    }
                }
                else
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] * int2dz[k +
                                                                        ze *
                                                                        int2dz_dim1]
                            + temp2[k + 1] * int2dz[k + 1 +
                                                    ze * int2dz_dim1] +
                            temp2[k + 2] * int2dz[k + 2 + ze * int2dz_dim1] +
                            temp2[k + 3] * int2dz[k + 3 + ze * int2dz_dim1];
                        k += 4;
                    }
                }
                i__ = i__ - nxyze + xe;
                nxyze = nxyze - se - 1;
/* L420: */
            }
/* L410: */
        }
/* L400: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 5  * */
/*                       ******************** */


  L5:
    xep = *nxyzet + 3;
    i__1 = *shellp;
    for (xe = 0; xe <= i__1; ++xe)
    {
        xep = xep + xe - 2;
        xemax = xe * *shellp;
        yeend = *shellp - xe;
        i__2 = *ngqexq;
        for (m = 1; m <= i__2; ++m)
        {
            temp1[m] = scale[m] * int2dx[m + xe * int2dx_dim1];
        }
        xyep = xep - xemax;
        i__2 = yeend;
        for (ye = 0; ye <= i__2; ++ye)
        {
            xye = xe + ye;
            --xyep;
            seend = max (*shella, xye);
            if (ye == 0)
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n];
                }
            }
            else
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n] * int2dy[n + ye * int2dy_dim1];
                }
            }
            i__ = xyep;
            nxyze = *nxyzp;
            i__3 = seend;
            for (se = *shellp; se >= i__3; --se)
            {
                ze = se - xye;
                if (ze == 0)
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] + temp2[k + 1]
                            + temp2[k + 2] + temp2[k + 3] + temp2[k + 4];
                        k += 5;
                    }
                }
                else
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] * int2dz[k +
                                                                        ze *
                                                                        int2dz_dim1]
                            + temp2[k + 1] * int2dz[k + 1 +
                                                    ze * int2dz_dim1] +
                            temp2[k + 2] * int2dz[k + 2 + ze * int2dz_dim1] +
                            temp2[k + 3] * int2dz[k + 3 + ze * int2dz_dim1] +
                            temp2[k + 4] * int2dz[k + 4 + ze * int2dz_dim1];
                        k += 5;
                    }
                }
                i__ = i__ - nxyze + xe;
                nxyze = nxyze - se - 1;
/* L520: */
            }
/* L510: */
        }
/* L500: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 6  * */
/*                       ******************** */


  L6:
    xep = *nxyzet + 3;
    i__1 = *shellp;
    for (xe = 0; xe <= i__1; ++xe)
    {
        xep = xep + xe - 2;
        xemax = xe * *shellp;
        yeend = *shellp - xe;
        i__2 = *ngqexq;
        for (m = 1; m <= i__2; ++m)
        {
            temp1[m] = scale[m] * int2dx[m + xe * int2dx_dim1];
        }
        xyep = xep - xemax;
        i__2 = yeend;
        for (ye = 0; ye <= i__2; ++ye)
        {
            xye = xe + ye;
            --xyep;
            seend = max (*shella, xye);
            if (ye == 0)
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n];
                }
            }
            else
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n] * int2dy[n + ye * int2dy_dim1];
                }
            }
            i__ = xyep;
            nxyze = *nxyzp;
            i__3 = seend;
            for (se = *shellp; se >= i__3; --se)
            {
                ze = se - xye;
                if (ze == 0)
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] + temp2[k + 1]
                            + temp2[k + 2] + temp2[k + 3] + temp2[k + 4]
                            + temp2[k + 5];
                        k += 6;
                    }
                }
                else
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] * int2dz[k +
                                                                        ze *
                                                                        int2dz_dim1]
                            + temp2[k + 1] * int2dz[k + 1 +
                                                    ze * int2dz_dim1] +
                            temp2[k + 2] * int2dz[k + 2 + ze * int2dz_dim1] +
                            temp2[k + 3] * int2dz[k + 3 + ze * int2dz_dim1] +
                            temp2[k + 4] * int2dz[k + 4 + ze * int2dz_dim1] +
                            temp2[k + 5] * int2dz[k + 5 + ze * int2dz_dim1];
                        k += 6;
                    }
                }
                i__ = i__ - nxyze + xe;
                nxyze = nxyze - se - 1;
/* L620: */
            }
/* L610: */
        }
/* L600: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 7  * */
/*                       ******************** */


  L7:
    xep = *nxyzet + 3;
    i__1 = *shellp;
    for (xe = 0; xe <= i__1; ++xe)
    {
        xep = xep + xe - 2;
        xemax = xe * *shellp;
        yeend = *shellp - xe;
        i__2 = *ngqexq;
        for (m = 1; m <= i__2; ++m)
        {
            temp1[m] = scale[m] * int2dx[m + xe * int2dx_dim1];
        }
        xyep = xep - xemax;
        i__2 = yeend;
        for (ye = 0; ye <= i__2; ++ye)
        {
            xye = xe + ye;
            --xyep;
            seend = max (*shella, xye);
            if (ye == 0)
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n];
                }
            }
            else
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n] * int2dy[n + ye * int2dy_dim1];
                }
            }
            i__ = xyep;
            nxyze = *nxyzp;
            i__3 = seend;
            for (se = *shellp; se >= i__3; --se)
            {
                ze = se - xye;
                if (ze == 0)
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] + temp2[k + 1]
                            + temp2[k + 2] + temp2[k + 3] + temp2[k + 4]
                            + temp2[k + 5] + temp2[k + 6];
                        k += 7;
                    }
                }
                else
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] * int2dz[k +
                                                                        ze *
                                                                        int2dz_dim1]
                            + temp2[k + 1] * int2dz[k + 1 +
                                                    ze * int2dz_dim1] +
                            temp2[k + 2] * int2dz[k + 2 + ze * int2dz_dim1] +
                            temp2[k + 3] * int2dz[k + 3 + ze * int2dz_dim1] +
                            temp2[k + 4] * int2dz[k + 4 + ze * int2dz_dim1] +
                            temp2[k + 5] * int2dz[k + 5 + ze * int2dz_dim1] +
                            temp2[k + 6] * int2dz[k + 6 + ze * int2dz_dim1];
                        k += 7;
                    }
                }
                i__ = i__ - nxyze + xe;
                nxyze = nxyze - se - 1;
/* L720: */
            }
/* L710: */
        }
/* L700: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 8  * */
/*                       ******************** */


  L8:
    xep = *nxyzet + 3;
    i__1 = *shellp;
    for (xe = 0; xe <= i__1; ++xe)
    {
        xep = xep + xe - 2;
        xemax = xe * *shellp;
        yeend = *shellp - xe;
        i__2 = *ngqexq;
        for (m = 1; m <= i__2; ++m)
        {
            temp1[m] = scale[m] * int2dx[m + xe * int2dx_dim1];
        }
        xyep = xep - xemax;
        i__2 = yeend;
        for (ye = 0; ye <= i__2; ++ye)
        {
            xye = xe + ye;
            --xyep;
            seend = max (*shella, xye);
            if (ye == 0)
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n];
                }
            }
            else
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n] * int2dy[n + ye * int2dy_dim1];
                }
            }
            i__ = xyep;
            nxyze = *nxyzp;
            i__3 = seend;
            for (se = *shellp; se >= i__3; --se)
            {
                ze = se - xye;
                if (ze == 0)
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] + temp2[k + 1]
                            + temp2[k + 2] + temp2[k + 3] + temp2[k + 4]
                            + temp2[k + 5] + temp2[k + 6] + temp2[k + 7];
                        k += 8;
                    }
                }
                else
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] * int2dz[k +
                                                                        ze *
                                                                        int2dz_dim1]
                            + temp2[k + 1] * int2dz[k + 1 +
                                                    ze * int2dz_dim1] +
                            temp2[k + 2] * int2dz[k + 2 + ze * int2dz_dim1] +
                            temp2[k + 3] * int2dz[k + 3 + ze * int2dz_dim1] +
                            temp2[k + 4] * int2dz[k + 4 + ze * int2dz_dim1] +
                            temp2[k + 5] * int2dz[k + 5 + ze * int2dz_dim1] +
                            temp2[k + 6] * int2dz[k + 6 + ze * int2dz_dim1] +
                            temp2[k + 7] * int2dz[k + 7 + ze * int2dz_dim1];
                        k += 8;
                    }
                }
                i__ = i__ - nxyze + xe;
                nxyze = nxyze - se - 1;
/* L820: */
            }
/* L810: */
        }
/* L800: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots = 9  * */
/*                       ******************** */


  L9:
    xep = *nxyzet + 3;
    i__1 = *shellp;
    for (xe = 0; xe <= i__1; ++xe)
    {
        xep = xep + xe - 2;
        xemax = xe * *shellp;
        yeend = *shellp - xe;
        i__2 = *ngqexq;
        for (m = 1; m <= i__2; ++m)
        {
            temp1[m] = scale[m] * int2dx[m + xe * int2dx_dim1];
        }
        xyep = xep - xemax;
        i__2 = yeend;
        for (ye = 0; ye <= i__2; ++ye)
        {
            xye = xe + ye;
            --xyep;
            seend = max (*shella, xye);
            if (ye == 0)
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n];
                }
            }
            else
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n] * int2dy[n + ye * int2dy_dim1];
                }
            }
            i__ = xyep;
            nxyze = *nxyzp;
            i__3 = seend;
            for (se = *shellp; se >= i__3; --se)
            {
                ze = se - xye;
                if (ze == 0)
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] + temp2[k + 1]
                            + temp2[k + 2] + temp2[k + 3] + temp2[k + 4]
                            + temp2[k + 5] + temp2[k + 6] + temp2[k + 7]
                            + temp2[k + 8];
                        k += 9;
                    }
                }
                else
                {
                    k = 1;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        batch[m + i__ * batch_dim1] = temp2[k] * int2dz[k +
                                                                        ze *
                                                                        int2dz_dim1]
                            + temp2[k + 1] * int2dz[k + 1 +
                                                    ze * int2dz_dim1] +
                            temp2[k + 2] * int2dz[k + 2 + ze * int2dz_dim1] +
                            temp2[k + 3] * int2dz[k + 3 + ze * int2dz_dim1] +
                            temp2[k + 4] * int2dz[k + 4 + ze * int2dz_dim1] +
                            temp2[k + 5] * int2dz[k + 5 + ze * int2dz_dim1] +
                            temp2[k + 6] * int2dz[k + 6 + ze * int2dz_dim1] +
                            temp2[k + 7] * int2dz[k + 7 + ze * int2dz_dim1] +
                            temp2[k + 8] * int2dz[k + 8 + ze * int2dz_dim1];
                        k += 9;
                    }
                }
                i__ = i__ - nxyze + xe;
                nxyze = nxyze - se - 1;
/* L920: */
            }
/* L910: */
        }
/* L900: */
    }
    return 0;


/*                       ******************** */
/*                       *  # of roots > 9  * */
/*                       ******************** */

/*             ...outer loops over x-contributions. No skipping of the */
/*                x-contributions of 0-type can be done here, since */
/*                the 2DX integrals carry the Rys weight! */


  L10:
    xep = *nxyzet + 3;
    i__1 = *shellp;
    for (xe = 0; xe <= i__1; ++xe)
    {
        xep = xep + xe - 2;
        xemax = xe * *shellp;
        yeend = *shellp - xe;
        i__2 = *ngqexq;
        for (m = 1; m <= i__2; ++m)
        {
            temp1[m] = scale[m] * int2dx[m + xe * int2dx_dim1];
        }


/*             ...middle loops over y-contributions. Skip multiplication */
/*                of y-contributions, if we have a 0-type, as then the */
/*                2DY integral is equal to 1. */


        xyep = xep - xemax;
        i__2 = yeend;
        for (ye = 0; ye <= i__2; ++ye)
        {
            xye = xe + ye;
            --xyep;
            seend = max (*shella, xye);
            if (ye == 0)
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n];
                }
            }
            else
            {
                i__3 = *ngqexq;
                for (n = 1; n <= i__3; ++n)
                {
                    temp2[n] = temp1[n] * int2dy[n + ye * int2dy_dim1];
                }
            }


/*             ...inner loops over E-pairs. Skip multiplication */
/*                of z-contributions, if we have a 0-type, as */
/*                then the 2DZ integral is equal to 1. */


            i__ = xyep;
            nxyze = *nxyzp;
            i__3 = seend;
            for (se = *shellp; se >= i__3; --se)
            {
                ze = se - xye;


/*             ...all info concerning all x-, y- and z-contributions */
/*                have been collected for all exponent quadruplets at */
/*                once. Sum up the 2D X,Y,Z integral products to the */
/*                appropriate place of the batch. */


                if (ze == 0)
                {
                    k = 0;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        sum = 0.;
                        i__5 = *ngqp;
                        for (n = 1; n <= i__5; ++n)
                        {
                            sum += temp2[k + n];
                        }
                        k += *ngqp;
                        batch[m + i__ * batch_dim1] = sum;
                    }
                }
                else
                {
                    k = 0;
                    i__4 = *nexq;
                    for (m = 1; m <= i__4; ++m)
                    {
                        sum = 0.;
                        i__5 = *ngqp;
                        for (n = 1; n <= i__5; ++n)
                        {
                            sum += temp2[k + n] * int2dz[k + n + ze *
                                                         int2dz_dim1];
                        }
                        k += *ngqp;
                        batch[m + i__ * batch_dim1] = sum;
                    }
                }


/*             ...next z-contribution. */


                i__ = i__ - nxyze + xe;
                nxyze = nxyze - se - 1;
/* L1020: */
            }


/*             ...next y- and x-contribution. */


/* L1010: */
        }
/* L1000: */
    }


/*             ...ready! */


    return 0;
}                               /* erd__int2d_to_e000__ */
