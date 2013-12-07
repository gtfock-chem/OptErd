/* erd__rys_roots_weights.f -- translated by f2c (version 20100827).
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
erd__rys_roots_weights_ (int * nt, int * ntgqp,
                          int * ngqp, int * nmom, double * tval,
                          double * ryszero, double * ftable,
                          int * mgrid, int * ngrid, double * tmax,
                          double * tstep, double * tvstep,
                          double * a, double * b, double * mom,
                          double * dia, double * off,
                          double * row1, double * row2,
                          double * rts, double * wts)
{
    int ftable_dim1, ftable_offset, i__1;

    int n;
    double t;
    double delta;
    int tgrid;

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__RYS_ROOTS_WEIGHTS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__RYS_1_ROOTS_WEIGHTS */
/*                ERD__RYS_2_ROOTS_WEIGHTS */
/*                ERD__RYS_3_ROOTS_WEIGHTS */
/*                ERD__RYS_4_ROOTS_WEIGHTS */
/*                ERD__RYS_5_ROOTS_WEIGHTS */
/*                ERD__RYS_X_ROOTS_WEIGHTS */
/*  DESCRIPTION : This routine calculates NGQP-point Gaussian quadrature */
/*                rules on [0,1] over the Rys weight functions: */


/*                                       exp(-T*x) */
/*                             W   (x) = --------- */
/*                              Rys      2*sqrt(x) */


/*                for a set of NT T-exponents. Special interpolation */
/*                routines are provided for low number of roots and */
/*                weigths (NGQP < 6). On exit, NT x NGQP = NTGQP roots */
/*                and weights have been produced. */


/*                  Input: */

/*                    NT           =  # of T-exponents */
/*                    NTGQP        =  # of roots times # of T-exponents */
/*                    NGQP         =  # of gaussian quadrature points */
/*                                    (roots) */
/*                    NMOM         =  # of necessary moment integrals */
/*                                    to calculate the quadrature roots */
/*                    TVAL         =  the T-exponents */
/*                    RYSZERO      =  will hold the zeroth Rys moments */
/*                                    for all T-exponents */
/*                    FTABLE       =  Fm (T) table for interpolation */
/*                                    in low T region */
/*                    MGRID        =  maximum m in Fm (T) table */
/*                    NGRID        =  # of T's for which Fm (T) table */
/*                                    was set up */
/*                    TMAX         =  maximum T in Fm (T) table */
/*                    TSTEP        =  difference between two consecutive */
/*                                    T's in Fm (T) table */
/*                    TVSTEP       =  Inverse of TSTEP */
/*                    A,B          =  will contain the recurrence */
/*                                    coefficients for the auxillary */
/*                                    polynomials */
/*                    MOM          =  will contain the normed auxillary */
/*                                    polynomial modified moments */
/*                    DIA,OFF      =  will contain the diagonal and */
/*                                    offdiagonal elements of the */
/*                                    tridiagonal symmetric terminal */
/*                                    matrix */
/*                    ROW1,ROW2    =  first,second row intermediates. */
/*                                    Will be used to evaluate the */
/*                                    tridiagonal elements of the */
/*                                    symmetric terminal matrix in an */
/*                                    efficient way using Sack and */
/*                                    Donovan's method */

/*                  Output: */

/*                    RTS          =  the roots array */
/*                    WTS          =  the weights array */


/*  AUTHOR      : Norbert Flocke */
/* ------------------------------------------------------------------------ */


/*             ...include files and declare variables. */




/* ------------------------------------------------------------------------ */


/*             ...proceed according to the number of roots and weights */
/*                required. */


    /* Parameter adjustments */
    --ryszero;
    --tval;
    --wts;
    --rts;
    --off;
    --dia;
    --row2;
    --row1;
    --mom;
    b -= 2;
    --a;
    ftable_dim1 = *mgrid - 0 + 1;
    ftable_offset = 0 + ftable_dim1 * 0;
    ftable -= ftable_offset;

    /* Function Body */
    switch (min (*ngqp, 6))
    {
    case 1:
        goto L1000;
    case 2:
        goto L2000;
    case 3:
        goto L3000;
    case 4:
        goto L4000;
    case 5:
        goto L5000;
    case 6:
        goto L6000;
    }


/*             ...# of roots and weights < 6 */


  L1000:
    erd__rys_1_roots_weights_ (nt, &tval[1], &rts[1], &wts[1]);
    return 0;
  L2000:
    erd__rys_2_roots_weights_ (nt, ntgqp, &tval[1], &rts[1], &wts[1]);
    return 0;
  L3000:
    erd__rys_3_roots_weights_ (nt, ntgqp, &tval[1], &rts[1], &wts[1]);
    return 0;
  L4000:
    erd__rys_4_roots_weights_ (nt, ntgqp, &tval[1], &rts[1], &wts[1]);
    return 0;
  L5000:
    erd__rys_5_roots_weights_ (nt, ntgqp, &tval[1], &rts[1], &wts[1]);
    return 0;


/*             ...# of roots and weights >= 6. Accumulate all zeroth */
/*                Rys moments and call the general routine. */


  L6000:
    i__1 = *nt;
    for (n = 1; n <= i__1; ++n)
    {
        t = tval[n];
        if (t == 0.)
        {
            ryszero[n] = 1.;
        }
        else if (t <= *tmax)
        {
            tgrid = (int) (t * *tvstep + .5);
            delta = tgrid * *tstep - t;
            ryszero[n] = (((((ftable[tgrid * ftable_dim1 + 6] * delta *
                              .166666666666667 + ftable[tgrid * ftable_dim1 +
                                                        5]) * delta * .2 +
                             ftable[tgrid * ftable_dim1 + 4]) * delta * .25 +
                            ftable[tgrid * ftable_dim1 +
                                   3]) * delta * .333333333333333 +
                           ftable[tgrid * ftable_dim1 + 2]) * delta * .5 +
                          ftable[tgrid * ftable_dim1 + 1]) * delta +
                ftable[tgrid * ftable_dim1];
        }
        else
        {
            ryszero[n] = sqrt (3.141592653589793 / t) * .5;
        }
    }
    erd__rys_x_roots_weights_ (nt, ntgqp, ngqp, nmom, &tval[1], &ryszero[1],
                                &a[1], &b[2], &mom[1], &dia[1], &off[1],
                                &row1[1], &row2[1], &rts[1], &wts[1]);


    return 0;
}
