/* erd__memory_eri_batch.f -- translated by f2c (version 20100827).
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

/* Common Block Declarations */


/* Table of constant values */

static integer c__32768 = 32768;
static integer c__110 = 110;

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
erd__memory_eri_batch_ (integer * nalpha, integer * ncoeff,
                         integer * ncgto1, integer * ncgto2, integer * ncgto3,
                         integer * ncgto4, integer * npgto1, integer * npgto2,
                         integer * npgto3, integer * npgto4, integer * shell1,
                         integer * shell2, integer * shell3, integer * shell4,
                         doublereal * x1, doublereal * y1, doublereal * z1,
                         doublereal * x2, doublereal * y2, doublereal * z2,
                         doublereal * x3, doublereal * y3, doublereal * z3,
                         doublereal * x4, doublereal * y4, doublereal * z4,
                         doublereal * alpha, doublereal * cc,
                         logical * spheric, integer * imin, integer * iopt,
                         integer * zmin, integer * zopt)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int erd__memory_csgto_ (integer *, integer *,
                                                     integer *, integer *,
                                                     integer *, integer *,
                                                     integer *, integer *,
                                                     integer *, integer *,
                                                     integer *, integer *,
                                                     integer *, integer *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *,
                                                     doublereal *, integer *,
                                                     integer *, logical *,
                                                     integer *, integer *,
                                                     integer *, integer *),
        erd__memory_1111_csgto_ (integer *, integer *, integer *, integer *,
                                  integer *, integer *, integer *, integer *,
                                  integer *, integer *, integer *, integer *,
                                  integer *, integer *, doublereal *,
                                  doublereal *, doublereal *, doublereal *,
                                  doublereal *, doublereal *, doublereal *,
                                  doublereal *, doublereal *, doublereal *,
                                  doublereal *, doublereal *, doublereal *,
                                  doublereal *, integer *, integer *,
                                  integer *, integer *, integer *, integer *);

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__MEMORY_ERI_BATCH */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__MEMORY_1111_CSGTO */
/*                ERD__MEMORY_CSGTO */
/*  DESCRIPTION : Main operation that calculates the minimum and optimum */
/*                integer/flp memory needed for evaluating a batch */
/*                of contracted electron repulsion integrals on up to */
/*                four different centers between cartesian or spherical */
/*                gaussian type shells. */


/*                  Input (x = 1,2,3 and 4): */

/*                    NALPHA       =  total # of exponents */
/*                    NCOEFF       =  total # of contraction coeffs */
/*                    NCGTOx       =  # of contractions for csh x */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for csh x */
/*                    SHELLx       =  the shell type for csh x */
/*                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers */
/*                                    y = 1,2,3 and 4 */
/*                    ALPHA        =  primitive exponents for csh */
/*                                    1,2,3,4 in that order */
/*                    CC           =  contraction coefficient for csh */
/*                                    1,2,3,4 in that order, for each */
/*                                    csh individually such that an */
/*                                    (I,J) element corresponds to the */
/*                                    I-th primitive and J-th contraction. */
/*                    SPHERIC      =  is true, if spherical integrals */
/*                                    are wanted, false if cartesian */
/*                                    ones are wanted */


/*                  Output: */

/*                    IMIN,IOPT    =  minimum/optimum integer memory */
/*                    ZMIN,ZOPT    =  minimum/optimum flp memory */


/*  AUTHOR      : Norbert Flocke */
/* ------------------------------------------------------------------------ */

/*             ...include files and declare variables. */


/* ------------------------------------------------------------------------ */
/*  INCLUDE FILE: ERD__TUNING */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  DESCRIPTION : Include file used to set all tuning parameters. */
/*                For optimum performance of the present two electron */
/*                integral package, there are 3 parameters that are */
/*                set and which allow for tuning of performance: */


/*                  L1CACHE   =   Size of level 1 cache in units of */
/*                                8 Byte */

/*                  TILE      =   Number of rows and columns in units */
/*                                of 8 Byte of level 1 cache square */
/*                                tile array used for performing optimum */
/*                                out-of-place matrix transposition */
/*                                operations */

/*                  NCTROW    =   Minimum number of rows that are */
/*                                still accepted for block contractions. */
/*                                This quantity controls extremly low */
/*                                # of row values due to small level */
/*                                1 cache size. */

/*                The value for TILE can be estimated as follows. Assume */
/*                the level 1 cache is 3/4 free. We have to fit 2 square */
/*                matrices (out-of-place transpositions!) of tile size */
/*                into it, giving: */

/*                       TILE = Int ( SQRT ((3/4)*L1CACHE/2) */


/*  AUTHOR      : Norbert Flocke */
/* ------------------------------------------------------------------------ */


/*             ...declare and set the tuning parameters. */




/* ------------------------------------------------------------------------ */


/*             ...call special memory routine for only s- and p-type */
/*                integrals. */


    /* Parameter adjustments */
    --alpha;
    --cc;

    /* Function Body */
/* Computing MAX */
    i__1 = max (*shell1, *shell2), i__1 = max (i__1, *shell3);
    if (max (i__1, *shell4) < 2)
    {
        erd__memory_1111_csgto_ (nalpha, ncoeff, ncgto1, ncgto2, ncgto3,
                                  ncgto4, npgto1, npgto2, npgto3, npgto4,
                                  shell1, shell2, shell3, shell4, x1, y1, z1,
                                  x2, y2, z2, x3, y3, z3, x4, y4, z4,
                                  &alpha[1], &cc[1], &c__32768, &c__110, imin,
                                  iopt, zmin, zopt);
    }
    else
    {
        erd__memory_csgto_ (nalpha, ncoeff, ncgto1, ncgto2, ncgto3, ncgto4,
                             npgto1, npgto2, npgto3, npgto4, shell1, shell2,
                             shell3, shell4, x1, y1, z1, x2, y2, z2, x3, y3,
                             z3, x4, y4, z4, &alpha[1], &cc[1], &c__32768,
                             &c__110, spheric, imin, iopt, zmin, zopt);
    }


/*             ...ready! */


    return 0;
}                               /* erd__memory_eri_batch__ */
