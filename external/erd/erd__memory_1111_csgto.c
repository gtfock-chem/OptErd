/* erd__memory_1111_csgto.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static integer c__0 = 0;

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
erd__memory_1111_csgto_ (integer * nalpha, integer *
                          ncoeff, integer * ncgto1, integer * ncgto2,
                          integer * ncgto3, integer * ncgto4,
                          integer * npgto1, integer * npgto2,
                          integer * npgto3, integer * npgto4,
                          integer * shell1, integer * shell2,
                          integer * shell3, integer * shell4, doublereal * x1,
                          doublereal * y1, doublereal * z1, doublereal * x2,
                          doublereal * y2, doublereal * z2, doublereal * x3,
                          doublereal * y3, doublereal * z3, doublereal * x4,
                          doublereal * y4, doublereal * z4,
                          doublereal * alpha, doublereal * cc,
                          integer * l1cache, integer * nctrow, integer * imin,
                          integer * iopt, integer * zmin, integer * zopt)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l;
    extern /* Subroutine */ int erd__1111_def_blocks_ (integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, logical *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *,
                                                        integer *, integer *);
    static integer nij, nkl, lcc1, lcc2, lcc3, lcc4, nctr, lexp1, lexp2,
        lexp3, lexp4, zout1, zout2, ineed, nxyz1, nxyz2, nxyz3, nxyz4;
    static logical atom12, atom23;
    static integer zneed;
    static logical atom34;
    static integer nxyzt;
    static logical equal12, atomic;
    static integer ncgto12;
    static logical equal34;
    static integer ncgto34, shellp, npgto12, shellt, npgto34, mnprim,
        dummyi[22];
    static logical memory;
    static integer mxprim;

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__MEMORY_1111_CSGTO */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__MEMORY_1111_BLOCKS */
/*  DESCRIPTION : This operation calculates the minimum and optimum */
/*                integer/flp memory needed for evaluating a batch */
/*                of contracted electron repulsion integrals on up to */
/*                four different centers involving s- and p-type shells */
/*                only! */


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
/*                    L1CACHE      =  Size of level 1 cache in units of */
/*                                    8 Byte */
/*                    NCTROW       =  minimum # of rows that are */
/*                                    accepted for blocked contractions */

/*                  Output: */

/*                    IMIN,IOPT    =  minimum/optimum integer memory */
/*                    ZMIN,ZOPT    =  minimum/optimum flp memory */


/*  AUTHOR      : Norbert Flocke */
/* ------------------------------------------------------------------------ */


/*             ...include files and declare variables. */




/* ------------------------------------------------------------------------ */


/*             ...set initial memory values. */


    /* Parameter adjustments */
    --alpha;
    --cc;

    /* Function Body */
    *imin = 0;
    *iopt = 0;
    *zmin = 0;
    *zopt = 0;


/*             ...simulate the cartesian contracted (12|34) batch */
/*                generation. */


    shellp = *shell1 + *shell2;
    shellt = shellp + *shell3 + *shell4;
    atom12 = *x1 == *x2 && *y1 == *y2 && *z1 == *z2;
    atom23 = *x2 == *x3 && *y2 == *y3 && *z2 == *z3;
    atom34 = *x3 == *x4 && *y3 == *y4 && *z3 == *z4;
    atomic = atom12 && atom34 && atom23;
    if (atomic && shellt % 2 == 1)
    {
        return 0;
    }
    lexp1 = 1;
    lexp2 = lexp1 + *npgto1;
    lexp3 = lexp2 + *npgto2;
    lexp4 = lexp3 + *npgto3;
    lcc1 = 1;
    lcc2 = lcc1 + *npgto1 * *ncgto1;
    lcc3 = lcc2 + *npgto2 * *ncgto2;
    lcc4 = lcc3 + *npgto3 * *ncgto3;


/*             ...determine csh equality between center pairs 1,2 */
/*                and 3,4 in increasing order of complexity: */

/*                 centers -> shells -> exponents -> ctr coefficients */


    equal12 = atom12;
    if (equal12)
    {
        equal12 = *shell1 == *shell2 && *npgto1 == *npgto2
            && *ncgto1 == *ncgto2;
        if (equal12)
        {
            k = lexp1 - 1;
            l = lexp2 - 1;
            i__1 = *npgto1;
            for (i__ = 1; i__ <= i__1; ++i__)
            {
                equal12 = equal12 && alpha[k + i__] == alpha[l + i__];
/* L100: */
            }
            if (equal12)
            {
                k = lcc1 - 1;
                l = lcc2 - 1;
                i__1 = *ncgto1;
                for (j = 1; j <= i__1; ++j)
                {
                    if (equal12)
                    {
                        i__2 = *npgto1;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            equal12 = equal12 && cc[k + i__] == cc[l + i__];
/* L120: */
                        }
                        k += *npgto1;
                        l += *npgto1;
                    }
/* L110: */
                }
            }
        }
    }
    equal34 = atom34;
    if (equal34)
    {
        equal34 = *shell3 == *shell4 && *npgto3 == *npgto4
            && *ncgto3 == *ncgto4;
        if (equal34)
        {
            k = lexp3 - 1;
            l = lexp4 - 1;
            i__1 = *npgto3;
            for (i__ = 1; i__ <= i__1; ++i__)
            {
                equal34 = equal34 && alpha[k + i__] == alpha[l + i__];
/* L130: */
            }
            if (equal34)
            {
                k = lcc3 - 1;
                l = lcc4 - 1;
                i__1 = *ncgto3;
                for (j = 1; j <= i__1; ++j)
                {
                    if (equal34)
                    {
                        i__2 = *npgto3;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            equal34 = equal34 && cc[k + i__] == cc[l + i__];
/* L150: */
                        }
                        k += *npgto3;
                        l += *npgto3;
                    }
/* L140: */
                }
            }
        }
    }


/*             ...calculate relevant data for the [12|34] batch of */
/*                integrals, such as dimensions, total # of integrals */
/*                to be expected, etc... */


    nxyz1 = *shell1 + *shell1 + 1;
    nxyz2 = *shell2 + *shell2 + 1;
    nxyz3 = *shell3 + *shell3 + 1;
    nxyz4 = *shell4 + *shell4 + 1;
    nxyzt = nxyz1 * nxyz2 * nxyz3 * nxyz4;
    if (equal12)
    {
        npgto12 = *npgto1 * (*npgto1 + 1) / 2;
        ncgto12 = *ncgto1 * (*ncgto1 + 1) / 2;
    }
    else
    {
        npgto12 = *npgto1 * *npgto2;
        ncgto12 = *ncgto1 * *ncgto2;
    }
    if (equal34)
    {
        npgto34 = *npgto3 * (*npgto3 + 1) / 2;
        ncgto34 = *ncgto3 * (*ncgto3 + 1) / 2;
    }
    else
    {
        npgto34 = *npgto3 * *npgto4;
        ncgto34 = *ncgto3 * *ncgto4;
    }
    nctr = ncgto12 * ncgto34;


/*             ...at this point we would determine the IJ and KL */
/*                exponent pairs necessay to evaluate the cartesian */
/*                contracted (12|34) batch after a possible screening */
/*                of the primitives. Since at the moment we do not */
/*                apply screening for memory determination, we use the */
/*                complete set of IJ and KL pairs. In any event this */
/*                will be changed in the future, this is where a memory */
/*                routine handling the IJ and KL pair determination */
/*                should be placed. */


    nij = npgto12;
    nkl = npgto34;
    zneed = npgto12 + npgto34;
    ineed = zneed << 1;
    *imin = max (*imin, ineed);
    *iopt = max (*iopt, ineed);
    *zmin = max (*zmin, zneed);
    *zopt = max (*zopt, zneed);


/*             ...determine minimum and optimum flp needs for the */
/*                unnormalized cartesian (12|34) contracted batch */
/*                generation. */


    memory = TRUE_;
    erd__1111_def_blocks_ (&c__0, npgto1, npgto2, npgto3, npgto4, &nij, &nkl,
                            &ncgto12, &ncgto34, &nctr, &nxyzt, l1cache,
                            nctrow, &memory, &zout1, &zout2, dummyi,
                            &dummyi[1], &dummyi[2], &mxprim, &mnprim,
                            &dummyi[3], &dummyi[4], &dummyi[5], &dummyi[6],
                            &dummyi[7], &dummyi[8], &dummyi[9], &dummyi[10],
                            &dummyi[11], &dummyi[12], &dummyi[13],
                            &dummyi[14], &dummyi[15], &dummyi[16],
                            &dummyi[17], &dummyi[18], &dummyi[19],
                            &dummyi[20], &dummyi[21]);
    ineed = ineed + (mxprim << 1) + mnprim;
    *imin = max (*imin, ineed);
    *iopt = max (*iopt, ineed);
    *zmin = max (*zmin, zout1);
    *zopt = max (*zopt, zout2);


/*             ...determine the integer/flp memory needs for the next */
/*                steps: */

/*                1) expanding the contraction indices (if any) */
/*                2) reordering the contraction indices (if any) */

/*                The space partitioning of the flp array will be */
/*                as follows: */


/*                         |  Zone 1  |  Zone 2  | */


/*                 Zone 1 and 2:  2 batches of final (12|34) size */


    nctr = *ncgto1 * *ncgto2 * *ncgto3 * *ncgto4;
    zneed = (nctr << 1) * nxyzt;
    *zmin = max (*zmin, zneed);
    *zopt = max (*zopt, zneed);


/*             ...ready! */


    return 0;
}                               /* erd__memory_1111_csgto__ */
