#include <stdio.h>
#include <stdlib.h>

#include "erd.h"


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__MEMORY_1111_CSGTO */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__MEMORY_1111_BLOCKS */
/*  DESCRIPTION : This operation calculates the minimum and optimum */
/*                int/flp memory needed for evaluating a batch */
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
/*                    IMIN,IOPT    =  minimum/optimum int memory */
/*                    ZMIN,ZOPT    =  minimum/optimum flp memory */
/* ------------------------------------------------------------------------ */
int erd__memory_1111_csgto (int npgto1, int npgto2,
                            int npgto3, int npgto4,
                            int shell1, int shell2,
                            int shell3, int shell4,
                            double x1, double y1, double z1,
                            double x2, double y2, double z2,
                            double x3, double y3, double z3,
                            double x4, double y4, double z4,
                            double *alpha, double *cc,
                            int *imin, int *iopt, int *zmin, int *zopt)
{
    int i;
    int lcc1, lcc2, lcc3, lcc4, lexp1, lexp2,
        lexp3, lexp4, zout2, ineed, nxyz1, nxyz2, nxyz3, nxyz4;
    int atom12, atom23;
    int zneed;
    int atom34;
    int nxyzt;
    int equal12, atomic;
    int equal34;
    int shellp, npgto12, shellt, npgto34;
    int mxprim;
    int mnprim;
    int npminrs;
    int npmintu;
    
    *iopt = 0;
    *zopt = 0;
   
    shellp = shell1 + shell2;
    shellt = shellp + shell3 + shell4;
    atom12 = ((x1 == x2) && (y1 == y2) && (z1 == z2));
    atom23 = ((x2 == x3) && (y2 == y3) && (z2 == z3));
    atom34 = ((x3 == x4) && (y3 == y4) && (z3 == z4));
    atomic = atom12 && atom34 && atom23;
    if (atomic && shellt % 2 == 1)
    {
        return 0;
    }
    
    lexp1 = 0;
    lexp2 = lexp1 + npgto1;
    lexp3 = lexp2 + npgto2;
    lexp4 = lexp3 + npgto3;
    lcc1 = 0;
    lcc2 = lcc1 + npgto1;
    lcc3 = lcc2 + npgto2;
    lcc4 = lcc3 + npgto3;

/*             ...determine csh equality between center pairs 1,2 */
/*                and 3,4 in increasing order of complexity: */
/*                 centers -> shells -> exponents -> ctr coefficients */
    equal12 = atom12;
    if (equal12)
    {
        equal12 = ((shell1 == shell2) &&
                   (npgto1 == npgto2));
        if (equal12)
        {
            for (i = 0; i < npgto1; ++i)
            {
                equal12 = (equal12 &&
                           (alpha[lexp1 + i] == alpha[lexp2 + i]));
            }
            if (equal12)
            {
                for (i = 0; i < npgto1; ++i)
                {
                    equal12 = (equal12 && 
                        (cc[lcc1 + i] == cc[lcc2 + i]));
                }
            }
        }
    }
    equal34 = atom34;
    if (equal34)
    {
        equal34 = ((shell3 == shell4) &&
                   (npgto3 == npgto4));
        if (equal34)
        {
            for (i = 0; i < npgto3; ++i)
            {
                equal34 = (equal34 &&
                           (alpha[lexp3 + i] == alpha[lexp4 + i]));
            }
            if (equal34)
            {
                for (i = 0; i < npgto3; ++i)
                {
                    equal34 = (equal34 &&
                        (cc[lcc3 + i] == cc[lcc4 + i]));
                }
            }
        }
    }


/*             ...calculate relevant data for the [12|34] batch of */
/*                integrals, such as dimensions, total # of integrals */
/*                to be expected, etc... */
    nxyz1 = shell1 + shell1 + 1;
    nxyz2 = shell2 + shell2 + 1;
    nxyz3 = shell3 + shell3 + 1;
    nxyz4 = shell4 + shell4 + 1;
    nxyzt = nxyz1 * nxyz2 * nxyz3 * nxyz4;
    if (equal12)
    {
        npgto12 = npgto1 * (npgto1 + 1) / 2;
    }
    else
    {
        npgto12 = npgto1 * npgto2;
    }
    if (equal34)
    {
        npgto34 = npgto3 * (npgto3 + 1) / 2;
    }
    else
    {
        npgto34 = npgto3 * npgto4;
    }

/*             ...at this point we would determine the IJ and KL */
/*                exponent pairs necessay to evaluate the cartesian */
/*                contracted (12|34) batch after a possible screening */
/*                of the primitives. Since at the moment we do not */
/*                apply screening for memory determination, we use the */
/*                complete set of IJ and KL pairs. In any event this */
/*                will be changed in the future, this is where a memory */
/*                routine handling the IJ and KL pair determination */
/*                should be placed. */
    zneed = npgto12 + npgto34;
    ineed = 2 * zneed;
    *iopt = MAX(*iopt, ineed);
    *zopt = MAX(*zopt, zneed);


/*             ...determine minimum and optimum flp needs for the */
/*                unnormalized cartesian (12|34) contracted batch */
/*                generation. */
    erd__1111_def_blocks (0, npgto1, npgto2, npgto3, npgto4,
                          npgto12, npgto34, nxyzt, 1,
                          &zout2, NULL, NULL, NULL, NULL,
                          NULL, NULL, NULL, NULL, NULL,
                          NULL, NULL, NULL, NULL, NULL,
                          NULL, NULL, NULL, NULL, NULL);

    mxprim = MAX(npgto1, npgto2);
    mxprim = MAX(mxprim, npgto3);
    mxprim = MAX(mxprim, npgto4);

    npminrs = MIN(npgto1, npgto2);
    npmintu = MIN(npgto1, npgto2);         
    mnprim = MAX(npminrs, npmintu);
    
    ineed = ineed + 2 * mxprim + mnprim;
    *iopt = MAX(*iopt, ineed);
    *zopt = MAX(*zopt, zout2);


/*             ...determine the int/flp memory needs for the next */
/*                steps: */
/*                1) expanding the contraction indices (if any) */
/*                2) reordering the contraction indices (if any) */
/*                The space partitioning of the flp array will be */
/*                as follows: */
/*                         |  Zone 1  |  Zone 2  | */
/*                 Zone 1 and 2:  2 batches of final (12|34) size */
    zneed = 2 * nxyzt;
    *zopt = MAX(*zopt, zneed);

    *zmin = *zopt;
    *imin = *iopt;

    return 0;
}
