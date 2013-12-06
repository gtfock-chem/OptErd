#include <stdio.h>
#include <stdlib.h>

#include "erd.h"


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__MEMORY_ERI_BATCH */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__MEMORY_1111_CSGTO */
/*                ERD__MEMORY_CSGTO */
/*  DESCRIPTION : Main operation that calculates the minimum and optimum */
/*                int/flp memory needed for evaluating a batch */
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
/*                    IMIN,IOPT    =  minimum/optimum int memory */
/*                    ZMIN,ZOPT    =  minimum/optimum flp memory */
/* ------------------------------------------------------------------------ */
int erd__memory_eri_batch_ (int * nalpha, int * ncoeff,
                            int * ncgto1, int * ncgto2, int * ncgto3,
                            int * ncgto4, int * npgto1, int * npgto2,
                            int * npgto3, int * npgto4, int * shell1,
                            int * shell2, int * shell3, int * shell4,
                            double * x1, double * y1, double * z1,
                            double * x2, double * y2, double * z2,
                            double * x3, double * y3, double * z3,
                            double * x4, double * y4, double * z4,
                            double * alpha, double * cc,
                            int * spheric, int * imin, int * iopt,
                            int * zmin, int * zopt)
{
    int maxshell;
    
    maxshell = MAX(*shell1, *shell2);
    maxshell = MAX(maxshell, *shell3);
    if (MAX(maxshell, *shell4) < 2)
    {
        erd__memory_1111_csgto (*npgto1, *npgto2, *npgto3, *npgto4,
                                *shell1, *shell2, *shell3, *shell4,
                                *x1, *y1, *z1, *x2, *y2, *z2,
                                *x3, *y3, *z3, *x4, *y4, *z4,
                                alpha, cc,
                                imin, iopt, zmin, zopt);
    }
    else
    {
        erd__memory_csgto (*npgto1, *npgto2, *npgto3, *npgto4,
                           *shell1, *shell2, *shell3, *shell4,
                           *x1, *y1, *z1, *x2, *y2, *z2,
                           *x3, *y3, *z3, *x4, *y4, *z4,
                           alpha, cc, *spheric, imin, iopt, zmin, zopt);
    }

    return 0;
}
