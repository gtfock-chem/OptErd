#include <stdio.h>
#include <stdlib.h>

#include "erd.h"
#include "erdutil.h"

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__MEMORY_CSGTO */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__SET_ABCD */
/*                ERD__E0F0_DEF_BLOCKS */
/*  DESCRIPTION : This operation calculates the minimum and optimum */
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
/*                    L1CACHE      =  Size of level 1 cache in units of */
/*                                    8 Byte */
/*                    NCTROW       =  minimum # of rows that are */
/*                                    accepted for blocked contractions */
/*                    SPHERIC      =  is true, if spherical integrals */
/*                                    are wanted, false if cartesian */
/*                                    ones are wanted */
/*                  Output: */
/*                    IMIN,IOPT    =  minimum/optimum int memory */
/*                    ZMIN,ZOPT    =  minimum/optimum flp memory */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD size_t erd__memory_csgto(uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
    uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    bool spheric)
{
    uint32_t nxyzhrr;
    uint32_t nrya, nryb, nryc, nryd;
    uint32_t nxyza, nxyzb, nxyzc, nxyzd, shella, shellb, shellc, shelld, npgtoa, npgtob,
        npgtoc, npgtod, mnprim, dummyi[15];
    double dummyr[12];
    uint32_t mxsize, nxyzet, nxyzft;
    uint32_t ncolhrr, nrothrr;

    /*
     * ...fix the A,B,C,D labels from the 1,2,3,4 ones.
     *    Calculate the relevant data for the A,B,C,D batch of integrals.
     *    Most of this data is not needed to evaluate the memory requirements and is dumped into 
     *    the DUMMYx arrays with x=I,R,L standing for int, real and int, respectively.
     */ 
    const bool atomic = false;
    bool empty, tr1234;
    erd__set_abcd(npgto1, npgto2, npgto3, npgto4,
        shell1, shell2, shell3, shell4,
        atomic,
        x1, y1, z1, x2, y2, z2,
        x3, y3, z3, x4, y4, z4, spheric,
        &npgtoa, &npgtob, &npgtoc, &npgtod,
        &shella, &shellb, &shellc, &shelld,
        &dummyr[0], &dummyr[1], &dummyr[2],
        &dummyr[3], &dummyr[4], &dummyr[5],
        &dummyr[6], &dummyr[7], &dummyr[8],
        &dummyr[9], &dummyr[10], &dummyr[11],
        &nxyza, &nxyzb,
        &nxyzc, &nxyzd,
        &nxyzet, &nxyzft,
        &nrya, &nryb, &nryc, &nryd,
        &dummyi[0], &dummyi[1],
        &ncolhrr, &nrothrr,
        &nxyzhrr, &empty, &tr1234);
    if (empty) {
        return 0;
    }

    const uint32_t nxyzt = nxyzet * nxyzft;

    return PAD_LEN(nxyzt) + 2 * nxyzhrr;
}
