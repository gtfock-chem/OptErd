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
ERD_OFFLOAD void erd__memory_csgto(uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
    uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    bool spheric, size_t *iopt, size_t *zopt)
{
    uint32_t nxyzhrr;
    uint32_t nrya, nryb, nryc, nryd;
    uint32_t nxyza, nxyzb, nxyzc, nxyzd, shella, shellb, shellc, shelld, npgtoa, npgtob,
        npgtoc, npgtod, mnprim, dummyi[15];
    int mxprim;
    double dummyr[12];
    uint32_t mxsize, nxyzet, nxyzft;
    uint32_t ncolhrr, nrothrr;

    *iopt = 0;
    *zopt = 0;
/*             ...fix the A,B,C,D labels from the 1,2,3,4 ones. */
/*                Calculate the relevant data for the A,B,C,D batch of */
/*                integrals. Most of this data is not needed to */
/*                evaluate the memory requirements and is dumped into */
/*                the DUMMYx arrays with x=I,R,L standing for int, */
/*                real and int, respectively. */ 
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
        return;
    }

/*             ...simulate the cartesian contracted (e0|f0) batch */
/*                generation. */
    const uint32_t npgtoab = npgtoa * npgtob;
    const uint32_t npgtocd = npgtoc * npgtod;
    const uint32_t nxyzt = nxyzet * nxyzft;

/*             ...at this point we would determine the IJ and KL */
/*                exponent pairs necessay to evaluate the cartesian */
/*                contracted (e0|f0) batch after a possible screening */
/*                of the primitives. Since at the moment we do not */
/*                apply screening for memory determination, we use the */
/*                complete set of IJ and KL pairs. In any event this */
/*                will be changed in the future, this is where a memory */
/*                routine handling the IJ and KL pair determination */
/*                should be placed. */
    const uint32_t nij = npgtoab;
    const uint32_t nkl = npgtocd;
    uint32_t zneed = PAD_LEN(npgtoab) + PAD_LEN(npgtocd);
    uint32_t ineed = 2 * PAD_LEN2(npgtoab) + 2 * PAD_LEN2(npgtocd);
    *iopt = MAX (*iopt, ineed);
    *zopt = MAX (*zopt, zneed);


/*             ...determine minimum and optimum flp needs for the */
/*                unnormalized cartesian (e0|f0) contracted batch */
/*                generation. */
    const uint32_t shellp = shella + shellb;
    const uint32_t shellq = shellc + shelld;
    const uint32_t shellt = shellp + shellq;
    const uint32_t ngqp = shellt / 2 + 1;
    const uint32_t nmom = (ngqp << 1) - 1;
    const uint32_t ngqscr = nmom * 5 + (ngqp << 1) - 2;
    const uint32_t mxshell = max4x32u(shell1, shell2, shell3, shell4);

    int zout2;
    erd__e0f0_def_blocks(0, npgtoa, npgtob, npgtoc, npgtod,
                          shellp, shellq, nij, nkl,
                          ngqp, ngqscr, nxyzt,
                          1, &zout2,
                          NULL);

    mnprim = MAX (MIN (npgtoa, npgtob), MIN (npgtoc, npgtod));    
    mxprim = MAX (npgtoa, npgtob);
    mxprim = MAX (mxprim, npgtoc);
    mxprim = MAX (mxprim, npgtod); 
    
    ineed = ineed + (mxprim << 1) + mnprim;
    *iopt = MAX (*iopt, ineed);
    *zopt = MAX (*zopt, zout2);

/*             ...determine the int/flp memory needs for the next */
/*                steps: */
/*                1) expanding the contraction indices (if any) */
/*                2) reordering the contraction indices (if any) */
/*                3) transposing the contraction indices (if any) */
/*                4) HRR contraction */
/*                5) cartesian -> spherical transformation or */
/*                   cartesian normalization */
/*                The space partitioning of the flp array will be */
/*                as follows: */
/*                 |  Zone 1  |  Zone 2  |  Zone 3  |  Zone 4  | */
/*                 Zone 1 and 2:  2 batches of HRR maximum size */
/*                       Zone 3:  cart -> spher transformation data */
/*                                             or */
/*                                cartesian normalization factors */
/*                       Zone 4:  HRR contraction data */
    mxsize = nxyzhrr;


/*             ...memory for Zone 1 and 2. */
    ineed = 0;
    zneed = mxsize + mxsize;


/*             ...memory for Zone 3. */
    if (spheric) {
        if (mxshell > 1) {
            if (shelld > 1) {
                const uint32_t nrowd = (shelld / 2 + 1) * (shelld / 2 + 2) / 2;
                const uint32_t nrotd = nrowd * nryd;
                ineed = ineed + nryd + nrotd;
                zneed = zneed + nrotd + nxyzd;
            }
            if (shellc > 1 && shellc != shelld) {
                const uint32_t nrowc = (shellc / 2 + 1) * (shellc / 2 + 2) / 2;
                const uint32_t nrotc = nrowc * nryc;
                ineed = ineed + nryc + nrotc;
                zneed = zneed + nrotc + nxyzc;
            }
            if (shellb > 1 && shellb != shellc && shellb != shelld) {
                const uint32_t nrowb = (shellb / 2 + 1) * (shellb / 2 + 2) / 2;
                const uint32_t nrotb = nrowb * nryb;
                ineed = ineed + nryb + nrotb;
                zneed = zneed + nrotb + nxyzb;
            }
            if (shella > 1 && shella != shellb && shella != shellc && shella != shelld) {
                const uint32_t nrowa = (shella / 2 + 1) * (shella / 2 + 2) / 2;
                const uint32_t nrota = nrowa * nrya;
                ineed = ineed + nrya + nrota;
                zneed = zneed + nrota + nxyza;
            }
        }
    } else {
        if (mxshell > 1) {
            zneed = zneed + mxshell + 1;
        }
    }


/*             ...memory for Zone 4. */
    ineed = ineed + (ncolhrr << 2) + (nrothrr << 1);
    zneed += nrothrr << 1;
    *iopt = MAX (*iopt, ineed);
    *zopt = MAX (*zopt, zneed);
}
