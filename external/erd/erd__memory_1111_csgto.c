#include <stdio.h>
#include <stdlib.h>

#include "erd.h"


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

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
void erd__memory_1111_csgto(uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
                            uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
                            double x1, double y1, double z1,
                            double x2, double y2, double z2,
                            double x3, double y3, double z3,
                            double x4, double y4, double z4,
                            size_t *iopt, size_t *zopt)
{
    *iopt = 0;
    *zopt = 0;
   
    const uint32_t shellp = shell1 + shell2;
    const uint32_t shellt = shellp + shell3 + shell4;
    const bool atom12 = ((x1 == x2) && (y1 == y2) && (z1 == z2));
    const bool atom23 = ((x2 == x3) && (y2 == y3) && (z2 == z3));
    const bool atom34 = ((x3 == x4) && (y3 == y4) && (z3 == z4));
    const bool atomic = atom12 && atom34 && atom23;
    if (atomic && shellt % 2 == 1) {
        return;
    }
/*             ...determine csh equality between center pairs 1,2 */
/*                and 3,4 in increasing order of complexity: */
/*                 centers -> shells -> exponents -> ctr coefficients */
/*             ...calculate relevant data for the [12|34] batch of */
/*                integrals, such as dimensions, total # of integrals */
/*                to be expected, etc... */
    const uint32_t nxyz1 = shell1 + shell1 + 1;
    const uint32_t nxyz2 = shell2 + shell2 + 1;
    const uint32_t nxyz3 = shell3 + shell3 + 1;
    const uint32_t nxyz4 = shell4 + shell4 + 1;
    const uint32_t nxyzt = nxyz1 * nxyz2 * nxyz3 * nxyz4;
    const uint32_t npgto12 = npgto1 * npgto2;
    const uint32_t npgto34 = npgto3 * npgto4;

/*             ...at this point we would determine the IJ and KL */
/*                exponent pairs necessay to evaluate the cartesian */
/*                contracted (12|34) batch after a possible screening */
/*                of the primitives. Since at the moment we do not */
/*                apply screening for memory determination, we use the */
/*                complete set of IJ and KL pairs. In any event this */
/*                will be changed in the future, this is where a memory */
/*                routine handling the IJ and KL pair determination */
/*                should be placed. */
    uint32_t zneed = PAD_LEN(npgto12) + PAD_LEN(npgto34);
    uint32_t ineed = 2 * PAD_LEN2(npgto12) + 2 * PAD_LEN2(npgto34);
    *iopt = MAX(*iopt, ineed);
    *zopt = MAX(*zopt, zneed);


/*             ...determine minimum and optimum flp needs for the */
/*                unnormalized cartesian (12|34) contracted batch */
/*                generation. */
    uint32_t zout2;
    erd__1111_def_blocks (0, npgto1, npgto2, npgto3, npgto4,
                          npgto12, npgto34, nxyzt, 1,
                          &zout2, NULL, NULL, NULL,
                          NULL, NULL, NULL, NULL, NULL,
                          NULL, NULL, NULL, NULL, NULL);  
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
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
