#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>

#include "erd.h"
#include "erdutil.h"

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__1111_CSGTO */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__SET_IJ_KL_PAIRS */
/*                ERD__1111_DEF_BLOCKS */
/*                ERD__PREPARE_CTR */
/*                ERD__SSSS_PCGTO_BLOCK */
/*                ERD__SSSP_PCGTO_BLOCK */
/*                ERD__SSPP_PCGTO_BLOCK */
/*                ERD__SPPP_PCGTO_BLOCK */
/*                ERD__PPPP_PCGTO_BLOCK */
/*                ERD__CTR_4INDEX_BLOCK */
/*                ERD__CTR_RS_EXPAND */
/*                ERD__CTR_TU_EXPAND */
/*                ERD__CTR_4INDEX_REORDER */
/*                ERD__MAP_IJKL_TO_IKJL */
/*  DESCRIPTION : This operation calculates a batch of contracted */
/*                electron repulsion integrals on up to four different */
/*                centers between spherical gaussian type shells. */
/*                Special fast routine for integrals involving s- and */
/*                p-type shells only! */
/*                  Input (x = 1,2,3 and 4): */
/*                    IMAX,ZMAX    =  maximum int,flp memory */
/*                    NALPHA       =  total # of exponents */
/*                    NCOEFF       =  total # of contraction coeffs */
/*                    NCSUM        =  total # of contractions */
/*                    NCGTOx       =  # of contractions for csh x */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for csh x */
/*                    SHELLx       =  the shell type for csh x */
/*                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers */
/*                                    y = 1,2,3 and 4 */
/*                    ALPHA        =  primitive exponents for csh */
/*                                    1,2,3,4 in that order */
/*                    CC           =  full set (including zeros) of */
/*                                    contraction coefficients for csh */
/*                                    1,2,3,4 in that order, for each */
/*                                    csh individually such that an */
/*                                    (I,J) element corresponds to the */
/*                                    I-th primitive and J-th contraction. */
/*                    CC(BEG)END   =  (lowest)highest nonzero primitive */
/*                                    index for contractions for csh */
/*                                    1,2,3,4 in that order. They are */
/*                                    different from (1)NPGTOx only for */
/*                                    segmented contractions */
/*                    FTABLE       =  Fm (T) table for interpolation */
/*                                    in low T region */
/*                    MGRID        =  maximum m in Fm (T) table */
/*                    NGRID        =  # of T's for which Fm (T) table */
/*                                    was set up */
/*                    TMAX         =  maximum T in Fm (T) table */
/*                    TSTEP        =  difference between two consecutive */
/*                                    T's in Fm (T) table */
/*                    TVSTEP       =  Inverse of TSTEP */
/*                    L1CACHE      =  Size of level 1 cache in units of */
/*                                    8 Byte */
/*                    TILE         =  Number of rows and columns in */
/*                                    units of 8 Byte of level 1 cache */
/*                                    square tile array used for */
/*                                    performing optimum matrix */
/*                                    transpositions */
/*                    NCTROW       =  minimum # of rows that are */
/*                                    accepted for blocked contractions */
/*                    SCREEN       =  is true, if screening will be */
/*                                    done at primitive integral level */
/*                    ICORE        =  int scratch space */
/*                    ZCORE (part) =  flp scratch space */
/*                  Output: */
/*                    NBATCH       =  # of integrals in batch */
/*                    NFIRST       =  first address location inside the */
/*                                    ZCORE array containing the first */
/*                                    integral */
/*                    ZCORE        =  full batch of contracted (12|34) */
/*                                    integrals over spherical gaussians */
/*                                    starting at ZCORE (NFIRST) */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__1111_csgto(
    uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
    uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
    bool atomic,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double alpha1[restrict static npgto1], const double alpha2[restrict static npgto2], const double alpha3[restrict static npgto3], const double alpha4[restrict static npgto4],
    const double cc1[restrict static npgto1], const double cc2[restrict static npgto2], const double cc3[restrict static npgto3], const double cc4[restrict static npgto4],
    const double norm1[restrict static npgto1], const double norm2[restrict static npgto2], const double norm3[restrict static npgto3], const double norm4[restrict static npgto4],
    uint32_t integrals_count[restrict static 1], double integrals_ptr[restrict static 81])
{

#ifdef __ERD_PROFILE__
    #ifdef _OPENMP
        const int tid = omp_get_thread_num();
    #else
        const int tid = 0;
    #endif
#endif
    ERD_PROFILE_START(erd__1111_csgto)

    const uint32_t shellp = shell1 + shell2;
    const uint32_t shellt = shellp + shell3 + shell4;
    if (atomic && ((shellt % 2) == 1)) {
        *integrals_count = 0;
        ERD_PROFILE_END(erd__1111_csgto)
        return;
    }
    /*
     * ...determine csh equality between center pairs 1,2
     *    and 3,4 in increasing order of complexity:
     *    centers -> shells -> exponents -> ctr coefficients
     * ...calculate relevant data for the [12|34] batch of
     *    integrals, such as dimensions, total # of integrals
     *    to be expected, relevant ij and kl primitive exponent
     *    pairs, etc... The integral prefactor PREFACT has been
     *    set as a parameter, its value being = 16 / sqrt(pi).
     *    Calculate here also the overall norm factor SPNORM due
     *    to presence of s- or p-type shells. The contribution
     *    to SPNORM is very simple: each s-type shell -> * 1.0,
     *    each p-type shell -> * 2.0.
     */
    const uint32_t nxyz1 = 2 * shell1 + 1;
    const uint32_t nxyz2 = 2 * shell2 + 1;
    const uint32_t nxyz3 = 2 * shell3 + 1;
    const uint32_t nxyz4 = 2 * shell4 + 1;
    const uint32_t nxyzt = nxyz1 * nxyz2 * nxyz3 * nxyz4;
    const double x12 = x1 - x2;
    const double y12 = y1 - y2;
    const double z12 = z1 - z2;
    const double rn12sq = x12 * x12 + y12 * y12 + z12 * z12;
    const double x34 = x3 - x4;
    const double y34 = y3 - y4;
    const double z34 = z3 - z4;
    const double rn34sq = x34 * x34 + y34 * y34 + z34 * z34;
    const uint32_t npgto12 = npgto1 * npgto2;
    const uint32_t npgto34 = npgto3 * npgto4;
    ERD_SIMD_ALIGN uint32_t prim1[PAD_LEN(npgto12)], prim2[PAD_LEN(npgto12)], prim3[PAD_LEN(npgto34)], prim4[PAD_LEN(npgto34)];
    const double spnorm = (double)(1 << (shell1 + shell2 + shell3 + shell4));

    ERD_PROFILE_START(erd__1111_set_ij_kl_pairs)
    ERD_SIMD_ALIGN double rhoab[PAD_LEN(npgto12)];
    ERD_SIMD_ALIGN double rhocd[PAD_LEN(npgto34)];
    uint32_t nij, nkl;
    erd__set_ij_kl_pairs(npgto1, npgto2, npgto3, npgto4,
                         x1, y1, z1, x2, y2, z2,
                         x3, y3, z3, x4, y4, z4,
                         rn12sq, rn34sq, PREFACT,
                         alpha1, alpha2, alpha3, alpha4, 
                         &nij, &nkl,
                         prim1, prim2, prim3, prim4,
                         rhoab, rhocd);
    ERD_PROFILE_END(erd__1111_set_ij_kl_pairs)
    if (nij * nkl == 0) {
        *integrals_count = 0;
        ERD_PROFILE_END(erd__1111_csgto)
        return;
    }
    /* 
     * ...decide on the primitive [12|34] block size and return array sizes and pointers for the primitive [12|34] generation.
     *    Perform also some preparation steps for contraction.
     */

    ERD_PROFILE_START(erd__1111_prepare_ctr)
    const double factor = PREFACT * spnorm;
    const uint32_t npmin = min4x32u(npgto1, npgto2, npgto3, npgto4);
    ERD_SIMD_ALIGN double norm[PAD_LEN(npmin)];
    if (npgto1 == npmin) {
        for (uint32_t i = 0; i < npgto1; i++) {
            norm[i] = factor * norm1[i];
        }
        norm1 = &norm[0];
    } else if (npgto2 == npmin) {
        for (uint32_t i = 0; i < npgto2; i++) {
            norm[i] = factor * norm2[i];
        }
        norm2 = &norm[0];
    } else if (npgto3 == npmin) {
        for (uint32_t i = 0; i < npgto3; i++) {
            norm[i] = factor * norm3[i];
        }
        norm3 = &norm[0];
    } else {
        for (uint32_t i = 0; i < npgto4; i++) {
            norm[i] = factor * norm4[i];
        }
        norm4 = &norm[0];
    }
    ERD_PROFILE_END(erd__1111_prepare_ctr)

    /* ...evaluate [12|34] in blocks over ij and kl pairs and add to final contracted (12|34) with full contracted index ranges r,s,t,u.
     *    The keyword REORDER indicates, if the primitive [12|34] blocks needs to be transposed before being contracted. */
    for (uint32_t i = 0; i < nxyzt; i++) {
        integrals_ptr[i] = 0.0;
    }
    switch (shellt) {
        case 0:
        {
            ERD_PROFILE_START(erd__ssss_pcgto_block)
            erd__ssss_pcgto_block(nij, nkl,
                x1, y1, z1,
                x2, y2, z2,
                x3, y3, z3,
                x4, y4, z4,
                alpha1, alpha2, alpha3, alpha4,
                cc1, cc2, cc3, cc4,
                prim1, prim2, prim3, prim4,
                norm1, norm2, norm3, norm4,
                rhoab, rhocd,
                integrals_ptr);
            ERD_PROFILE_END(erd__ssss_pcgto_block)
            break;
        }
        case 1:
        {
            ERD_PROFILE_START(erd__sssp_pcgto_block)
            erd__sssp_pcgto_block(nij, nkl,
                shell1, shell3, shellp, 
                x1, y1, z1,
                x2, y2, z2,
                x3, y3, z3,
                x4, y4, z4,
                alpha1, alpha2, alpha3, alpha4,
                cc1, cc2, cc3, cc4,
                prim1, prim2, prim3, prim4,
                norm1, norm2, norm3, norm4,
                rhoab, rhocd,
                integrals_ptr);
            ERD_PROFILE_END(erd__sssp_pcgto_block)
            break;
        }
        case 2:
        {
            ERD_PROFILE_START(erd__sspp_pcgto_block)
            erd__sspp_pcgto_block(nij, nkl,
                shell1, shell3, shellp, 
                x1, y1, z1,
                x2, y2, z2,
                x3, y3, z3,
                x4, y4, z4,
                alpha1, alpha2, alpha3, alpha4,
                cc1, cc2, cc3, cc4,
                prim1, prim2, prim3, prim4,
                norm1, norm2, norm3, norm4,
                rhoab, rhocd,
                integrals_ptr);
            ERD_PROFILE_END(erd__sspp_pcgto_block)
            break;
        }
        case 3:
        {
            ERD_PROFILE_START(erd__sppp_pcgto_block)
            erd__sppp_pcgto_block(nij, nkl,
                shell1, shell3, shellp, 
                x1, y1, z1,
                x2, y2, z2,
                x3, y3, z3,
                x4, y4, z4,
                alpha1, alpha2, alpha3, alpha4,
                cc1, cc2, cc3, cc4,
                prim1, prim2, prim3, prim4,
                norm1, norm2, norm3, norm4,
                rhoab, rhocd,
                integrals_ptr);
            ERD_PROFILE_END(erd__sppp_pcgto_block)
            break;
        }
        case 4:
        {
            ERD_PROFILE_START(erd__pppp_pcgto_block)
            erd__pppp_pcgto_block(nij, nkl,
                x1, y1, z1, x2, y2, z2,
                x3, y3, z3, x4, y4, z4,
                alpha1, alpha2, alpha3, alpha4,
                cc1, cc2, cc3, cc4,
                prim1, prim2, prim3, prim4,
                norm1, norm2, norm3, norm4,
                rhoab, rhocd,
                integrals_ptr);
            ERD_PROFILE_END(erd__pppp_pcgto_block)
            break;
        }
    }

/*             ...expand the contraction indices (if necessary): */
/*                   batch (nxyzt,r>=s,t>=u) --> batch (nxyzt,r,s,t,u) */
/*                and reorder the contraction indices (if necessary): */
/*                   batch (nxyzt,r,s,t,u) --> batch (nxyzt,i,j,k,l) */
/*                such that they are in final correspondence: */
/*                                    i --> 1 */
/*                                    j --> 2 */
/*                                    k --> 3 */
/*                                    l --> 4 */
/*             ...reorder contracted (12|34) batch: */
/*                      batch (nxyz1,nxyz2,nxyz3,nxyz4,rstu) --> */
/*                               batch (nxyz1,r,nxyz2,s,nxyz3,t,nxyz4,u) */
/*                Do this in three steps (if necessary): */
/*                   i) batch (nxyz1,nxyz2,nxyz3,nxyz4,rstu) --> */
/*                               batch (nxyz1,nxyz2,nxyz3,rst,nxyz4,u) */
/*                  ii) batch (nxyz1,nxyz2,nxyz3,rst,nxyz4,u) --> */
/*                               batch (nxyz1,nxyz2,rs,nxyz3,t,nxyz4,u) */
/*                 iii) batch (nxyz1,nxyz2,rs,nxyz3,t,nxyz4,u) --> */
/*                               batch (nxyz1,r,nxyz2,s,nxyz3,t,nxyz4,u) */    
/*             ...set final pointer to integrals in ZCORE array. */
    *integrals_count = nxyzt;

    ERD_PROFILE_END(erd__1111_csgto)
}
