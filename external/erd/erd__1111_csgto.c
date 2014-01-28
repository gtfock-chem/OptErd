#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


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
int erd__1111_csgto (int zmax, int npgto1, int npgto2,
                     int npgto3, int npgto4,
                     int shell1, int shell2,
                     int shell3, int shell4,
                     double x1, double y1, double z1,
                     double x2, double y2, double z2,
                     double x3, double y3, double z3,
                     double x4, double y4, double z4,
                     double *alpha, double *cc, int screen,
                     int *icore, int *nbatch, int *nfirst, double *zcore)
{    
    int i;
    double x12, y12, z12, x34, y34, z34;
    int zp, zq;
    int nij, nkl, zpx, lcc1, lcc2, lcc3, lcc4, zpy,
        zpz, zqx, zqy, zqz;
    int lexp1, lexp2, lexp3, lexp4;
    int nxyz1, nxyz2, nxyz3, nxyz4;
    int atom12, atom23;
    int atom34;
    int zrho12;
    double rn12sq;
    int zrho34;
    double rn34sq;
    int empty;
    int nxyzt, iprim1, iprim2, iprim3, iprim4, zscpk2, zscqk2;
    int znorm1, znorm2, znorm3, znorm4;
    int atomic;
    int shellp, npgto12, npgto34;
    int shellt;
    double spnorm;
    int zcbatch;
#ifdef __ERD_PROFILE__    
    uint64_t start_clock, end_clock;
    int tid = omp_get_thread_num();
#endif

    shellp = shell1 + shell2;
    shellt = shellp + shell3 + shell4;
    atom12 = ((x1 == x2) && (y1 == y2) && (z1 == z2));
    atom23 = ((x2 == x3) && (y2 == y3) && (z2 == z3));
    atom34 = ((x3 == x4) && (y3 == y4) && (z3 == z4));
    atomic = atom12 && atom34 && atom23;
    if (atomic && shellt % 2 == 1)
    {
        *nbatch = 0;
        return 0;
    }
/*             ...set the pointers to the alpha exponents, contraction */
/*                coefficients and segmented contraction boundaries. */
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
/*             ...calculate relevant data for the [12|34] batch of */
/*                integrals, such as dimensions, total # of integrals */
/*                to be expected, relevant ij and kl primitive exponent */
/*                pairs, etc... The integral prefactor PREFACT has been */
/*                set as a parameter, its value being = 16 / sqrt(pi). */
/*                Calculate here also the overall norm factor SPNORM due */
/*                to presence of s- or p-type shells. The contribution */
/*                to SPNORM is very simple: each s-type shell -> * 1.0, */
/*                each p-type shell -> * 2.0. */
    nxyz1 = shell1 + shell1 + 1;
    nxyz2 = shell2 + shell2 + 1;
    nxyz3 = shell3 + shell3 + 1;
    nxyz4 = shell4 + shell4 + 1;
    nxyzt = nxyz1 * nxyz2 * nxyz3 * nxyz4;
    x12 = x1 - x2;
    y12 = y1 - y2;
    z12 = z1 - z2;
    rn12sq = x12 * x12 + y12 * y12 + z12 * z12;
    x34 = x3 - x4;
    y34 = y3 - y4;
    z34 = z3 - z4;
    rn34sq = x34 * x34 + y34 * y34 + z34 * z34;
    npgto12 = npgto1 * npgto2;
    npgto34 = npgto3 * npgto4;
    iprim1 = 0;
    iprim2 = iprim1 + npgto12;
    iprim3 = iprim2 + npgto12;
    iprim4 = iprim3 + npgto34;
    spnorm = 1.0;
    if (shell1 == 1)
    {
        spnorm += spnorm;
    }
    if (shell2 == 1)
    {
        spnorm += spnorm;
    }
    if (shell3 == 1)
    {
        spnorm += spnorm;
    }
    if (shell4 == 1)
    {
        spnorm += spnorm;
    }

#ifdef __ERD_PROFILE__
    start_clock = __rdtsc();
#endif
    erd__set_ij_kl_pairs  (npgto1, npgto2, npgto3, npgto4,
                           x1, y1, z1, x2, y2, z2,
                           x3, y3, z3, x4, y4, z4,
                           rn12sq, rn34sq, PREFACT,
                           &alpha[lexp1], &alpha[lexp2],
                           &alpha[lexp3], &alpha[lexp4],
                           screen, &empty, &nij, &nkl,
                           &icore[iprim1], &icore[iprim2],
                           &icore[iprim3], &icore[iprim4], &zcore[0]);
#ifdef __ERD_PROFILE__
    end_clock = __rdtsc();
    erd_ticks[tid][erd__set_ij_kl_pairs_ticks_1111] += (end_clock - start_clock);
#endif
    if (empty)
    {
        *nbatch = 0;
        return 0;
    }
/*             ...decide on the primitive [12|34] block size and */
/*                return array sizes and pointers for the primitive */
/*                [12|34] generation. Perform also some preparation */
/*                steps for contraction. */
    erd__1111_def_blocks (zmax, npgto1, npgto2, npgto3, npgto4,
                          nij, nkl, nxyzt, 0,
                          &zcbatch, &znorm1, &znorm2,
                          &znorm3, &znorm4,
                          &zrho12, &zrho34,
                          &zp, &zpx, &zpy, &zpz, &zscpk2,
                          &zq, &zqx, &zqy, &zqz, &zscqk2);
#ifdef __ERD_PROFILE__
    start_clock = __rdtsc();
#endif
    erd__prepare_ctr (npgto1, npgto2, npgto3, npgto4,
                      shell1, shell2, shell3, shell4,
                      &alpha[lexp1], &alpha[lexp2],
                      &alpha[lexp3], &alpha[lexp4],
                      spnorm, &zcore[znorm1], &zcore[znorm2],
                      &zcore[znorm3], &zcore[znorm4]);
#ifdef __ERD_PROFILE__
    end_clock = __rdtsc();
    erd_ticks[tid][erd__prepare_ctr_ticks_1111] += (end_clock - start_clock);
#endif

/*             ...evaluate [12|34] in blocks over ij and kl pairs */
/*                and add to final contracted (12|34) with full */
/*                contracted index ranges r,s,t,u. The keyword REORDER */
/*                indicates, if the primitive [12|34] blocks needs to */
/*                be transposed before being contracted. */
    for (i = 0; i < nxyzt; i++)
    {
        zcore[zcbatch + i] = 0.0;   
    }
    if (shellt == 0)
    {
    #ifdef __ERD_PROFILE__
        start_clock = __rdtsc();
    #endif
        erd__ssss_pcgto_block (nij, nkl,
                               x1, y1, z1, x2, y2, z2,
                               x3, y3, z3, x4, y4, z4,
                               &alpha[lexp1], &alpha[lexp2],
                               &alpha[lexp3], &alpha[lexp4],
                               &cc[lcc1], &cc[lcc2],
                               &cc[lcc3], &cc[lcc4],
                               &icore[iprim1],
                               &icore[iprim2],
                               &icore[iprim3],
                               &icore[iprim4],
                               &zcore[znorm1], &zcore[znorm2],
                               &zcore[znorm3], &zcore[znorm4],
                               &zcore[zrho12], &zcore[zrho34],
                               &zcore[zp], &zcore[zpx],
                               &zcore[zpy], &zcore[zpz], &zcore[zscpk2],
                               &zcore[zq], &zcore[zqx],
                               &zcore[zqy], &zcore[zqz], &zcore[zscqk2],
                               &zcore[zcbatch]);
    #ifdef __ERD_PROFILE__
        end_clock = __rdtsc();
        erd_ticks[tid][erd__ssss_pcgto_block_ticks] += (end_clock - start_clock);
    #endif
    }
    else if (shellt == 1)
    {
    #ifdef __ERD_PROFILE__
        start_clock = __rdtsc();
    #endif
        erd__sssp_pcgto_block (nij, nkl,
                               shell1, shell3, shellp, 
                               x1, y1, z1, x2, y2, z2,
                               x3, y3, z3, x4, y4, z4,
                               &alpha[lexp1], &alpha[lexp2],
                               &alpha[lexp3], &alpha[lexp4],
                               &cc[lcc1], &cc[lcc2],
                               &cc[lcc3], &cc[lcc4],
                               &icore[iprim1],
                               &icore[iprim2],
                               &icore[iprim3],
                               &icore[iprim4],
                               &zcore[znorm1], &zcore[znorm2],
                               &zcore[znorm3], &zcore[znorm4],
                               &zcore[zrho12], &zcore[zrho34],
                               &zcore[zp], &zcore[zpx],
                               &zcore[zpy], &zcore[zpz], &zcore[zscpk2],
                               &zcore[zq], &zcore[zqx],
                               &zcore[zqy], &zcore[zqz], &zcore[zscqk2],
                               &zcore[zcbatch]);
    #ifdef __ERD_PROFILE__
        end_clock = __rdtsc();
        erd_ticks[tid][erd__sssp_pcgto_block_ticks] += (end_clock - start_clock);
    #endif
    }
    else if (shellt == 2)
    {
    #ifdef __ERD_PROFILE__
        start_clock = __rdtsc();
    #endif
        erd__sspp_pcgto_block (nij, nkl,
                               shell1, shell3, shellp, 
                               x1, y1, z1, x2, y2, z2,
                               x3, y3, z3, x4, y4, z4,
                               &alpha[lexp1], &alpha[lexp2],
                               &alpha[lexp3], &alpha[lexp4],
                               &cc[lcc1], &cc[lcc2],
                               &cc[lcc3], &cc[lcc4],
                               &icore[iprim1],
                               &icore[iprim2],
                               &icore[iprim3],
                               &icore[iprim4],
                               &zcore[znorm1], &zcore[znorm2],
                               &zcore[znorm3], &zcore[znorm4],
                               &zcore[zrho12], &zcore[zrho34],
                               &zcore[zp], &zcore[zpx],
                               &zcore[zpy], &zcore[zpz], &zcore[zscpk2],
                               &zcore[zq], &zcore[zqx],
                               &zcore[zqy], &zcore[zqz], &zcore[zscqk2],
                               &zcore[zcbatch]);
    #ifdef __ERD_PROFILE__
        end_clock = __rdtsc();
        erd_ticks[tid][erd__sspp_pcgto_block_ticks] += (end_clock - start_clock);
    #endif
    }
    else if (shellt == 3)
    {
    #ifdef __ERD_PROFILE__
        start_clock = __rdtsc();
    #endif
        erd__sppp_pcgto_block (nij, nkl,
                               shell1, shell3, shellp, 
                               x1, y1, z1, x2, y2, z2,
                               x3, y3, z3, x4, y4, z4,
                               &alpha[lexp1], &alpha[lexp2],
                               &alpha[lexp3], &alpha[lexp4],
                               &cc[lcc1], &cc[lcc2],
                               &cc[lcc3], &cc[lcc4],
                               &icore[iprim1],
                               &icore[iprim2],
                               &icore[iprim3],
                               &icore[iprim4],
                               &zcore[znorm1], &zcore[znorm2],
                               &zcore[znorm3], &zcore[znorm4],
                               &zcore[zrho12], &zcore[zrho34],
                               &zcore[zp], &zcore[zpx],
                               &zcore[zpy], &zcore[zpz], &zcore[zscpk2],
                               &zcore[zq], &zcore[zqx],
                               &zcore[zqy], &zcore[zqz], &zcore[zscqk2],
                               &zcore[zcbatch]);
    #ifdef __ERD_PROFILE__
        end_clock = __rdtsc();
        erd_ticks[tid][erd__sppp_pcgto_block_ticks] += (end_clock - start_clock);
    #endif
    }
    else
    {
    #ifdef __ERD_PROFILE__
        start_clock = __rdtsc();
    #endif
        erd__pppp_pcgto_block (nij, nkl,
                               x1, y1, z1, x2, y2, z2,
                               x3, y3, z3, x4, y4, z4,
                               &alpha[lexp1], &alpha[lexp2],
                               &alpha[lexp3], &alpha[lexp4],
                               &cc[lcc1], &cc[lcc2],
                               &cc[lcc3], &cc[lcc4],
                               &icore[iprim1],
                               &icore[iprim2],
                               &icore[iprim3],
                               &icore[iprim4],
                               &zcore[znorm1], &zcore[znorm2],
                               &zcore[znorm3], &zcore[znorm4],
                               &zcore[zrho12], &zcore[zrho34],
                               &zcore[zp], &zcore[zpx],
                               &zcore[zpy], &zcore[zpz], &zcore[zscpk2],
                               &zcore[zq], &zcore[zqx],
                               &zcore[zqy], &zcore[zqz], &zcore[zscqk2],
                               &zcore[zcbatch]);
    #ifdef __ERD_PROFILE__
        end_clock = __rdtsc();
        erd_ticks[tid][erd__pppp_pcgto_block_ticks] += (end_clock - start_clock);
    #endif
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
    *nbatch = nxyzt;
    *nfirst = zcbatch + 1;

    return 0;
}
