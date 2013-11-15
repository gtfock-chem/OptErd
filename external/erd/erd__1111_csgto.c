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
int erd__1111_csgto (int imax, int zmax,
                     int nalpha, int ncoeff, int ncsum,
                     int ncgto1, int ncgto2,
                     int ncgto3, int ncgto4,
                     int npgto1, int npgto2,
                     int npgto3, int npgto4,
                     int shell1, int shell2,
                     int shell3, int shell4,
                     double x1, double y1, double z1,
                     double x2, double y2, double z2,
                     double x3, double y3, double z3,
                     double x4, double y4, double z4,
                     double *alpha, double *cc,
                     int *ccbeg, int *ccend,
                     double *ftable, int mgrid, int ngrid,
                     double tmax, double tstep, double tvstep,
                     int l1cache, int tile, int nctrow,int screen,
                     int *icore, int *nbatch, int *nfirst, double *zcore)
{
    int ftable_dim1;
    int ftable_offset;
    
    int i, j, k, l, in;
    double x12, y12, z12, x34, y34, z34;

    int zp, zq;
    int mij, nij, mkl, nkl, out, zpx, lcc1, lcc2, lcc3, lcc4, zpy,
        zpz, zqx, zqy, zqz, nctr;
    int lexp1, lexp2, lexp3, lexp4, nxyz1, nxyz2, nxyz3, nxyz4;
    int atom12, atom23;
    int mijkl;
    int atom34;
    int ixoff[4], zrho12;
    double rn12sq;
    int zrho34;
    double rn34sq;
    int empty;
    int zwork;
    int nxyzt, iprim1, iprim2, iprim3, iprim4, zscpk2, zscqk2;
    int znorm1, znorm2, znorm3, znorm4;
    int nijbeg, nklbeg, nijend, nijblk;
    int equal12;
    int nklend;
    int atomic;
    int ncgto12;
    int equal34;
    int nklblk, ncgto34, ippair, ipsave, shellp, indexr, indexs,
        indext, indexu, ctmove, ipused, ncgtor, ncgtos, mnprim, ncgtot,
        ncgtou, npgto12, npgto34;
    int npsize, ncsize, shellt, mxprim;
    double spnorm;
    int swaprs;
    int nwsize, lccseg1, xtmove;
    int swaptu;
    int lccseg2, lccseg3, lccseg4;
    int blocked;
    int zcbatch, zpbatch;
    int reorder;

    --icore;
    --zcore;
    --alpha;
    --cc;
    --ccend;
    --ccbeg;
    ftable_dim1 = mgrid - 0 + 1;
    ftable_offset = 0 + ftable_dim1 * 0;
    ftable -= ftable_offset;

    shellp = shell1 + shell2;
    shellt = shellp + shell3 + shell4;
    atom12 = x1 == x2 && y1 == y2 && z1 == z2;
    atom23 = x2 == x3 && y2 == y3 && z2 == z3;
    atom34 = x3 == x4 && y3 == y4 && z3 == z4;
    atomic = atom12 && atom34 && atom23;
    if (atomic && shellt % 2 == 1)
    {
        *nbatch = 0;
        return 0;
    }

/*             ...set the pointers to the alpha exponents, contraction */
/*                coefficients and segmented contraction boundaries. */
    lexp1 = 1;
    lexp2 = lexp1 + npgto1;
    lexp3 = lexp2 + npgto2;
    lexp4 = lexp3 + npgto3;
    lcc1 = 1;
    lcc2 = lcc1 + npgto1 * ncgto1;
    lcc3 = lcc2 + npgto2 * ncgto2;
    lcc4 = lcc3 + npgto3 * ncgto3;
    lccseg1 = 1;
    lccseg2 = lccseg1 + ncgto1;
    lccseg3 = lccseg2 + ncgto2;
    lccseg4 = lccseg3 + ncgto3;


/*             ...determine csh equality between center pairs 1,2 */
/*                and 3,4 in increasing order of complexity: */
/*                 centers -> shells -> exponents -> ctr coefficients */
    equal12 = atom12;
    if (equal12)
    {
        equal12 = ((shell1 == shell2) &&
                   (npgto1 == npgto2) &&
                   (ncgto1 == ncgto2));
        if (equal12)
        {
            k = lexp1 - 1;
            l = lexp2 - 1;
            for (i = 1; i <= npgto1; ++i)
            {
                equal12 = (equal12 &&
                           (alpha[k + i] == alpha[l + i]));
            }
            if (equal12)
            {
                k = lcc1 - 1;
                l = lcc2 - 1;
                for (j = 1; j <= ncgto1; ++j)
                {
                    if (equal12)
                    {
                        for (i = 1; i <= npgto1; ++i)
                        {
                            equal12 = (equal12 &&
                                       (cc[k + i] == cc[l + i]));
                        }
                        k += npgto1;
                        l += npgto1;
                    }
                }
            }
        }
    }
    equal34 = atom34;
    if (equal34)
    {
        equal34 = ((shell3 == shell4) &&
                   (npgto3 == npgto4) &&
                   (ncgto3 == ncgto4));
        if (equal34)
        {
            k = lexp3 - 1;
            l = lexp4 - 1;
            for (i = 1; i <= npgto3; ++i)
            {
                equal34 = (equal34 &&
                           (alpha[k + i] == alpha[l + i]));
            }
            if (equal34)
            {
                k = lcc3 - 1;
                l = lcc4 - 1;
                for (j = 1; j <= ncgto3; ++j)
                {
                    if (equal34)
                    {
                        for (i = 1; i <= npgto3; ++i)
                        {
                            equal34 = (equal34 &&
                                       (cc[k + i] == cc[l + i]));
                        }
                        k += npgto3;
                        l += npgto3;
                    }
                }
            }
        }
    }

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
    if (!atom12)
    {
        x12 = x1 - x2;
        y12 = y1 - y2;
        z12 = z1 - z2;
        rn12sq = x12 * x12 + y12 * y12 + z12 * z12;
    }
    else
    {
        x12 = 0.0;
        y12 = 0.0;
        z12 = 0.0;
        rn12sq = 0.0;
    }
    if (!atom34)
    {
        x34 = x3 - x4;
        y34 = y3 - y4;
        z34 = z3 - z4;
        rn34sq = x34 * x34 + y34 * y34 + z34 * z34;
    }
    else
    {
        x34 = 0.0;
        y34 = 0.0;
        z34 = 0.0;
        rn34sq = 0.0;
    }
    if (equal12)
    {
        npgto12 = npgto1 * (npgto1 + 1) / 2;
        ncgto12 = ncgto1 * (ncgto1 + 1) / 2;
    }
    else
    {
        npgto12 = npgto1 * npgto2;
        ncgto12 = ncgto1 * ncgto2;
    }
    if (equal34)
    {
        npgto34 = npgto3 * (npgto3 + 1) / 2;
        ncgto34 = ncgto3 * (ncgto3 + 1) / 2;
    }
    else
    {
        npgto34 = npgto3 * npgto4;
        ncgto34 = ncgto3 * ncgto4;
    }
    nctr = ncgto12 * ncgto34;
    iprim1 = 1;
    iprim2 = iprim1 + npgto12;
    iprim3 = iprim2 + npgto12;
    iprim4 = iprim3 + npgto34;
    swaprs = (npgto1 > npgto2);
    swaptu = (npgto3 > npgto4);
    spnorm = 1.;
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
    #if 0
    erd__set_ij_kl_pairs  (npgto1, npgto2, npgto3, npgto4,
                           npgto12, npgto34,
                           atom12, atom34, equal12, equal34,
                           swaprs, swaptu,
                           x1, y1, z1, x2, y2, z2,
                           x3, y3, z3, x4, y4, z4,
                           rn12sq, rn34sq, PREFACT,
                           &alpha[lexp1], &alpha[lexp2],
                           &alpha[lexp3], &alpha[lexp4],
                           &ftable[ftable_offset], mgrid, ngrid,
                           tmax, tstep, tvstep, screen,
                           &empty, &nij, &nkl,
                           &icore[iprim1], &icore[iprim2],
                           &icore[iprim3], &icore[iprim4], &zcore[1]);
    #else
    double aa = PREFACT;
    erd__set_ij_kl_pairs_  (&npgto1, &npgto2, &npgto3, &npgto4,
                           &npgto12, &npgto34,
                           &atom12, &atom34, &equal12, &equal34,
                           &swaprs, &swaptu,
                           &x1, &y1, &z1, &x2, &y2, &z2,
                           &x3, &y3, &z3, &x4, &y4, &z4,
                           &rn12sq, &rn34sq, &aa,
                           &alpha[lexp1], &alpha[lexp2],
                           &alpha[lexp3], &alpha[lexp4],
                           &ftable[ftable_offset], &mgrid, &ngrid,
                           &tmax, &tstep, &tvstep, &screen,
                           &empty, &nij, &nkl,
                           &icore[iprim1], &icore[iprim2],
                           &icore[iprim3], &icore[iprim4], &zcore[1]);
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
                          nij, nkl, ncgto12, ncgto34, nctr, nxyzt,
                          l1cache, nctrow, 0,
                          &nijblk, &nklblk, &npsize, &ncsize, &nwsize,
                          &mxprim, &mnprim,
                          &zcbatch, &zpbatch, &zwork,
                          &znorm1, &znorm2, &znorm3, &znorm4,
                          &zrho12, &zrho34,
                          &zp, &zpx, &zpy, &zpz, &zscpk2,
                          &zq, &zqx, &zqy, &zqz, &zscqk2);
    blocked = 0;

    erd__prepare_ctr (ncsize, nij, nkl,
                      npgto1, npgto2, npgto3, npgto4,
                      shell1, shell2, shell3, shell4,
                      &alpha[lexp1], &alpha[lexp2],
                      &alpha[lexp3], &alpha[lexp4], PREFACT,
                      spnorm, equal12, equal34, blocked,
                      &zcore[1], &zcore[znorm1], &zcore[znorm2],
                      &zcore[znorm3], &zcore[znorm4],
                      &zcore[zrho12], &zcore[zrho34], &zcore[zcbatch]);
    ipused = iprim4 + npgto34;
    ipsave = ipused + mnprim;
    ippair = ipsave + mxprim;

/*             ...evaluate [12|34] in blocks over ij and kl pairs */
/*                and add to final contracted (12|34) with full */
/*                contracted index ranges r,s,t,u. The keyword REORDER */
/*                indicates, if the primitive [12|34] blocks needs to */
/*                be transposed before being contracted. */
    reorder = 0;
    nijbeg = 1;
    nijend = nij;
    mij = nij;
    nklbeg = 1;
    nklend = nkl;
    mkl = nkl;
    mijkl = mij * mkl;
    npsize = nxyzt * mijkl;

    if (shellt == 0)
    {
        erd__ssss_pcgto_block (npsize, atomic, atom12, atom34,
                               mij, mkl, nij, nijbeg, nijend,
                               nkl, nklbeg, nklend,
                               npgto1, npgto2, npgto3, npgto4,
                               x1, y1, z1, x2, y2, z2,
                               x3, y3, z3, x4, y4, z4,
                               x12, y12, z12, x34, y34, z34,
                               &alpha[lexp1], &alpha[lexp2],
                               &alpha[lexp3], &alpha[lexp4],
                               &ftable[ftable_offset], mgrid, ngrid,
                               tmax, tstep, tvstep,
                               &icore[iprim1 + nijbeg - 1],
                               &icore[iprim2 + nijbeg - 1],
                               &icore[iprim3 + nklbeg - 1],
                               &icore[iprim4 + nklbeg - 1],
                               &zcore[znorm1], &zcore[znorm2],
                               &zcore[znorm3], &zcore[znorm4],
                               &zcore[zrho12], &zcore[zrho34],
                               &zcore[zp], &zcore[zpx],
                               &zcore[zpy], &zcore[zpz], &zcore[zscpk2],
                               &zcore[zq], &zcore[zqx],
                               &zcore[zqy], &zcore[zqz], &zcore[zscqk2],
                               &zcore[zpbatch]);
    }
    else if (shellt == 1)
    {
        erd__sssp_pcgto_block (npsize, atomic, atom12, atom34,
                               mij, mkl, nij, nijbeg, nijend,
                               nkl, nklbeg, nklend,
                               npgto1, npgto2, npgto3, npgto4,
                               shell1, shell3, shellp, 
                               x1, y1, z1, x2, y2, z2,
                               x3, y3, z3, x4, y4, z4,
                               x12, y12, z12, x34, y34, z34,
                               &alpha[lexp1], &alpha[lexp2],
                               &alpha[lexp3], &alpha[lexp4],
                               &ftable[ftable_offset], mgrid, ngrid,
                               tmax, tstep, tvstep,
                               &icore[iprim1 + nijbeg - 1],
                               &icore[iprim2 + nijbeg - 1],
                               &icore[iprim3 + nklbeg - 1],
                               &icore[iprim4 + nklbeg - 1],
                               &zcore[znorm1], &zcore[znorm2],
                               &zcore[znorm3], &zcore[znorm4],
                               &zcore[zrho12], &zcore[zrho34],
                               &zcore[zp], &zcore[zpx],
                               &zcore[zpy], &zcore[zpz], &zcore[zscpk2],
                               &zcore[zq], &zcore[zqx],
                               &zcore[zqy], &zcore[zqz], &zcore[zscqk2],
                               &zcore[zpbatch]);
    }
    else if (shellt == 2)
    {
        erd__sspp_pcgto_block (npsize, atomic, atom12, atom34,
                               mij, mkl, nij, nijbeg, nijend,
                               nkl, nklbeg, nklend,
                               npgto1, npgto2, npgto3, npgto4,
                               shell1, shell3, shellp, 
                               x1, y1, z1, x2, y2, z2,
                               x3, y3, z3, x4, y4, z4,
                               x12, y12, z12, x34, y34, z34,
                               &alpha[lexp1], &alpha[lexp2],
                               &alpha[lexp3], &alpha[lexp4],
                               &ftable[ftable_offset], mgrid, ngrid,
                               tmax, tstep, tvstep,
                               &icore[iprim1 + nijbeg - 1],
                               &icore[iprim2 + nijbeg - 1],
                               &icore[iprim3 + nklbeg - 1],
                               &icore[iprim4 + nklbeg - 1],
                               &zcore[znorm1], &zcore[znorm2],
                               &zcore[znorm3], &zcore[znorm4],
                               &zcore[zrho12], &zcore[zrho34],
                               &zcore[zp], &zcore[zpx],
                               &zcore[zpy], &zcore[zpz], &zcore[zscpk2],
                               &zcore[zq], &zcore[zqx],
                               &zcore[zqy], &zcore[zqz], &zcore[zscqk2],
                               &zcore[zpbatch]);
    }
    else if (shellt == 3)
    {
        erd__sppp_pcgto_block (npsize, atomic, atom12, atom34,
                               mij, mkl, nij, nijbeg, nijend,
                               nkl, nklbeg, nklend,
                               npgto1, npgto2, npgto3, npgto4,
                               shell1, shell3, shellp, 
                               x1, y1, z1, x2, y2, z2,
                               x3, y3, z3, x4, y4, z4,
                               x12, y12, z12, x34, y34, z34,
                               &alpha[lexp1], &alpha[lexp2],
                               &alpha[lexp3], &alpha[lexp4],
                               &ftable[ftable_offset], mgrid, ngrid,
                               tmax, tstep, tvstep,
                               &icore[iprim1 + nijbeg - 1],
                               &icore[iprim2 + nijbeg - 1],
                               &icore[iprim3 + nklbeg - 1],
                               &icore[iprim4 + nklbeg - 1],
                               &zcore[znorm1], &zcore[znorm2],
                               &zcore[znorm3], &zcore[znorm4],
                               &zcore[zrho12], &zcore[zrho34],
                               &zcore[zp], &zcore[zpx],
                               &zcore[zpy], &zcore[zpz], &zcore[zscpk2],
                               &zcore[zq], &zcore[zqx],
                               &zcore[zqy], &zcore[zqz], &zcore[zscqk2],
                               &zcore[zpbatch]);
    }
    else
    {
        erd__pppp_pcgto_block (npsize, atomic, atom12, atom34,
                               mij, mkl, nij, nijbeg, nijend,
                               nkl, nklbeg, nklend,
                               npgto1, npgto2, npgto3, npgto4,
                               x1, y1, z1, x2, y2, z2,
                               x3, y3, z3, x4, y4, z4,
                               x12, y12, z12, x34, y34, z34,
                               &alpha[lexp1], &alpha[lexp2],
                               &alpha[lexp3], &alpha[lexp4],
                               &ftable[ftable_offset], mgrid, ngrid,
                               tmax, tstep, tvstep,
                               &icore[iprim1 + nijbeg - 1],
                               &icore[iprim2 + nijbeg - 1],
                               &icore[iprim3 + nklbeg - 1],
                               &icore[iprim4 + nklbeg - 1],
                               &zcore[znorm1], &zcore[znorm2],
                               &zcore[znorm3], &zcore[znorm4],
                               &zcore[zrho12], &zcore[zrho34],
                               &zcore[zp], &zcore[zpx],
                               &zcore[zpy], &zcore[zpz], &zcore[zscpk2],
                               &zcore[zq], &zcore[zqx],
                               &zcore[zqy], &zcore[zqz], &zcore[zscqk2],
                               &zcore[zpbatch]);
    }

    erd__ctr_4index_block (npsize, ncsize, nwsize, nxyzt, mijkl,
                           mij, mkl,
                           ncgto12, ncgto34,
                           npgto1, npgto2, npgto3, npgto4,
                           ncgto1, ncgto2, ncgto3, ncgto4,
                           mxprim, mnprim,
                           &cc[lcc1], &cc[lcc2], &cc[lcc3], &cc[lcc4],
                           &ccbeg[lccseg1], &ccbeg[lccseg2],
                           &ccbeg[lccseg3], &ccbeg[lccseg4],
                           &ccend[lccseg1], &ccend[lccseg2],
                           &ccend[lccseg3], &ccend[lccseg4],
                           &icore[iprim1 + nijbeg - 1],
                           &icore[iprim2 + nijbeg - 1],
                           &icore[iprim3 + nklbeg - 1],
                           &icore[iprim4 + nklbeg - 1],
                           l1cache, tile, nctrow,
                           equal12, equal34, swaprs, swaptu,
                           reorder, blocked,
                           &icore[ipused], &icore[ipsave], &icore[ippair],
                           &zcore[zpbatch], &zcore[zwork], &zcore[zcbatch]);

/*             ...expand the contraction indices (if necessary): */
/*                   batch (nxyzt,r>=s,t>=u) --> batch (nxyzt,r,s,t,u) */
/*                and reorder the contraction indices (if necessary): */
/*                   batch (nxyzt,r,s,t,u) --> batch (nxyzt,i,j,k,l) */
/*                such that they are in final correspondence: */
/*                                    i --> 1 */
/*                                    j --> 2 */
/*                                    k --> 3 */
/*                                    l --> 4 */
    nctr = ncgto1 * ncgto2 * ncgto3 * ncgto4;
    *nbatch = nxyzt * nctr;
    in = zcbatch;
    out = zcbatch + *nbatch;
    #if 0
    if (equal12 && ncgto12 > 1)
    {
        erd__ctr_rs_expand_ (&nxyzt, &ncgto12, &ncgto34, ncgto1, ncgto2,
                              &zcore[in], &zcore[out]);
        i = in;
        in = out;
        out = i;
    }
    ncgto12 = ncgto1 * ncgto2;
    if (equal34 && ncgto34 > 1)
    {
        int i1;
        i1 = nxyzt * ncgto12;
        erd__ctr_tu_expand_ (&i1, &ncgto34, ncgto3, ncgto4, &zcore[in],
                              &zcore[out]);
        i = in;
        in = out;
        out = i;
    }
    ncgto34 = ncgto3 * ncgto4;
    if ((swaprs || swaptu) && nctr > 1)
    {
        ixoff[0] = 1;
        ixoff[1] = ncgto1;
        ixoff[2] = ncgto1 * ncgto2;
        ixoff[3] = ncgto1 * ncgto2 * ncgto3;
        if (swaprs)
        {
            indexr = 2;
            indexs = 1;
            ncgtor = ncgto2;
            ncgtos = ncgto1;
        }
        else
        {
            indexr = 1;
            indexs = 2;
            ncgtor = ncgto1;
            ncgtos = ncgto2;
        }
        if (swaptu)
        {
            indext = 4;
            indexu = 3;
            ncgtot = ncgto4;
            ncgtou = ncgto3;
        }
        else
        {
            indext = 3;
            indexu = 4;
            ncgtot = ncgto3;
            ncgtou = ncgto4;
        }
        erd__ctr_4index_reorder_ (&nxyzt, &nctr, &ncgtor, &ncgtos, &ncgtot,
                                   &ncgtou, &ixoff[indexr - 1],
                                   &ixoff[indexs - 1], &ixoff[indext - 1],
                                   &ixoff[indexu - 1], &zcore[in],
                                   &zcore[out]);
        i = in;
        in = out;
        out = i;
    }
    #endif

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
    xtmove = nxyz2 * nxyz3 * nxyz4;
    ctmove = ncgto12 * ncgto3;
    if (xtmove > 1 || ctmove > 1)
    {
        if (nxyz4 > 1)
        {
            erd__map_ijkl_to_ikjl (nxyz1 * nxyz2 * nxyz3,
                                   nxyz4, ctmove, ncgto4,
                                   &zcore[in], &zcore[out]);
            i = in;
            in = out;
            out = i;
        }
        if (nxyz3 > 1 && ncgto12 > 1)
        {
            erd__map_ijkl_to_ikjl (nxyz1 * nxyz2, nxyz3,
                                   ncgto12, nxyz4 * ncgto34,
                                   &zcore[in], &zcore[out]);
            i = in;
            in = out;
            out = i;
        }
        if (nxyz2 > 1 && ncgto1 > 1)
        {
            erd__map_ijkl_to_ikjl (nxyz1, nxyz2, ncgto1,
                                   nxyz3 * nxyz4 * ncgto2 * ncgto34,
                                   &zcore[in], &zcore[out]);
            i = in;
            in = out;
            out = i;
        }
    }
/*             ...set final pointer to integrals in ZCORE array. */
    *nfirst = in;

    return 0;
}


int erd__1111_csgto_ (int *imax, int *zmax,
                      int *nalpha, int * ncoeff, int * ncsum,
                      int * ncgto1, int * ncgto2, int * ncgto3,
                      int * ncgto4, int * npgto1, int * npgto2,
                      int * npgto3, int * npgto4, int * shell1,
                      int * shell2, int * shell3, int * shell4,
                      double * x1, double * y1, double * z1,
                      double * x2, double * y2, double * z2,
                      double * x3, double * y3, double * z3,
                      double * x4, double * y4, double * z4,
                      double * alpha, double * cc, int * ccbeg,
                      int * ccend, double * ftable, int * mgrid,
                      int * ngrid, double * tmax, double * tstep,
                      double * tvstep, int * l1cache, int * tile,
                      int * nctrow, int * screen, int * icore,
                      int * nbatch, int * nfirst, double * zcore)
{
    erd__1111_csgto (*imax, *zmax,
                     *nalpha, *ncoeff, *ncsum,
                     *ncgto1, *ncgto2,
                     *ncgto3, *ncgto4,
                     *npgto1, *npgto2,
                     *npgto3, *npgto4,
                     *shell1, *shell2,
                     *shell3, *shell4,
                     *x1, *y1, *z1,
                     *x2, *y2, *z2,
                     *x3, *y3, *z3,
                     *x4, *y4, *z4,
                     alpha, cc, ccbeg, ccend,
                     ftable, *mgrid, *ngrid,
                     *tmax, *tstep, *tvstep,
                     *l1cache, *tile, *nctrow, *screen,
                     icore, nbatch, nfirst, zcore);
    
    return 0;
}