#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__CSGTO */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__SET_ABCD */
/*                ERD__SET_IJ_KL_PAIRS */
/*                ERD__E0F0_DEF_BLOCKS */
/*                ERD__PREPARE_CTR */
/*                ERD__E0F0_PCGTO_BLOCK */
/*                ERD__CTR_4INDEX_BLOCK */
/*                ERD__CTR_RS_EXPAND */
/*                ERD__CTR_TU_EXPAND */
/*                ERD__CTR_4INDEX_REORDER */
/*                ERD__TRANSPOSE_BATCH */
/*                ERD__XYZ_TO_RY_ABCD */
/*                ERD__CARTESIAN_NORMS */
/*                ERD__HRR_MATRIX */
/*                ERD__HRR_TRANSFORM */
/*                ERD__SPHERICAL_TRANSFORM */
/*                ERD__NORMALIZE_CARTESIAN */
/*                ERD__MOVE_RY */
/*  DESCRIPTION : This operation calculates a batch of contracted */
/*                electron repulsion integrals on up to four different */
/*                centers between spherical or cartesian gaussian type */
/*                shells. */
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
/*                                    I-th primitive and J-th contraction */
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
/*                    SPHERIC      =  is true, if spherical integrals */
/*                                    are wanted, false if cartesian */
/*                                    ones are wanted */
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
/*                                    integrals over cartesian or */
/*                                    spherical gaussians starting at */
/*                                    ZCORE (NFIRST) */
/*              --- NOTES ABOUT THE OVERALL INTEGRAL PREFACTOR --- */
/*                The overal integral prefactor is defined here as */
/*                follows. Consider the normalization factors for a */
/*                primitive cartesian GTO and for a spherical GTO */
/*                belonging to the angular momentum L = l+m+n: */
/*                    lmn                        l m n           2 */
/*                 GTO  (x,y,z) = N (l,m,n,a) * x y z  * exp (-ar ) */
/*                    a */
/*                    LM                     L    LM                2 */
/*                 GTO  (r,t,p) = N (L,a) * r  * Y  (t,p) * exp (-ar ) */
/*                    a */
/*                where a = alpha exponent, t = theta and p = phi and */
/*                N (l,m,n,a) and N (L,a) denote the respective */
/*                cartesian and spherical normalization factors such */
/*                that: */
/*                              lmn            lmn */
/*                 integral {GTO  (x,y,z) * GTO   (x,y,z) dx dy dz} = 1 */
/*                              a              a */
/*                              LM             LM */
/*                 integral {GTO  (r,t,p) * GTO  (r,t,p) dr dt dp} = 1 */
/*                              a              a */
/*                The normalization constants have then the following */
/*                values, assuming the spherical harmonics are */
/*                normalized: */
/*                              _____________________________________ */
/*                             /      2^(2L+1+1/2) * a^((2L+3)/2) */
/*            N (l,m,n,a) =   / ---------------------------------------- */
/*                          \/ (2l-1)!!(2m-1)!!(2n-1)!! * pi * sqrt (pi) */
/*                                   ____________________________ */
/*                                  / 2^(2L+3+1/2) * a^((2L+3)/2) */
/*                     N (L,a) =   / ----------------------------- */
/*                               \/     (2L+1)!! * sqrt (pi) */
/*                Note, that the extra pi under the square root in */
/*                N (l,m,n,a) belongs to the normalization of the */
/*                spherical harmonic functions and therefore does not */
/*                appear in the expression for N (L,a). The common */
/*                L-,l-,m-,n- and a-independent part of the cartesian */
/*                norm is a scalar quantity needed for all integrals */
/*                no matter what L-,l-,m-,n- and a-values they have: */
/*                                       _____________ */
/*                                      /  2^(1+1/2) */
/*                     N (0,0,0,0) =   / -------------- */
/*                                   \/  pi * sqrt (pi) */
/*                Also every ERI integral has a factor of 2*pi**(5/2) */
/*                associated with it, hence the overall common factor */
/*                for all integrals will be N(0,0,0,0)**4 times */
/*                2*pi**(5/2), which is equal to: */
/*                            PREFACT = 16 / sqrt(pi) */
/*                and is set as a parameter inside the present routine. */
/*                The alpha exponent dependent part of the norms: */
/*                                   ______________ */
/*                                 \/ a^((2L+3)/2) */
/*                will be calculated separately (see below) and their */
/*                inclusion in evaluating the primitive cartesian */
/*                [E0|F0] integrals will be essential for numerical */
/*                stability during contraction. */
/* ------------------------------------------------------------------------ */
int erd__csgto (int imax, int zmax,
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
                int l1cache, int tile, int nctrow,
                int spheric, int screen, int *icore,
                int *nbatch, int *nfirst, double *zcore)
{
    int ftable_dim1, ftable_offset;

    int nxyzhrr;   
    int in;
    double xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd;
    int zp, zq;
    double abx;
    int zb00, zb01, zb10;
    double aby;
    int mij, nij;
    double abz, cdx;
    int mkl, nkl;
    double cdy, cdz;
    int out, zpx, lcc1, lcc2, lcc3, lcc4, zpy, zpz, zqx, zqy, zqz,
        pos1, pos2, lcca, lccb, lccc, lccd;
    int tr1234;
    int zc00x, ngqp, move, nctr, nmom, nrya, nryb, nryc, nryd,
        temp, zc00y, zc00z, zpax, zpay, zpaz, zqcx, zqcy, zqcz, zd00x,
        zd00y, zd00z;
    int zrts, zwts;     
    int lexp1, lexp2, lexp3, lexp4;
    int zbase;
    int mijkl, ihscr, iused, lexpa, lexpb, lexpc;
    int swap12;
    int lexpd, ixoff[4];
    int swap34;
    int nrota, nrotb, nrotc, nrotd, ihrow, nrowa, zused, nrowb,
        nrowc, nrowd;
    int empty;
    int ztval, zhrot, nxyza, nxyzb, nxyzc, nxyzd, zwork, nint2d,
        nxyzp, nxyzq, nxyzt, zscpk2, zscqk2, nijbeg, nklbeg;
    int atomab, atomcd;
    int indexa, indexb, indexc, indexd;
    int atomic;
    int ncgtoa, ncgtob, ncgtoc, iprima, iprimb, iprimc, iprimd,
        ippair, ncgtod, ipsave, nijblk, indexr, indexs, indext, indexu,
        ipused, ncgtor, ncgtos, ncgtot, mnprim, ncgtou, isrowa, ihnrow,
        isrowb;
    int memory;
    int isrowc, isrowd, ngqscr, mxprim, nklblk;
    int swaprs;
    int mxsize, nijend, nklend;
    int swaptu;
    int npgtoa, npgtob, npgtoc, npgtod, npsize, ncsize, nwsize,
        shella, shellb, shellc, shelld, shellp, nxyzet, nxyzft, shellq,
        shellt, znorma, znormb, znormc, znormd, zcnorm, zrhoab, zrhocd,
        zgqscr, zsrota, zsrotb, zsrotc, zsrotd;
    double rnabsq, rncdsq, spnorm;
    int zscpqk4, lccsega, lccsegb;
    int blocked;
    int lccsegc, lccsegd, zint2dx, zint2dy, zint2dz;
    int equalab;
    int ncgtoab;
    int equalcd;
    int zcbatch, ncgtocd, nabcoor, ncdcoor, npgtoab, zpbatch,
        mgqijkl, npgtocd;
    int reorder;
    int ncolhrr;
    int mxshell, isnrowa, isnrowb, isnrowc, isnrowd, zpinvhf,
        notmove, zqinvhf, nrothrr, nrowhrr;
    int zpqpinv;

    --icore;
    --zcore;
    --alpha;
    --cc;
    --ccend;
    --ccbeg;
    ftable_dim1 = mgrid - 0 + 1;
    ftable_offset = 0 + ftable_dim1 * 0;
    ftable -= ftable_offset;

/*             ...fix the A,B,C,D labels from the 1,2,3,4 ones. */
/*                Calculate the relevant data for the A,B,C,D batch of */
/*                integrals. */
    lexp1 = 1;
    lexp2 = lexp1 + npgto1;
    lexp3 = lexp2 + npgto2;
    lexp4 = lexp3 + npgto3;
    lcc1 = 1;
    lcc2 = lcc1 + npgto1 * ncgto1;
    lcc3 = lcc2 + npgto2 * ncgto2;
    lcc4 = lcc3 + npgto3 * ncgto3;
    erd__set_abcd_ (&ncgto1, &ncgto2, &ncgto3, &ncgto4,
                    &npgto1, &npgto2, &npgto3, &npgto4,
                    &shell1, &shell2, &shell3, &shell4,
                    &x1, &y1, &z1, &x2, &y2, &z2,
                    &x3, &y3, &z3, &x4, &y4, &z4,
                    &alpha[lexp1], &alpha[lexp2],
                    &alpha[lexp3], &alpha[lexp4],
                    &cc[lcc1], &cc[lcc2],
                    &cc[lcc3], &cc[lcc4],
                    &spheric, &ncgtoa, &ncgtob,
                    &ncgtoc, &ncgtod,
                    &npgtoa, &npgtob, &npgtoc, &npgtod,
                    &shella, &shellb, &shellc, &shelld,
                    &shellp, &shellq, &shellt, &mxshell,
                    &xa, &ya, &za, &xb, &yb, &zb,
                    &xc, &yc, &zc, &xd, &yd, &zd,
                    &atomic, &atomab, &atomcd,
                    &equalab, &equalcd,
                    &abx, &aby, &abz, &cdx, &cdy, &cdz,
                    &nabcoor, &ncdcoor, &rnabsq, &rncdsq,
                    &spnorm, &nxyza, &nxyzb, &nxyzc, &nxyzd,
                    &nxyzet, &nxyzft, &nxyzp, &nxyzq,
                    &nrya, &nryb, &nryc, &nryd,
                    &indexa, &indexb, &indexc, &indexd,
                    &swap12, &swap34, &swaprs, &swaptu, &tr1234,
                    &lexpa, &lexpb, &lexpc, &lexpd,
                    &lcca, &lccb, &lccc, &lccd,
                    &lccsega, &lccsegb, &lccsegc, &lccsegd,
                    &nxyzhrr, &ncolhrr, &nrothrr, &empty);
    if (empty)
    {
        *nbatch = 0;
        return 0;
    }

/*             ...enter the cartesian contracted (e0|f0) batch */
/*                generation. Set the ij and kl primitive exponent */
/*                pairs and the corresponding exponential prefactors. */
    if (equalab)
    {
        npgtoab = npgtoa * (npgtoa + 1) / 2;
        ncgtoab = ncgtoa * (ncgtoa + 1) / 2;
    }
    else
    {
        npgtoab = npgtoa * npgtob;
        ncgtoab = ncgtoa * ncgtob;
    }
    if (equalcd)
    {
        npgtocd = npgtoc * (npgtoc + 1) / 2;
        ncgtocd = ncgtoc * (ncgtoc + 1) / 2;
    }
    else
    {
        npgtocd = npgtoc * npgtod;
        ncgtocd = ncgtoc * ncgtod;
    }
    nctr = ncgtoab * ncgtocd;
    nxyzt = nxyzet * nxyzft;
    iprima = 1;
    iprimb = iprima + npgtoab;
    iprimc = iprimb + npgtoab;
    iprimd = iprimc + npgtocd;
    erd__set_ij_kl_pairs (npgtoa, npgtob, npgtoc, npgtod,
                          atomab, atomcd, equalab, equalcd,
                          swaprs, swaptu, xa, ya, za, xb, yb, zb,
                          xc, yc, zc, xd, yd, zd,
                          rnabsq, rncdsq, PREFACT,
                          &alpha[lexpa], &alpha[lexpb],
                          &alpha[lexpc], &alpha[lexpd],
                          &ftable[ftable_offset], mgrid, ngrid,
                          tmax, tstep, tvstep,
                          screen, &empty, &nij, &nkl,
                          &icore[iprima], &icore[iprimb], &icore[iprimc],
                          &icore[iprimd], &zcore[1]);
    if (empty)
    {
        *nbatch = 0;
        return 0;
    }

/*             ...decide on the primitive [e0|f0] block size and */
/*                return array sizes and pointers for the primitive */
/*                [e0|f0] generation. Perform also some preparation */
/*                steps for contraction. */
    ngqp = shellt / 2 + 1;
    nmom = (ngqp << 1) - 1;
    ngqscr = nmom * 5 + (ngqp << 1) - 2;
    memory = 0;
    erd__e0f0_def_blocks (zmax, npgtoa, npgtob, npgtoc, npgtod,
                          shellp, shellq, nij, nkl,
                          ncgtoab, ncgtocd, nctr, ngqp, ngqscr,
                          nxyzt, l1cache, nctrow, memory,
                          &nijblk, &nklblk, &npsize, &ncsize, &nwsize,
                          &nint2d, &mxprim, &mnprim, &zcbatch, &zpbatch,
                          &zwork, &znorma, &znormb, &znormc, &znormd,
                          &zrhoab, &zrhocd, &zp, &zpx, &zpy, &zpz, &zpax,
                          &zpay, &zpaz, &zpinvhf, &zscpk2,
                          &zq, &zqx, &zqy, &zqz, &zqcx, &zqcy,
                          &zqcz, &zqinvhf, &zscqk2,
                          &zrts, &zwts, &zgqscr, &ztval,
                          &zpqpinv, &zscpqk4, &zb00, &zb01, &zb10,
                          &zc00x, &zc00y, &zc00z,
                          &zd00x, &zd00y, &zd00z,
                          &zint2dx, &zint2dy, &zint2dz);
    blocked = 0;

    erd__prepare_ctr (ncsize, nij, nkl,
                      npgtoa, npgtob, npgtoc, npgtod,
                      shella, shellb, shellc, shelld,
                      &alpha[lexpa], &alpha[lexpb],
                      &alpha[lexpc], &alpha[lexpd],
                      PREFACT, spnorm, equalab, equalcd,
                      blocked, &zcore[1],
                      &zcore[znorma], &zcore[znormb],
                      &zcore[znormc], &zcore[znormd],
                      &zcore[zrhoab], &zcore[zrhocd], &zcore[zcbatch]);
    ipused = iprimd + npgtocd;
    ipsave = ipused + mnprim;
    ippair = ipsave + mxprim;


/*             ...evaluate unnormalized rescaled [e0|f0] in blocks */
/*                over ij and kl pairs and add to final contracted */
/*                (e0|f0). The keyword REORDER indicates, if the */
/*                primitive [e0|f0] blocks need to be transposed */
/*                before being contracted. */
    reorder = 1;
    nijbeg = 1;
    nijend = nij;
    mij = nij;
    nklbeg = 1;
    nklend = nkl;
    mkl = nkl;
    mijkl = mij * mkl;
    mgqijkl = ngqp * mijkl;
    erd__e0f0_pcgto_block_ (&npsize, &nint2d, &atomic, &atomab, &atomcd,
                            &mij, &mkl, &mijkl, &nij, &nijbeg, &nijend, &nkl,
                            &nklbeg, &nklend, &ngqp, &nmom, &ngqscr,
                            &mgqijkl, &npgtoa, &npgtob, &npgtoc, &npgtod,
                            &nxyzet, &nxyzft, &nxyzp, &nxyzq, &shella,
                            &shellp, &shellc, &shellq, &xa, &ya, &za, &xb,
                            &yb, &zb, &xc, &yc, &zc, &xd, &yd, &zd, &abx,
                            &aby, &abz, &cdx, &cdy, &cdz, &alpha[lexpa],
                            &alpha[lexpb], &alpha[lexpc], &alpha[lexpd],
                            &ftable[ftable_offset], &mgrid, &ngrid, &tmax,
                            &tstep, &tvstep, &icore[iprima + nijbeg - 1],
                            &icore[iprimb + nijbeg - 1],
                            &icore[iprimc + nklbeg - 1],
                            &icore[iprimd + nklbeg - 1], &zcore[znorma],
                            &zcore[znormb], &zcore[znormc], &zcore[znormd],
                            &zcore[zrhoab], &zcore[zrhocd], &zcore[zp],
                            &zcore[zpx], &zcore[zpy], &zcore[zpz],
                            &zcore[zpax], &zcore[zpay], &zcore[zpaz],
                            &zcore[zpinvhf], &zcore[zscpk2], &zcore[zq],
                            &zcore[zqx], &zcore[zqy], &zcore[zqz],
                            &zcore[zqcx], &zcore[zqcy], &zcore[zqcz],
                            &zcore[zqinvhf], &zcore[zscqk2], &zcore[zrts],
                            &zcore[zwts], &zcore[zgqscr], &zcore[ztval],
                            &zcore[zpqpinv], &zcore[zscpqk4], &zcore[zb00],
                            &zcore[zb01], &zcore[zb10], &zcore[zc00x],
                            &zcore[zc00y], &zcore[zc00z], &zcore[zd00x],
                            &zcore[zd00y], &zcore[zd00z], &zcore[zint2dx],
                            &zcore[zint2dy], &zcore[zint2dz],
                            &zcore[zpbatch]);

    erd__ctr_4index_block (npsize, ncsize, nwsize, nxyzt, mijkl,
                           mij, mkl, ncgtoab, ncgtocd,
                           npgtoa, npgtob, npgtoc, npgtod,
                           ncgtoa, ncgtob, ncgtoc, ncgtod,
                           mxprim, mnprim,
                           &cc[lcca], &cc[lccb], &cc[lccc],&cc[lccd], 
                           &ccbeg[lccsega], &ccbeg[lccsegb],
                           &ccbeg[lccsegc], &ccbeg[lccsegd],
                           &ccend[lccsega], &ccend[lccsegb],
                           &ccend[lccsegc], &ccend[lccsegd],
                           &icore[iprima + nijbeg - 1],
                           &icore[iprimb + nijbeg - 1],
                           &icore[iprimc + nklbeg - 1],
                           &icore[iprimd + nklbeg - 1],
                           l1cache, tile, nctrow, equalab, equalcd,
                           swaprs, swaptu, reorder, blocked,
                           &icore[ipused], &icore[ipsave],
                           &icore[ippair], &zcore[zpbatch],
                           &zcore[zwork], &zcore[zcbatch]);

/*             ...the unnormalized cartesian (e0|f0) contracted batch is */
/*                ready. Expand the contraction indices (if necessary): */
/*                   batch (nxyzt,r>=s,t>=u) --> batch (nxyzt,r,s,t,u) */
/*                and reorder the contraction index part (if necessary): */
/*                   batch (nxyzt,r,s,t,u) --> batch (nxyzt,1,2,3,4) */
/*                The array IXOFF (x) indicates the total # of indices */
/*                to the left of x without including the nxyzt-part. */
/*                For the left most IXOFF value it is convenient to */
/*                set it equal to 1 instead of 0. Note, that the IXOFF */
/*                array indicates the true # of indices to the left */
/*                after! the batch has been transposed (see below) and */
/*                can be used as initial values when moving the */
/*                ry-components later on during the HRR and cartesian -> */
/*                spherical transformation procedure. */
/*                For efficient application of the HRR contraction */
/*                scheme we need the batch elements ordered as: */
/*                          batch (1,2,3,4,nxyzt) */
/*                hence we transpose the batch after the reordering. */
/*                The space partitioning of the flp array for all of */
/*                these steps will be as follows: */
/*                          |  Zone 1  |  Zone 2  | */
/*                in which Zone 1 and 2 are 2 batches of HRR maximum */
/*                size. This can be done because we have always the */
/*                following dimension inequality: */
/*                              NXYZT =< NXYZHRR */

    ixoff[0] = 1;
    ixoff[1] = ncgto1;
    ixoff[2] = ncgto1 * ncgto2;
    ixoff[3] = ixoff[2] * ncgto3;
    nctr = ixoff[3] * ncgto4;
    mxsize = nctr * nxyzhrr;
    in = zcbatch;
    out = in + mxsize;
#if 0    
    if (equalab && ncgtoab > 1)
    {
        erd__ctr_rs_expand__ (&nxyzt, &ncgtoab, &ncgtocd, &ncgtoa, &ncgtob,
                              &zcore[in], &zcore[out]);
        temp = in;
        in = out;
        out = temp;
    }
    if (equalcd && ncgtocd > 1)
    {
        int ii;
        ii = nxyzt * ncgtoa * ncgtob;
        erd__ctr_tu_expand__ (&ii, &ncgtocd, &ncgtoc, &ncgtod, &zcore[in],
                              &zcore[out]);
        temp = in;
        in = out;
        out = temp;
    }
    reorder = tr1234 || swap12 != swaprs || swap34 != swaptu;
    if (reorder && nctr > 1)
    {
        if (swaprs)
        {
            indexr = indexb;
            indexs = indexa;
            ncgtor = ncgtob;
            ncgtos = ncgtoa;
        }
        else
        {
            indexr = indexa;
            indexs = indexb;
            ncgtor = ncgtoa;
            ncgtos = ncgtob;
        }
        if (swaptu)
        {
            indext = indexd;
            indexu = indexc;
            ncgtot = ncgtod;
            ncgtou = ncgtoc;
        }
        else
        {
            indext = indexc;
            indexu = indexd;
            ncgtot = ncgtoc;
            ncgtou = ncgtod;
        }
        erd__ctr_4index_reorder__ (&nxyzt, &nctr, &ncgtor, &ncgtos, &ncgtot,
                                   &ncgtou, &ixoff[indexr - 1],
                                   &ixoff[indexs - 1], &ixoff[indext - 1],
                                   &ixoff[indexu - 1], &zcore[in],
                                   &zcore[out]);
        temp = in;
        in = out;
        out = temp;
    }
#endif
    
    if (nxyzt > 1 && nctr > 1)
    {
        erd__transpose_batch (nxyzt, nctr, &zcore[in], &zcore[out]);
        temp = in;
        in = out;
        out = temp;
    }
    
/*             ...enter the HRR contraction and cartesian -> spherical */
/*                transformation / cartesian normalization section. */
/*                The sequence of events is to apply the HRR followed */
/*                by cartesian -> spherical transformations or cartesian */
/*                normalizations and immediate final positioning of */
/*                the finished parts to correspond with the contraction */
/*                indices i,j,k,l. First we do the f0-part followed */
/*                by the e0-part, which hence gives rise to the sequence */
/*                (where ' means spherical or cartesian normalization */
/*                and [] means the indices are in correspondence): */
/*                   batch (ijkl,e0,f0) --> batch (ijkl,e0,cd) */
/*                   batch (ijkl,e0,c,d) --> batch (ijkl,e0,c,d') */
/*                   batch (ijkl,e0,c,d') --> batch (ijkl[d'],e0,c) */
/*                   batch (ijkl[d'],e0,c) --> batch (ijkl[d'],e0,c') */
/*                   batch (ijkl[d'],e0,c') --> batch (ijkl[c'd'],e0) */
/*                   batch (ijkl[c'd'],e0) --> batch (ijkl[c'd'],ab) */
/*                   batch (ijkl[c'd'],a,b) --> batch (ijkl[c'd'],a,b') */
/*                   batch (ijkl[c'd'],a,b') --> batch (ijkl[b'c'd'],a) */
/*                   batch (ijkl[b'c'd'],a) --> batch (ijkl[b'c'd'],a') */
/*                   batch (ijkl[b'c'd'],a') --> batch (ijkl[a'b'c'd']) */
/*                The space partitioning of the flp array will be */
/*                as follows: */
/*                  |  Zone 1  |  Zone 2  |  Zone 3  |  Zone 4  | */
/*                 Zone 1 and 2:  2 batches of MXSIZE maximum size */
/*                                (set previously) */
/*                       Zone 3:  cart -> spher transformation data */
/*                                             or */
/*                                cartesian normalization factors */
/*                       Zone 4:  HRR contraction data */
/*                Determine memory allocation offsets for the entire HRR */
/*                procedure + cartesian -> spherical transformations or */
/*                cartesian normalizations and generate the transformation */
/*                matrices + associated data for those shells > p-shell. */
/*                The offsets are as follows (x=A,B,C,D): */
/*                    IN = offset for input HRR batch */
/*                   OUT = offset for output HRR batch */
/*                ZSROTx = offset for x-part transformation matrix */
/*               ISNROWx = offset for # of non-zero XYZ contribution row */
/*                         labels for x-part transformation matrix */
/*                ISROWx = offset for non-zero XYZ contribution row */
/*                         labels for x-part transformation matrix */
/*                 ZHROT = offset for HRR transformation matrix */
/*                IHNROW = offset for # of nonzero row labels for */
/*                         each HRR matrix column */
/*                 IHROW = offset for nonzero row labels for the HRR */
/*                         matrix */
/*                 IHSCR = int scratch space for HRR matrix */
/*                         assembly */
/*                In case of s- or p-shells no transformation matrix is */
/*                generated, hence if we have s- and/or p-shells, then */
/*                no call to the cartesian -> spherical transformation */
/*                or cartesian normalization routines needs to be done. */
/*                All integrals have already been multiplied by a factor */
/*                SPNORM, which has the following value for each s- and */
/*                p-shell: */
/*                       For s-type shell  =  1 */
/*                       For p-type shell  =  2 * norm for s-type */
/*                This factor was introduced together with the overall */
/*                prefactor during evaluation of the primitive integrals */
/*                in order to save multiplications. */
    zbase = MAX (in, out) + mxsize;
    if (spheric)
    {
        if (mxshell > 1)
        {
            int c__1 = 1;
            erd__xyz_to_ry_abcd (nxyza, nxyzb, nxyzc, nxyzd,
                                 nrya, nryb, nryc, nryd,
                                 shella, shellb, shellc, shelld,
                                 1, zbase,
                                 &nrowa, &nrowb, &nrowc, &nrowd,
                                 &nrota, &nrotb, &nrotc, &nrotd,
                                 &zsrota, &zsrotb, &zsrotc, &zsrotd,
                                 &isnrowa, &isnrowb, &isnrowc, &isnrowd,
                                 &isrowa, &isrowb, &isrowc, &isrowd,
                                 &iused, &zused, &icore[1], &zcore[1]);
        }
        else
        {
            iused = 0;
            zused = 0;
        }
    }
    else
    {
        if (mxshell > 1)
        {
            zcnorm = zbase;
            erd__cartesian_norms_ (&mxshell, &zcore[zcnorm]);
            iused = 0;
            zused = mxshell + 1;
        }
        else
        {
            iused = 0;
            zused = 0;
        }
    }
    ihnrow = iused + 1;
    ihrow = ihnrow + ncolhrr + ncolhrr;
    ihscr = ihrow + nrothrr + nrothrr;
    zhrot = zbase + zused;


/*             ...do the first stage of processing the integrals: */
/*                   batch (ijkl,e0,f0) --> batch (ijkl,e0,cd) */
/*                   batch (ijkl,e0,c,d) --> batch (ijkl,e0,c,d') */
/*                   batch (ijkl,e0,c,d') --> batch (ijkl[d'],e0,c) */
/*                   batch (ijkl[d'],e0,c) --> batch (ijkl[d'],e0,c') */
/*                   batch (ijkl[d'],e0,c') --> batch (ijkl[c'd'],e0) */
    if (shelld != 0)
    {
        erd__hrr_matrix (nrothrr, ncolhrr, nxyzft, nxyzc, nxyzq,
                         shellc, shelld, shellq,
                         ncdcoor, cdx, cdy, cdz,
                         &icore[ihscr], &pos1, &pos2, &nrowhrr,
                         &icore[ihnrow], &icore[ihrow], &zcore[zhrot]);
        erd__hrr_transform (nctr * nxyzet, nrowhrr, nxyzft,
                            nxyzc * nxyzd, nxyzc, nxyzd,
                            &icore[ihnrow + pos1 - 1],
                            &icore[ihrow + pos2 - 1],
                            &zcore[zhrot + pos2 - 1], &zcore[in],
                            &zcore[out]);
        temp = in;
        in = out;
        out = temp;
        if (shelld > 1)
        {
            if (spheric)
            {
                erd__spherical_transform (nctr * nxyzet * nxyzc,
                                          nrowd, nxyzd, nryd,
                                          &icore[isnrowd], &icore[isrowd],
                                          &zcore[zsrotd], &zcore[in],
                                          &zcore[out]);
                temp = in;
                in = out;
                out = temp;
            }
            else
            {
                int i__1;
                i__1 = nctr * nxyzet * nxyzc;
                erd__normalize_cartesian_ (&i__1, &nxyzd, &shelld,
                                           &zcore[zcnorm], &zcore[in]);
            }
        }
    }
    if (nryd > 1)
    {
        *nbatch = nctr * nxyzet * nxyzc * nryd;
        notmove = ixoff[indexd - 1];
        move = *nbatch / (notmove * nryd);
        if (move > 1)
        {
            erd__move_ry (*nbatch, 4, notmove, move, nryd, indexd,
                          &zcore[in], ixoff, &zcore[out]);
            temp = in;
            in = out;
            out = temp;
        }
    }
    if (shellc > 1)
    {
        if (spheric)
        {
            erd__spherical_transform (nctr * nxyzet * nryd,
                                      nrowc, nxyzc, nryc,
                                      &icore[isnrowc], &icore[isrowc],
                                      &zcore[zsrotc], &zcore[in],
                                      &zcore[out]);
            temp = in;
            in = out;
            out = temp;
        }
        else
        {
            int i__1;
            i__1 = nctr * nxyzet * nryd;
            erd__normalize_cartesian_ (&i__1, &nxyzc, &shellc,
                                        &zcore[zcnorm], &zcore[in]);
        }
    }
    if (nryc > 1)
    {
        *nbatch = nctr * nryd * nxyzet * nryc;
        notmove = ixoff[indexc - 1];
        move = *nbatch / (notmove * nryc);
        if (move > 1)
        {
            erd__move_ry (*nbatch, 4, notmove, move, nryc, indexc,
                           &zcore[in], ixoff, &zcore[out]);
            temp = in;
            in = out;
            out = temp;
        }
    }

/*             ...do the second stage of processing the integrals: */
/*                   batch (ijkl[c'd'],e0) --> batch (ijkl[c'd'],ab) */
/*                   batch (ijkl[c'd'],a,b) --> batch (ijkl[c'd'],a,b') */
/*                   batch (ijkl[c'd'],a,b') --> batch (ijkl[b'c'd'],a) */
/*                   batch (ijkl[b'c'd'],a) --> batch (ijkl[b'c'd'],a') */
/*                   batch (ijkl[b'c'd'],a') --> batch (ijkl[a'b'c'd']) */
    if (shellb != 0)
    {
        erd__hrr_matrix (nrothrr, ncolhrr, nxyzet, nxyza, nxyzp,
                         shella, shellb, shellp,
                         nabcoor, abx, aby, abz,
                         &icore[ihscr], &pos1, &pos2, &nrowhrr,
                         &icore[ihnrow], &icore[ihrow], &zcore[zhrot]);
        erd__hrr_transform (nctr * nryc * nryd, nrowhrr, nxyzet,
                            nxyza * nxyzb, nxyza, nxyzb,
                            &icore[ihnrow + pos1 - 1],
                            &icore[ihrow + pos2 - 1],
                            &zcore[zhrot + pos2 - 1], &zcore[in],
                            &zcore[out]);
        temp = in;
        in = out;
        out = temp;
        if (shellb > 1)
        {
            if (spheric)
            {
                erd__spherical_transform (nctr * nryc * nryd * nxyza,
                                          nrowb, nxyzb, nryb,
                                          &icore[isnrowb], &icore[isrowb],
                                          &zcore[zsrotb], &zcore[in],
                                          &zcore[out]);
                temp = in;
                in = out;
                out = temp;
            }
            else
            {
                int i__1;
                i__1 = nctr * nryc * nryd * nxyza;
                erd__normalize_cartesian_ (&i__1, &nxyzb, &shellb,
                                           &zcore[zcnorm], &zcore[in]);
            }
        }
    }
    if (nryb > 1)
    {
        *nbatch = nctr * nryc * nryd * nxyza * nryb;
        notmove = ixoff[indexb - 1];
        move = *nbatch / (notmove * nryb);
        if (move > 1)
        {
            erd__move_ry (*nbatch, 4, notmove, move, nryb, indexb,
                          &zcore[in], ixoff, &zcore[out]);
            temp = in;
            in = out;
            out = temp;
        }
    }
    if (shella > 1)
    {
        if (spheric)
        {
            erd__spherical_transform (nctr * nryb * nryc * nryd,
                                      nrowa, nxyza, nrya,
                                      &icore[isnrowa], &icore[isrowa],
                                      &zcore[zsrota], &zcore[in],
                                      &zcore[out]);
            temp = in;
            in = out;
            out = temp;
        }
        else
        {
            int i__1;
            i__1 = nctr * nryb * nryc * nryd;
            erd__normalize_cartesian_ (&i__1, &nxyza, &shella,
                                        &zcore[zcnorm], &zcore[in]);
        }
    }
    *nbatch = nctr * nryb * nryc * nryd * nrya;
    if (nrya > 1)
    {
        notmove = ixoff[indexa - 1];
        move = *nbatch / (notmove * nrya);
        if (move > 1)
        {
            erd__move_ry (*nbatch, 4, notmove, move, nrya, indexa,
                          &zcore[in], ixoff, &zcore[out]);
            temp = in;
            in = out;
            out = temp;
        }
    }


/*             ...set final pointer to integrals in ZCORE array. */
    *nfirst = in;

    return 0;
}


int erd__csgto_ (int * imax, int * zmax, int * nalpha,
              int * ncoeff, int * ncsum, int * ncgto1,
              int * ncgto2, int * ncgto3, int * ncgto4,
              int * npgto1, int * npgto2, int * npgto3,
              int * npgto4, int * shell1, int * shell2,
              int * shell3, int * shell4, double * x1,
              double * y1, double * z1, double * x2,
              double * y2, double * z2, double * x3,
              double * y3, double * z3, double * x4,
              double * y4, double * z4, double * alpha,
              double * cc, int * ccbeg, int * ccend,
              double * ftable, int * mgrid, int * ngrid,
              double * tmax, double * tstep, double * tvstep,
              int * l1cache, int * tile, int * nctrow,
              int * spheric, int * screen, int * icore,
              int * nbatch, int * nfirst, double * zcore)
{
    erd__csgto (*imax, *zmax,
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
                alpha, cc,
                ccbeg, ccend,
                ftable, *mgrid, *ngrid,
                *tmax, *tstep, *tvstep,
                *l1cache, *tile, *nctrow,
                *spheric, *screen,
                icore, nbatch, nfirst, zcore);
                
    return 0;
}
