#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "erd.h"


#pragma offload_attribute(push, target(mic))

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
int erd__csgto (int zmax, int npgto1, int npgto2,
                int npgto3, int npgto4,
                int shell1, int shell2,
                int shell3, int shell4,
                double x1, double y1, double z1,
                double x2, double y2, double z2,
                double x3, double y3, double z3,
                double x4, double y4, double z4,
                double *alpha1, double *alpha2,
                double *alpha3, double *alpha4,
                double *cc1, double *cc2,
                double *cc3, double *cc4,
                int spheric, int screen, int *icore,
                int *nbatch, int *nfirst, double *zcore)
{
    int nxyzhrr;   
    int in;
    double xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd;
    int zp, zq;
    int zb00, zb01, zb10;
    int nij;
    int nkl;
    int out, zpx, zpy, zpz, zqx, zqy, zqz,
        pos1, pos2, lcca, lccb, lccc, lccd;
    int zc00x, ngqp, move, nmom, nrya, nryb, nryc, nryd,
        temp, zc00y, zc00z, zpax, zpay, zpaz, zqcx, zqcy, zqcz, zd00x,
        zd00y, zd00z;
    int zrts, zwts;
    int zbase;
    int ihscr, iused, lexpa, lexpb, lexpc;
    int lexpd, ixoff[4];
    int nrota, nrotb, nrotc, nrotd, ihrow, nrowa, zused, nrowb,
        nrowc, nrowd;
    int empty;
    int ztval, zhrot, nxyza, nxyzb, nxyzc, nxyzd, nint2d,
        nxyzp, nxyzq, nxyzt, zscpk2, zscqk2;
    int indexa, indexb, indexc, indexd;
    int iprima, iprimb, iprimc, iprimd,
        isrowa, ihnrow,
        isrowb;
    int isrowc, isrowd, ngqscr;
    int mxsize;
    int npgtoa, npgtob, npgtoc, npgtod,
        shella, shellb, shellc, shelld, shellp, nxyzet, nxyzft, shellq,
        shellt, znorma, znormb, znormc, znormd, zcnorm, zrhoab, zrhocd,
        zgqscr, zsrota, zsrotb, zsrotc, zsrotd;
    double rnabsq, rncdsq, spnorm;
    int zscpqk4;
    int zint2dx, zint2dy, zint2dz;
    int zcbatch, nabcoor, ncdcoor, npgtoab,
        npgtocd;
    int ncolhrr;
    int mxshell, isnrowa, isnrowb, isnrowc, isnrowd, zpinvhf,
        notmove, zqinvhf, nrothrr, nrowhrr;
    int zpqpinv;
    double abx, aby, abz, cdx, cdy, cdz;
    double *alphaa;
    double *alphab;
    double *alphac;
    double *alphad;
    double *cca;
    double *ccb;
    double *ccc;
    double *ccd;
    
    
#ifdef __ERD_PROFILE__
    uint64_t start_clock, end_clock;
    int tid = omp_get_thread_num();
#endif    
    --icore;
    --zcore;
/*             ...fix the A,B,C,D labels from the 1,2,3,4 ones. */
/*                Calculate the relevant data for the A,B,C,D batch of */
/*                integrals. */
    int tr1234;
    erd__set_abcd (npgto1, npgto2, npgto3, npgto4,
                   shell1, shell2, shell3, shell4,
                   x1, y1, z1, x2, y2, z2,
                   x3, y3, z3, x4, y4, z4, spheric,
                   &npgtoa, &npgtob, &npgtoc, &npgtod,
                   &shella, &shellb, &shellc, &shelld,
                   &xa, &ya, &za, &xb, &yb, &zb,
                   &xc, &yc, &zc, &xd, &yd, &zd,
                   &nxyza, &nxyzb, &nxyzc, &nxyzd,
                   &nxyzet, &nxyzft,
                   &nrya, &nryb, &nryc, &nryd,
                   &indexa, &indexb, &indexc, &indexd,
                   &lexpa, &lexpb, &lexpc, &lexpd,
                   &lcca, &lccb, &lccc, &lccd,
                   &nabcoor, &ncdcoor,
                   &ncolhrr, &nrothrr, &nxyzhrr, &empty, &tr1234);
    if (empty)
    {
        *nbatch = 0;
        return 0;
    }

    int swap12 = (shell1 < shell2);
    int swap34 = (shell3 < shell4);
    
    if (!tr1234)
    {
        if (!swap12)
        {
            alphaa = alpha1;
            alphab = alpha2;
            cca = cc1;
            ccb = cc2;
        }
        else
        {
            alphaa = alpha2;
            alphab = alpha1;
            cca = cc2;
            ccb = cc1;
        }
        if (!swap34)
        {
            alphac = alpha3;
            alphad = alpha4;
            ccc = cc3;
            ccd = cc4;
        }
        else
        {
            alphac = alpha4;
            alphad = alpha3;
            ccc = cc4;
            ccd = cc3;
        }
    }
    else
    {
        if (!swap12)
        {
            alphac = alpha1;
            alphad = alpha2;
            ccc = cc1;
            ccd = cc2;
        }
        else
        {
            alphac = alpha2;
            alphad = alpha1;
            ccc = cc2;
            ccd = cc1;
        }
        if (!swap34)
        {
            alphaa = alpha3;
            alphab = alpha4;
            cca = cc3;
            ccb = cc4;
        }
        else
        {
            alphaa = alpha4;
            alphab = alpha3;
            cca = cc4;
            ccb = cc3;
        }
    }
    
    // initialize values
    abx = xa - xb;
    aby = ya - yb;
    abz = za - zb;
    cdx = xc - xd;
    cdy = yc - yd;
    cdz = zc - zd;
/*             ...the new A,B,C,D shells are set. Calculate the */
/*                following info: 1) control variables to be used */
/*                during contraction, 2) total shell values for */
/*                electrons 1 and 2 in [AB|CD], 3) their corresponding */
/*                cartesian monomial sizes and 4) the overall norm */
/*                factor SPNORM due to presence of s- or p-type shells. */
/*                The latter is necessary, because for such shells */
/*                there will be no calls to the cartesian normalization */
/*                or spherical transformation routines. The contribution */
/*                to SPNORM is very simple: each s-type shell -> * 1.0, */
/*                each p-type shell -> * 2.0. */
    shellp = shella + shellb;
    shellq = shellc + shelld;
    shellt = shellp + shellq;
    mxshell = MAX(shell1, shell2);
    mxshell = MAX(mxshell, shell3);
    mxshell = MAX(mxshell, shell4);
    nxyzp = (shellp + 1) * (shellp + 2) / 2;
    nxyzq = (shellq + 1) * (shellq + 2) / 2;
    // spnorm
    spnorm = 1.;
    if (shella == 1)
    {
        spnorm += spnorm;
    }
    if (shellb == 1)
    {
        spnorm += spnorm;
    }
    if (shellc == 1)
    {
        spnorm += spnorm;
    }
    if (shelld == 1)
    {
        spnorm += spnorm;
    }
    rnabsq = abx * abx + aby * aby + abz * abz;
    rncdsq = cdx * cdx + cdy * cdy + cdz * cdz;

   
/*             ...enter the cartesian contracted (e0|f0) batch */
/*                generation. Set the ij and kl primitive exponent */
/*                pairs and the corresponding exponential prefactors. */
    npgtoab = npgtoa * npgtob;
    npgtocd = npgtoc * npgtod;
    nxyzt = nxyzet * nxyzft;
    iprima = 1;
    iprimb = iprima + npgtoab;
    iprimc = iprimb + npgtoab;
    iprimd = iprimc + npgtocd;

#ifdef __ERD_PROFILE__
    start_clock = __rdtsc();
#endif
    erd__set_ij_kl_pairs (npgtoa, npgtob, npgtoc, npgtod,
                          xa, ya, za, xb, yb, zb,
                          xc, yc, zc, xd, yd, zd,
                          rnabsq, rncdsq, PREFACT,
                          alphaa, alphab,
                          alphac, alphad,
                          screen, &empty, &nij, &nkl,
                          &icore[iprima], &icore[iprimb], &icore[iprimc],
                          &icore[iprimd], &zcore[1]);
#ifdef __ERD_PROFILE__
    end_clock = __rdtsc(); 
    erd_ticks[tid][erd__set_ij_kl_pairs_ticks] += (end_clock - start_clock);
#endif

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
    erd__e0f0_def_blocks (zmax, npgtoa, npgtob, npgtoc, npgtod,
                          shellp, shellq, nij, nkl,
                          ngqp, ngqscr, nxyzt, 0,
                          &nint2d, &zcbatch,
                          &znorma, &znormb, &znormc, &znormd,
                          &zrhoab, &zrhocd, &zp, &zpx, &zpy, &zpz, &zpax,
                          &zpay, &zpaz, &zpinvhf, &zscpk2,
                          &zq, &zqx, &zqy, &zqz, &zqcx, &zqcy,
                          &zqcz, &zqinvhf, &zscqk2,
                          &zrts, &zwts, &zgqscr, &ztval,
                          &zpqpinv, &zscpqk4, &zb00, &zb01, &zb10,
                          &zc00x, &zc00y, &zc00z,
                          &zd00x, &zd00y, &zd00z,
                          &zint2dx, &zint2dy, &zint2dz);
#ifdef __ERD_PROFILE__
    start_clock = __rdtsc();
#endif    
    erd__prepare_ctr (npgtoa, npgtob, npgtoc, npgtod,
                      shella, shellb, shellc, shelld,
                      alphaa, alphab,
                      alphac, alphad, spnorm,
                      &zcore[znorma], &zcore[znormb],
                      &zcore[znormc], &zcore[znormd]);
#ifdef __ERD_PROFILE__
    end_clock = __rdtsc();
    erd_ticks[tid][erd__prepare_ctr_ticks] += (end_clock - start_clock);
#endif
/*             ...evaluate unnormalized rescaled [e0|f0] in blocks */
/*                over ij and kl pairs and add to final contracted */
/*                (e0|f0). The keyword REORDER indicates, if the */
/*                primitive [e0|f0] blocks need to be transposed */
/*                before being contracted. */
#ifdef __ERD_PROFILE__
    start_clock = __rdtsc();
#endif
    erd__e0f0_pcgto_block (nij, nkl, ngqp, nmom,
                           nxyzet, nxyzft, nxyzp, nxyzq,
                           shella, shellp, shellc, shellq,
                           xa, ya, za, xb, yb, zb,
                           xc, yc, zc, xd, yd, zd,
                           alphaa, alphab, alphac, alphad, 
                           cca, ccb, ccc, ccd, 
                           &icore[iprima], &icore[iprimb],
                           &icore[iprimc], &icore[iprimd],
                           &zcore[znorma], &zcore[znormb],
                           &zcore[znormc], &zcore[znormd],
                           &zcore[zrhoab], &zcore[zrhocd],
                           &zcore[zp], &zcore[zpx], &zcore[zpy], &zcore[zpz],
                           &zcore[zpax], &zcore[zpay], &zcore[zpaz],
                           &zcore[zpinvhf], &zcore[zscpk2],
                           &zcore[zq], &zcore[zqx], &zcore[zqy], &zcore[zqz],
                           &zcore[zqcx], &zcore[zqcy], &zcore[zqcz],
                           &zcore[zqinvhf], &zcore[zscqk2],
                           &zcore[zrts], &zcore[zwts],
                           &zcore[zgqscr], &zcore[ztval],
                           &zcore[zpqpinv], &zcore[zscpqk4],
                           &zcore[zb00], &zcore[zb01], &zcore[zb10],
                           &zcore[zc00x], &zcore[zc00y], &zcore[zc00z],
                           &zcore[zd00x], &zcore[zd00y], &zcore[zd00z],
                           &zcore[zint2dx], &zcore[zint2dy], &zcore[zint2dz],
                           &zcore[zcbatch]);
#ifdef __ERD_PROFILE__
    end_clock = __rdtsc();
    erd_ticks[tid][erd__e0f0_pcgto_block_ticks] += (end_clock - start_clock);
#endif
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
    ixoff[1] = 1;
    ixoff[2] = 1;
    ixoff[3] = 1;
    mxsize = nxyzhrr;
    in = zcbatch;
    out = in + mxsize;
   
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
        #ifdef __ERD_PROFILE__
            start_clock = __rdtsc();
        #endif
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
        #ifdef __ERD_PROFILE__
            end_clock = __rdtsc();
            erd_ticks[tid][erd__xyz_to_ry_abcd_ticks] += (end_clock - start_clock);
        #endif
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
            erd__cartesian_norms (mxshell, &zcore[zcnorm]);
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
    #ifdef __ERD_PROFILE__
        start_clock = __rdtsc();
    #endif
        erd__hrr_matrix (nrothrr, ncolhrr, nxyzft, nxyzc, nxyzq,
                         shellc, shelld, shellq,
                         ncdcoor, cdx, cdy, cdz,
                         &icore[ihscr], &pos1, &pos2, &nrowhrr,
                         &icore[ihnrow], &icore[ihrow], &zcore[zhrot]);
    #ifdef __ERD_PROFILE__
        end_clock = __rdtsc();
        erd_ticks[tid][erd__hrr_matrix_ticks] += (end_clock - start_clock);
    #endif

    #ifdef __ERD_PROFILE__
        start_clock = __rdtsc();
    #endif
        erd__hrr_transform (nxyzet, nrowhrr, nxyzc, nxyzd,
                            &icore[ihnrow + pos1 - 1],
                            &icore[ihrow + pos2 - 1],
                            &zcore[zhrot + pos2 - 1], &zcore[in],
                            &zcore[out]);
    #ifdef __ERD_PROFILE__
        end_clock = __rdtsc();
        erd_ticks[tid][erd__hrr_transform_ticks] += (end_clock - start_clock);
    #endif
    
        temp = in;
        in = out;
        out = temp;
        if (shelld > 1)
        {
            if (spheric)
            {
            #ifdef __ERD_PROFILE__
                start_clock = __rdtsc();
            #endif
                erd__spherical_transform (nxyzet * nxyzc,
                                          nrowd, nryd,
                                          &icore[isnrowd], &icore[isrowd],
                                          &zcore[zsrotd], &zcore[in],
                                          &zcore[out]);
            #ifdef __ERD_PROFILE__
                end_clock = __rdtsc();
                erd_ticks[tid][erd__spherical_transform_ticks] += (end_clock - start_clock);
            #endif
                temp = in;
                in = out;
                out = temp;
            }
            else
            {
                erd__normalize_cartesian (nxyzet * nxyzc, shelld,
                                          &zcore[zcnorm], &zcore[in]);
            }
        }
    }
    if (nryd > 1)
    {
        *nbatch = nxyzet * nxyzc * nryd;
        notmove = ixoff[indexd - 1];
        move = *nbatch / (notmove * nryd);
        if (move > 1)
        {
        #ifdef __ERD_PROFILE__
            start_clock = __rdtsc();
        #endif
            erd__move_ry (4, notmove, move, nryd, indexd - 1,
                          &zcore[in], ixoff, &zcore[out]);
        #ifdef __ERD_PROFILE__
            end_clock = __rdtsc();
            erd_ticks[tid][erd__move_ry_ticks] += (end_clock - start_clock);
        #endif
            temp = in;
            in = out;
            out = temp;
        }
    } 
    if (shellc > 1)
    {
        if (spheric)
        {
        #ifdef __ERD_PROFILE__
                start_clock = __rdtsc();
        #endif
            erd__spherical_transform (nxyzet * nryd,
                                      nrowc, nryc,
                                      &icore[isnrowc], &icore[isrowc],
                                      &zcore[zsrotc], &zcore[in],
                                      &zcore[out]);
        #ifdef __ERD_PROFILE__
                end_clock = __rdtsc();
                erd_ticks[tid][erd__spherical_transform_ticks] += (end_clock - start_clock);
        #endif
            temp = in;
            in = out;
            out = temp;
        }
        else
        {
            erd__normalize_cartesian (nxyzet * nryd, shellc,
                                      &zcore[zcnorm], &zcore[in]);
        }
    }
    if (nryc > 1)
    {
        *nbatch = nryd * nxyzet * nryc;
        notmove = ixoff[indexc - 1];
        move = *nbatch / (notmove * nryc);
        if (move > 1)
        {
        #ifdef __ERD_PROFILE__
            start_clock = __rdtsc();
        #endif
            erd__move_ry (4, notmove, move, nryc, indexc - 1,
                          &zcore[in], ixoff, &zcore[out]);
        #ifdef __ERD_PROFILE__
            end_clock = __rdtsc();
            erd_ticks[tid][erd__move_ry_ticks] += (end_clock - start_clock);
        #endif
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
    #ifdef __ERD_PROFILE__
        start_clock = __rdtsc();
    #endif
        erd__hrr_matrix (nrothrr, ncolhrr, nxyzet, nxyza, nxyzp,
                         shella, shellb, shellp,
                         nabcoor, abx, aby, abz,
                         &icore[ihscr], &pos1, &pos2, &nrowhrr,
                         &icore[ihnrow], &icore[ihrow], &zcore[zhrot]);
    #ifdef __ERD_PROFILE__
        end_clock = __rdtsc();
        erd_ticks[tid][erd__hrr_matrix_ticks] += (end_clock - start_clock);
    #endif

    #ifdef __ERD_PROFILE__
        start_clock = __rdtsc();
    #endif
        erd__hrr_transform (nryc * nryd, nrowhrr, nxyza, nxyzb,
                            &icore[ihnrow + pos1 - 1],
                            &icore[ihrow + pos2 - 1],
                            &zcore[zhrot + pos2 - 1], &zcore[in],
                            &zcore[out]);
    #ifdef __ERD_PROFILE__
        end_clock = __rdtsc();
        erd_ticks[tid][erd__hrr_transform_ticks] += (end_clock - start_clock);
    #endif
        temp = in;
        in = out;
        out = temp;
        if (shellb > 1)
        {
            if (spheric)
            {
            #ifdef __ERD_PROFILE__
                start_clock = __rdtsc();
            #endif
                erd__spherical_transform (nryc * nryd * nxyza,
                                          nrowb, nryb,
                                          &icore[isnrowb], &icore[isrowb],
                                          &zcore[zsrotb], &zcore[in],
                                          &zcore[out]);
            #ifdef __ERD_PROFILE__
                end_clock = __rdtsc();
                erd_ticks[tid][erd__spherical_transform_ticks] += (end_clock - start_clock);
            #endif
                temp = in;
                in = out;
                out = temp;
            }
            else
            {
                erd__normalize_cartesian (nryc * nryd * nxyza, shellb,
                                          &zcore[zcnorm], &zcore[in]);
            }
        }
    }
    if (nryb > 1)
    {
        *nbatch = nryc * nryd * nxyza * nryb;
        notmove = ixoff[indexb - 1];
        move = *nbatch / (notmove * nryb);
        if (move > 1)
        {
        #ifdef __ERD_PROFILE__
            start_clock = __rdtsc();
        #endif
            erd__move_ry (4, notmove, move, nryb, indexb - 1,
                          &zcore[in], ixoff, &zcore[out]);
        #ifdef __ERD_PROFILE__
            end_clock = __rdtsc();
            erd_ticks[tid][erd__move_ry_ticks] += (end_clock - start_clock);
        #endif
            temp = in;
            in = out;
            out = temp;
        }
    }
    if (shella > 1)
    {
        if (spheric)
        {
        #ifdef __ERD_PROFILE__
            start_clock = __rdtsc();
        #endif
            erd__spherical_transform (nryb * nryc * nryd,
                                      nrowa, nrya,
                                      &icore[isnrowa], &icore[isrowa],
                                      &zcore[zsrota], &zcore[in],
                                      &zcore[out]);

        #ifdef __ERD_PROFILE__
            end_clock = __rdtsc();
            erd_ticks[tid][erd__spherical_transform_ticks] += (end_clock - start_clock);
        #endif
            temp = in;
            in = out;
            out = temp;
        }
        else
        {
            erd__normalize_cartesian (nryb * nryc * nryd, shella,
                                      &zcore[zcnorm], &zcore[in]);
        }
    }
    *nbatch = nryb * nryc * nryd * nrya;
    if (nrya > 1)
    {
        notmove = ixoff[indexa - 1];
        move = *nbatch / (notmove * nrya);
        if (move > 1)
        {
        #ifdef __ERD_PROFILE__
            start_clock = __rdtsc();
        #endif
            erd__move_ry (4, notmove, move, nrya, indexa - 1,
                          &zcore[in], ixoff, &zcore[out]);
        #ifdef __ERD_PROFILE__
            end_clock = __rdtsc();
            erd_ticks[tid][erd__move_ry_ticks] += (end_clock - start_clock);
        #endif
            temp = in;
            in = out;
            out = temp;
        }
    }


/*             ...set final pointer to integrals in ZCORE array. */
    *nfirst = in;

    return 0;
}

#pragma offload_attribute(pop)
