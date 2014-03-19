#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "erd.h"
#include "erdutil.h"

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
/*                    ICORE        =  int output_buffer space */
/*                    ZCORE (part) =  flp output_buffer space */
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
ERD_OFFLOAD void erd__csgto(
    bool atomic,
    uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
    uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double alpha1[restrict static npgto1], const double alpha2[restrict static npgto2], const double alpha3[restrict static npgto3], const double alpha4[restrict static npgto4],
    const double cc1[restrict static npgto1], const double cc2[restrict static npgto2], const double cc3[restrict static npgto3], const double cc4[restrict static npgto4],
    const double norm1[restrict static npgto1], const double norm2[restrict static npgto2], const double norm3[restrict static npgto3], const double norm4[restrict static npgto4],
    int **vrrtab, int ldvrrtab,
    bool spheric,
    uint32_t buffer_capacity, uint32_t output_length[restrict static 1], double output_buffer[restrict static 1])
{
#ifdef __ERD_PROFILE__
    #ifdef _OPENMP
    const int tid = omp_get_thread_num();
    #else
    const int tid = 0;
    #endif
#endif
    ERD_PROFILE_START(erd__csgto)

/*             ...fix the A,B,C,D labels from the 1,2,3,4 ones. */
/*                Calculate the relevant data for the A,B,C,D batch of */
/*                integrals. */
    bool empty, tr1234;
    double xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd;
    uint32_t npgtoa, npgtob, npgtoc, npgtod;
    uint32_t shella, shellb, shellc, shelld;
    uint32_t nxyza, nxyzb, nxyzc, nxyzd;
    uint32_t nxyzet, nxyzft;
    uint32_t nrya, nryb, nryc, nryd;
    uint32_t nabcoor, ncdcoor;
    uint32_t ncolhrr, nrothrr, nxyzhrr;

    erd__set_abcd(npgto1, npgto2, npgto3, npgto4,
        shell1, shell2, shell3, shell4,
        atomic,
        x1, y1, z1, x2, y2, z2,
        x3, y3, z3, x4, y4, z4, spheric,
        &npgtoa, &npgtob, &npgtoc, &npgtod,
        &shella, &shellb, &shellc, &shelld,
        &xa, &ya, &za, &xb, &yb, &zb,
        &xc, &yc, &zc, &xd, &yd, &zd,
        &nxyza, &nxyzb, &nxyzc, &nxyzd,
        &nxyzet, &nxyzft,
        &nrya, &nryb, &nryc, &nryd,
        &nabcoor, &ncdcoor,
        &ncolhrr, &nrothrr, &nxyzhrr, &empty, &tr1234);
    if (empty) {
        *output_length = 0;
        ERD_PROFILE_END(erd__csgto)
        return;
    }

    uint32_t indexa = 0, indexb = 1, indexc = 2, indexd = 3;
    if (tr1234) {
        ERD_SWAP(alpha1, alpha3);
        ERD_SWAP(alpha2, alpha4);
        ERD_SWAP(cc1, cc3);
        ERD_SWAP(cc2, cc4);
        ERD_SWAP(norm1, norm3);
        ERD_SWAP(norm2, norm4);
        ERD_SWAP(indexa, indexc);
        ERD_SWAP(indexb, indexd);
        ERD_SWAP(shell1, shell3);
        ERD_SWAP(shell2, shell4);
    }
    
    const double *restrict alphaa = alpha1, *restrict alphab = alpha2, *restrict alphac = alpha3, *restrict alphad = alpha4;
    const double *restrict cca = cc1, *restrict ccb = cc2, *restrict ccc = cc3, *restrict ccd = cc4;
    const double *restrict norma = norm1, *restrict normb = norm2, *restrict normc = norm3, *restrict normd = norm4;
    if (shell1 < shell2) {
        ERD_SWAP(alphaa, alphab);
        ERD_SWAP(cca, ccb);
        ERD_SWAP(norma, normb);
        ERD_SWAP(indexa, indexb);
    }
    if (shell3 < shell4) {
        ERD_SWAP(alphac, alphad);
        ERD_SWAP(ccc, ccd);
        ERD_SWAP(normc, normd);
        ERD_SWAP(indexc, indexd);
    }
    
    // initialize values
    const double abx = xa - xb;
    const double aby = ya - yb;
    const double abz = za - zb;
    const double cdx = xc - xd;
    const double cdy = yc - yd;
    const double cdz = zc - zd;
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
    const uint32_t shellp = shella + shellb;
    const uint32_t shellq = shellc + shelld;
    const uint32_t shellt = shellp + shellq;
    const uint32_t mxshell = max4x32u(shella, shellb, shellc, shelld);
    const uint32_t nxyzp = (shellp + 1) * (shellp + 2) / 2;
    const uint32_t nxyzq = (shellq + 1) * (shellq + 2) / 2;
    // spnorm
    double spnorm = 1.0;
    if (shella == 1) {
        spnorm += spnorm;
    }
    if (shellb == 1) {
        spnorm += spnorm;
    }
    if (shellc == 1) {
        spnorm += spnorm;
    }
    if (shelld == 1) {
        spnorm += spnorm;
    }
    const double rnabsq = abx * abx + aby * aby + abz * abz;
    const double rncdsq = cdx * cdx + cdy * cdy + cdz * cdz;

   
/*             ...enter the cartesian contracted (e0|f0) batch */
/*                generation. Set the ij and kl primitive exponent */
/*                pairs and the corresponding exponential prefactors. */
    const uint32_t npgtoab = npgtoa * npgtob;
    const uint32_t npgtocd = npgtoc * npgtod;
    const uint32_t nxyzt = nxyzet * nxyzft;

    uint32_t prima[npgtoab], primb[npgtoab], primc[npgtocd], primd[npgtocd];
    ERD_SIMD_ALIGN double rhoab[PAD_LEN(npgtoab)];
    ERD_SIMD_ALIGN double rhocd[PAD_LEN(npgtocd)];
    uint32_t nij, nkl;
    ERD_PROFILE_START(erd__set_ij_kl_pairs)
    erd__set_ij_kl_pairs(npgtoa, npgtob, npgtoc, npgtod,
        xa, ya, za,
        xb, yb, zb,
        xc, yc, zc,
        xd, yd, zd,
        rnabsq, rncdsq, PREFACT,
        alphaa, alphab, alphac, alphad,
        &nij, &nkl,
        prima, primb, primc, primd,
        rhoab, rhocd);
    ERD_PROFILE_END(erd__set_ij_kl_pairs)

    if (nij * nkl == 0) {
        *output_length = 0;
        ERD_PROFILE_END(erd__csgto)
        return;
    }

/*             ...decide on the primitive [e0|f0] block size and */
/*                return array sizes and pointers for the primitive */
/*                [e0|f0] generation. Perform also some preparation */
/*                steps for contraction. */
    const uint32_t ngqp = shellt / 2 + 1;
    const uint32_t nmom = (ngqp << 1) - 1;
    const uint32_t ngqscr = nmom * 5 + (ngqp << 1) - 2;

    ERD_PROFILE_START(erd__prepare_ctr)
    const double factor = PREFACT * spnorm;
    const uint32_t npmin = min4x32u(npgtoa, npgtob, npgtoc, npgtod);
    double norm[npmin];
    if (npgtoa == npmin) {
        for (uint32_t i = 0; i < npgtoa; i++) {
            norm[i] = factor * norma[i];
        }
        norma = &norm[0];
    } else if (npgtob == npmin) {
        for (uint32_t i = 0; i < npgtob; i++) {
            norm[i] = factor * normb[i];
        }
        normb = &norm[0];
    } else if (npgtoc == npmin) {
        for (uint32_t i = 0; i < npgtoc; i++) {
            norm[i] = factor * normc[i];
        }
        normc = &norm[0];
    } else {
        for (uint32_t i = 0; i < npgtod; i++) {
            norm[i] = factor * normd[i];
        }
        normd = &norm[0];
    }
    ERD_PROFILE_END(erd__prepare_ctr)

/*             ...evaluate unnormalized rescaled [e0|f0] in blocks */
/*                over ij and kl pairs and add to final contracted */
/*                (e0|f0). The keyword REORDER indicates, if the */
/*                primitive [e0|f0] blocks need to be transposed */
/*                before being contracted. */
    ERD_PROFILE_START(erd__e0f0_pcgto_block)
    erd__e0f0_pcgto_block(nij, nkl, ngqp, nmom,
                           nxyzet, nxyzft, nxyzp, nxyzq,
                           shella, shellp, shellc, shellq,
                           xa, ya, za, xb, yb, zb,
                           xc, yc, zc, xd, yd, zd,
                           alphaa, alphab, alphac, alphad, 
                           cca, ccb, ccc, ccd,
                           vrrtab, ldvrrtab,
                           prima, primb, primc, primd,
                           norma, normb, normc, normd,
                           rhoab, rhocd,
                           output_buffer);
    ERD_PROFILE_END(erd__e0f0_pcgto_block)
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

    uint32_t in = 0;
    uint32_t out = nxyzhrr;

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
/*                 IHSCR = int output_buffer space for HRR matrix */
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
    const uint32_t nrowmx = (mxshell / 2 + 1) * (mxshell / 2 + 2) / 2;
    double fbuffer[spheric ? nrowmx * 4 : 0];
    uint32_t ibuffer[spheric ? nrowmx * 4 + nrya + nryb + nryc + nryd : 0];

    double cartnorm[spheric ? 0 : mxshell + 1];
    uint32_t isrowa, isrowb, isrowc, isrowd;
    uint32_t nrota, nrotb, nrotc, nrotd, nrowa, nrowb, nrowc, nrowd;
    uint32_t zsrota, zsrotb, zsrotc, zsrotd;
    uint32_t isnrowa, isnrowb, isnrowc, isnrowd;
    if (mxshell > 1) {
        if (spheric) {
            ERD_PROFILE_START(erd__xyz_to_ry_abcd)
            erd__xyz_to_ry_abcd(nxyza, nxyzb, nxyzc, nxyzd,
                                 nrya, nryb, nryc, nryd,
                                 shella, shellb, shellc, shelld,
                                 &nrowa, &nrowb, &nrowc, &nrowd,
                                 &nrota, &nrotb, &nrotc, &nrotd,
                                 &zsrota, &zsrotb, &zsrotc, &zsrotd,
                                 &isnrowa, &isnrowb, &isnrowc, &isnrowd,
                                 &isrowa, &isrowb, &isrowc, &isrowd,
                                 ibuffer, fbuffer);
            ERD_PROFILE_END(erd__xyz_to_ry_abcd)
        } else {
            erd__cartesian_norms(mxshell, cartnorm);
        }
    }

/*             ...do the first stage of processing the integrals: */
/*                   batch (ijkl,e0,f0) --> batch (ijkl,e0,cd) */
/*                   batch (ijkl,e0,c,d) --> batch (ijkl,e0,c,d') */
/*                   batch (ijkl,e0,c,d') --> batch (ijkl[d'],e0,c) */
/*                   batch (ijkl[d'],e0,c) --> batch (ijkl[d'],e0,c') */
/*                   batch (ijkl[d'],e0,c') --> batch (ijkl[c'd'],e0) */
    uint32_t ixoff[4] = { 1, 1, 1, 1 };
    if (shelld != 0) {
        uint32_t pos1, pos2, nrowhrr;
        ERD_SIMD_ALIGN double t[nrothrr*2];
        ERD_SIMD_ALIGN uint32_t row[nrothrr*2];
        ERD_SIMD_ALIGN uint32_t nrow[ncolhrr*2];
        ERD_PROFILE_START(erd__hrr_matrix)
        erd__hrr_matrix(nrothrr, ncolhrr, nxyzft, nxyzc, nxyzq,
                         shellc, shelld, shellq,
                         ncdcoor, cdx, cdy, cdz,
                         &pos1, &pos2, &nrowhrr,
                         nrow, row, t);
        ERD_PROFILE_END(erd__hrr_matrix)

        ERD_PROFILE_START(erd__hrr_transform)
        erd__hrr_transform(nxyzet, nrowhrr, nxyzc, nxyzd,
                            &nrow[pos1],
                            &row[pos2],
                            &t[pos2], &output_buffer[in],
                            &output_buffer[out]);
        ERD_PROFILE_END(erd__hrr_transform)

        ERD_SWAP(in, out);
        if (shelld > 1) {
            if (spheric) {
                ERD_PROFILE_START(erd__spherical_transform)
                erd__spherical_transform(nxyzet * nxyzc, nrowd, nryd,
                                          (uint32_t*)&ibuffer[isnrowd], (uint32_t*)&ibuffer[isrowd],
                                          &fbuffer[zsrotd], &output_buffer[in],
                                          &output_buffer[out]);
                ERD_PROFILE_END(erd__spherical_transform)

                ERD_SWAP(in, out);
            } else {
                erd__normalize_cartesian(nxyzet * nxyzc, shelld, cartnorm, &output_buffer[in]);
            }
        }
    }
    if (nryd > 1) {
        const uint32_t batch_size = nxyzet * nxyzc * nryd;
        const uint32_t notmove = ixoff[indexd];
        const uint32_t move = batch_size / (notmove * nryd);
        if (move > 1) {
            ERD_PROFILE_START(erd__move_ry)
            erd__move_ry(4, notmove, move, nryd, indexd, &output_buffer[in], ixoff, &output_buffer[out]);
            ERD_PROFILE_END(erd__move_ry)

            ERD_SWAP(in, out);
        }
    }
    if (shellc > 1) {
        if (spheric) {
            ERD_PROFILE_START(erd__spherical_transform)
            erd__spherical_transform(nxyzet * nryd, nrowc, nryc,
                                      (uint32_t*)&ibuffer[isnrowc], (uint32_t*)&ibuffer[isrowc],
                                      &fbuffer[zsrotc], &output_buffer[in],
                                      &output_buffer[out]);
            ERD_PROFILE_END(erd__spherical_transform)

            ERD_SWAP(in, out);
        } else {
            erd__normalize_cartesian(nxyzet * nryd, shellc, cartnorm, &output_buffer[in]);
        }
    }
    if (nryc > 1) {
        const uint32_t batch_size = nryd * nxyzet * nryc;
        const uint32_t notmove = ixoff[indexc];
        const uint32_t move = batch_size / (notmove * nryc);
        if (move > 1) {
            ERD_PROFILE_START(erd__move_ry)
            erd__move_ry(4, notmove, move, nryc, indexc, &output_buffer[in], ixoff, &output_buffer[out]);
            ERD_PROFILE_END(erd__move_ry)

            ERD_SWAP(in, out);
        }
    }

    /*
     * ...do the second stage of processing the integrals:
     *    batch (ijkl[c'd'],e0) --> batch (ijkl[c'd'],ab)
     *    batch (ijkl[c'd'],a,b) --> batch (ijkl[c'd'],a,b')
     *    batch (ijkl[c'd'],a,b') --> batch (ijkl[b'c'd'],a)
     *    batch (ijkl[b'c'd'],a) --> batch (ijkl[b'c'd'],a')
     *    batch (ijkl[b'c'd'],a') --> batch (ijkl[a'b'c'd'])
     */
    if (shellb != 0) {
        uint32_t pos1, pos2, nrowhrr;
        ERD_SIMD_ALIGN double t[nrothrr*2];
        ERD_SIMD_ALIGN uint32_t row[nrothrr*2];
        ERD_SIMD_ALIGN uint32_t nrow[ncolhrr*2];
        ERD_PROFILE_START(erd__hrr_matrix)
        erd__hrr_matrix(nrothrr, ncolhrr, nxyzet, nxyza, nxyzp,
                         shella, shellb, shellp,
                         nabcoor, abx, aby, abz,
                         &pos1, &pos2, &nrowhrr,
                         nrow, row, t);
        ERD_PROFILE_END(erd__hrr_matrix)

        ERD_PROFILE_START(erd__hrr_transform)
        erd__hrr_transform(nryc * nryd, nrowhrr, nxyza, nxyzb,
                            &nrow[pos1],
                            &row[pos2],
                            &t[pos2], &output_buffer[in],
                            &output_buffer[out]);
        ERD_PROFILE_END(erd__hrr_transform)

        ERD_SWAP(in, out);
        if (shellb > 1) {
            if (spheric) {
                ERD_PROFILE_START(erd__spherical_transform)
                erd__spherical_transform(nryc * nryd * nxyza, nrowb, nryb,
                                          (uint32_t*)&ibuffer[isnrowb], (uint32_t*)&ibuffer[isrowb],
                                          &fbuffer[zsrotb], &output_buffer[in], &output_buffer[out]);
                ERD_PROFILE_END(erd__spherical_transform)

                ERD_SWAP(in, out);
            } else {
                erd__normalize_cartesian(nryc * nryd * nxyza, shellb, cartnorm, &output_buffer[in]);
            }
        }
    }
    if (nryb > 1) {
        const uint32_t batch_size = nryc * nryd * nxyza * nryb;
        const uint32_t notmove = ixoff[indexb];
        const uint32_t move = batch_size / (notmove * nryb);
        if (move > 1) {
            ERD_PROFILE_START(erd__move_ry)
            erd__move_ry(4, notmove, move, nryb, indexb, &output_buffer[in], ixoff, &output_buffer[out]);
            ERD_PROFILE_END(erd__move_ry)

            ERD_SWAP(in, out);
        }
    }
    if (shella > 1) {
        if (spheric) {
            ERD_PROFILE_START(erd__spherical_transform)
            erd__spherical_transform(nryb * nryc * nryd, nrowa, nrya,
                                      &ibuffer[isnrowa], (uint32_t*)&ibuffer[isrowa],
                                      &fbuffer[zsrota], &output_buffer[in],
                                      &output_buffer[out]);
            ERD_PROFILE_END(erd__spherical_transform)

            ERD_SWAP(in, out);
        } else {
            erd__normalize_cartesian(nryb * nryc * nryd, shella, cartnorm, &output_buffer[in]);
        }
    }
    const uint32_t batch_size = nryb * nryc * nryd * nrya;
    if (nrya > 1) {
        const uint32_t notmove = ixoff[indexa];
        const uint32_t move = batch_size / (notmove * nrya);
        if (move > 1) {
            ERD_PROFILE_START(erd__move_ry)
            erd__move_ry(4, notmove, move, nrya, indexa, &output_buffer[in], ixoff, &output_buffer[out]);
            ERD_PROFILE_END(erd__move_ry)

            ERD_SWAP(in, out);
        }
    }

    if (in != 0) {
        memcpy(output_buffer, &output_buffer[in], batch_size * sizeof(double));
    }

    /* ...set final pointer to integrals in ZCORE array. */
    *output_length = batch_size;
    ERD_PROFILE_END(erd__csgto)
}
