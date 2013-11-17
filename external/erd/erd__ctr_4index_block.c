#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__CTR_4INDEX_BLOCK */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__TRANSPOSE_BATCH */
/*                ERD__CTR_1ST_HALF */
/*                ERD__MAP_IJKL_TO_IKJL */
/*                ERD__CTR_2ND_HALF_NEW */
/*                ERD__CTR_2ND_HALF_UPDATE */
/*  DESCRIPTION : This operation performs a four-indexed contraction */
/*                on a block of primitive integrals and updates all */
/*                contracted integrals: */
/*                  sum  [ab|cd]      -->  (ab|cd) */
/*                 ijkl         ijkl              rstu */
/*                                             r = 1,NCR */
/*                                             s = 1,NCS */
/*                                             t = 1,NCT */
/*                                             u = 1,NCU */
/*                For optimum performance of the contraction procedure, */
/*                the contraction routine is written to operate on the */
/*                block of primitive integrals ordered in the following */
/*                form (i.e. primitive indices to the far right): */
/*                              pbatch (nxyzt,kl,ij) */
/*                where ij and kl are the partial primitive index pairs */
/*                being treated here. If on entry to this routine the */
/*                primitive integrals are ordered in the form: */
/*                              pbatch (kl,ij,nxyzt) */
/*                a transposition has to preceed the actual contraction */
/*                procedure. This is triggered by the keyword PTRANS, */
/*                which has to be set true, if such a transposition is */
/*                needed. The complete set of cartesian monomial */
/*                quadruplets is kept together for efficiency through */
/*                the entire contraction sequence. */
/*                STRATEGY: */
/*                Perform overall contraction in 2 half contraction */
/*                steps, using an intermediate reordering if necessary: */
/*                1) Contract over ij primitives: */
/*                  (nxyzt,kl,rs) = ccr (r,i) * ccs (s,j) * (nxyzt,kl,ij) */
/*                2) Reorder kl <-> rs (if needed): */
/*                  (nxyzt,kl,rs) --> (nxyzt,rs,kl) */
/*                3) Contract over kl primitives and add result to the */
/*                   final contraction vector: */
/*                  (nxyzt,rs,tu) = (nxyzt,rs,tu) */
/*                                + cct (t,k) * ccu (u,l) * (nxyzt,rs,kl) */
/*                  Input: */
/*                    NxSIZE       =  size of the primitive integral */
/*                                    block (x=P), contracted integrals */
/*                                    (x=C) and working (x=W) arrays */
/*                    NXYZT        =  total # of cartesian monomial */
/*                                    quadruplets */
/*                    MIJKL        =  total # of partial primitive */
/*                                    index quadruplets */
/*                    MIJ(KL)      =  # of partial ij(kl) primitive */
/*                                    pairs to be transformed */
/*                    NRS(TU)      =  # of rs(tu) contraction pairs */
/*                    NPx          =  # of respective i,j,k,l primitives */
/*                                    for contractions x=R,S,T,U */
/*                    NCx          =  # of contractions for x=R,S,T,U */
/*                    MXPRIM       =  the maximum # of primitives */
/*                                    between all i,j,k,l primitives, */
/*                                    i.e. = max (i,j,k,l) */
/*                    MNPRIM       =  the minimum # of primitives */
/*                                    between i and j primitives and */
/*                                    k and l primitives and form the */
/*                                    maximum between these two values, */
/*                                    i.e. = max (min(i,j),min(k,l)) */
/*                    CCx          =  full set (including zeros) of */
/*                                    contraction coefficients for */
/*                                    x=R,S,T,U contractions */
/*                    CCBEGx       =  lowest nonzero primitive i,j,k,l */
/*                                    index for x=R,S,T,U contractions */
/*                    CCENDx       =  highest nonzero primitive i,j,k,l */
/*                                    index for x=R,S,T,U contractions */
/*                    PRIMx        =  primitive i,j,k,l indices for */
/*                                    the x=R,S,T,U contractions */
/*                    L1CACHE      =  Size of level 1 cache in units of */
/*                                    8 Byte */
/*                    TILE         =  Number of rows and columns in */
/*                                    units of 8 Byte of level 1 cache */
/*                                    square tile array used for */
/*                                    performing optimum matrix */
/*                                    transpositions */
/*                    NCTROW       =  minimum # of rows that are */
/*                                    accepted for blocked contractions */
/*                    EQUALRS(TU)  =  is true, if only the lower */
/*                                    triangle of ij(kl) primitive */
/*                                    indices is present and consequently */
/*                                    only the lower triangle of rs(tu) */
/*                                    contractions needs to be evaluated */
/*                    SWAPRS(TU)   =  if this is true, the 1st quarter */
/*                                    transformation will be over R(T) */
/*                                    followed by the 2nd over S(U). */
/*                                    If false, the order will be */
/*                                    reversed: 1st over S(U) then 2nd */
/*                                    over R(T) */
/*                    PTRANS       =  if true, a necessary primitive */
/*                                    integral transposition needs to */
/*                                    be done in order to bring the */
/*                                    ijkl primitive indices to the */
/*                                    far right position */
/*                    BLOCKED      =  if false, only one call will be */
/*                                    made to the present contraction */
/*                                    routine. The contraction batch */
/*                                    has not! been initialized to zero */
/*                                    and there is no need to perform */
/*                                    an update of the contracted batch. */
/*                    Pxxxx        =  intermediate storage arrays for */
/*                                    primitive labels to bundle */
/*                                    contraction steps in do loops */
/*                                    (xxxx = USED,SAVE,PAIR) */
/*                    PBATCH       =  the batch of primitive integrals */
/*                                    to be contracted */
/*                    WORK         =  the working array for intermediate */
/*                                    storage */
/*                    CBATCH       =  the batch of contracted integrals */
/*                                    before contraction update */
/*                  Output: */
/*                    CBATCH       =  the update batch of contracted */
/*                                    integrals after contraction */
/* ------------------------------------------------------------------------ */
int erd__ctr_4index_block (int npsize, int ncsize, int nwsize,
                           int nxyzt, int mijkl,
                           int mij, int mkl, int nrs, int ntu,
                           int npr, int nps, int npt, int npu,
                           int ncr, int ncs, int nct, int ncu,
                           int mxprim, int mnprim,
                           double *ccr, double *ccs,
                           double *cct, double *ccu,
                           int *ccbegr, int *ccbegs,
                           int *ccbegt, int *ccbegu,
                           int *ccendr, int *ccends,
                           int *ccendt, int *ccendu,
                           int *primr, int *prims,
                           int *primt, int *primu,
                           int l1cache, int tile, int nctrow,
                           int equalrs, int equaltu,
                           int swaprs, int swaptu,
                           int ptrans, int blocked,
                           int *pused, int *psave,
                           int *ppair, double *pbatch,
                           double *work, double *cbatch)
{
    int ccr_dim1, ccr_offset, ccs_dim1, ccs_offset, cct_dim1, cct_offset,
        ccu_dim1, ccu_offset;

    /* Local variables */
    int n;    
    int npmin;
    int npmax;
    int wused;
    int inwork;

    --pbatch;
    --cbatch;
    --work;
    --prims;
    --primr;
    --primu;
    --primt;
    --ccendr;
    --ccbegr;
    ccr_dim1 = npr;
    ccr_offset = 1 + ccr_dim1 * 1;
    ccr -= ccr_offset;
    --ccends;
    --ccbegs;
    ccs_dim1 = nps;
    ccs_offset = 1 + ccs_dim1 * 1;
    ccs -= ccs_offset;
    --ccendt;
    --ccbegt;
    cct_dim1 = npt;
    cct_offset = 1 + cct_dim1 * 1;
    cct -= cct_offset;
    --ccendu;
    --ccbegu;
    ccu_dim1 = npu;
    ccu_offset = 1 + ccu_dim1 * 1;
    ccu -= ccu_offset;
    --ppair;
    --psave;
    --pused;

    if (ptrans && mijkl > 1 && nxyzt > 1)
    {
        erd__transpose_batch (mijkl, nxyzt, &(pbatch[1]), &(work[1]));
        inwork = 1;
    }
    else
    {
        inwork = 0;
    }

/*             ...prepare for contraction over ij. Logical variable */
/*                INWORK controls where the actual significant data */
/*                (i.e. the integrals to be contracted) is residing. */
/*                If INWORK is true, then array WORK contains the */
/*                integrals, if false, they are in array PBATCH. */
/*                One of the key variables to be determined here */
/*                is the blocking size of the invariant indices N */
/*                such that the cache is used efficiently. The */
/*                blocking size has to be adapted also to the size */
/*                of the working space available. L1USED contains */
/*                the amount of data that will occupy the cache */
/*                besides the three main big arrays containing the */
/*                initial, quarter transformed and halftransformed */
/*                integrals. This extra data is the contraction */
/*                coefficients, their segmentation limits and the */
/*                primitive indices. L1FREE indicates the size of */
/*                the level 1 cache that is finally available for */
/*                the three big integral arrays. */
    npmax = MAX (npr, nps);
    npmin = MIN (npr, nps);

/*             ...do contraction over ij. */
    n = nxyzt * mkl;
    if (inwork)
    {
        wused = n * mij;
        erd__ctr_1st_half (n, npmax, npmin, mij, nrs, n, ncr, ncs,
                           npr, nps, &(ccr[ccr_offset]), &(ccs[ccs_offset]),
                           &ccbegr[1], &ccbegs[1], &ccendr[1], &(ccends[1]),
                           &primr[1], &prims[1], equalrs, swaprs, &pused[1],
                           &psave[1], &ppair[1], &work[1], &work[wused + 1],
                           &pbatch[1]);
        inwork = 0;
    }
    else
    {
        wused = n * nrs;
        erd__ctr_1st_half (n, npmax, npmin, mij, nrs, n, ncr, ncs,
                           npr, nps, &ccr[ccr_offset], &ccs[ccs_offset],
                           &ccbegr[1], &ccbegs[1], &ccendr[1], &ccends[1],
                           &primr[1], &prims[1], equalrs, swaprs, &pused[1],
                           &psave[1], &ppair[1], &pbatch[1],
                           &work[wused + 1], &work[1]);
        inwork = 1;
    }


/*             ...reorder rs <-> kl , if needed. */
    if (mkl > 1 && nrs > 1)
    {
        if (inwork)
        {
            erd__map_ijkl_to_ikjl (nxyzt, mkl, nrs, 1,
                                   &work[1], &pbatch[1]);
            inwork = 0;
        }
        else
        {
            erd__map_ijkl_to_ikjl (nxyzt, mkl, nrs, 1,
                                   &pbatch[1], &work[1]);
            inwork = 1;
        }
    }


/*             ...prepare for contraction over kl. Same procedure */
/*                as for ij (see comments above). The contracted */
/*                result will be added (updated) directly to the */
/*                final contracted CBATCH array. Do not perform a */
/*                contracted batch update, if no contraction blocking */
/*                is necessary. */

/*             ...do contraction over kl. */
    n = nxyzt * nrs;
    if (inwork)
    {
        wused = n * mkl;
        erd__ctr_2nd_half_new (n, npmax, npmin, mkl, ntu, n, nct,
                               ncu, npt, npu, &cct[cct_offset],
                               &ccu[ccu_offset], &ccbegt[1], &ccbegu[1],
                               &ccendt[1], &ccendu[1], &primt[1], &primu[1],
                               equaltu, swaptu, &pused[1], &psave[1],
                               &ppair[1], &work[1], &work[wused + 1],
                               &cbatch[1]);
    }
    else
    {
        wused = 0;
        erd__ctr_2nd_half_new (n, npmax, npmin, mkl, ntu, n, nct,
                               ncu, npt, npu, &cct[cct_offset],
                               &ccu[ccu_offset], &ccbegt[1], &ccbegu[1],
                               &ccendt[1], &ccendu[1], &primt[1], &primu[1],
                               equaltu, swaptu, &pused[1], &psave[1],
                               &ppair[1], &pbatch[1], &work[1], &cbatch[1]);
    }

    return 0;
}