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
int erd__ctr_4index_block (int nxyzt, int mijkl,
                           int mij, int mkl, int nrs, int ntu,
                           int npr, int nps, int npt, int npu,
                           int ncr, int ncs, int nct, int ncu,
                           double *ccr, double *ccs,
                           double *cct, double *ccu,
                           int *ccbegr, int *ccbegs,
                           int *ccbegt, int *ccbegu,
                           int *ccendr, int *ccends,
                           int *ccendt, int *ccendu,
                           int *primr, int *prims,
                           int *primt, int *primu,
                           int equalrs, int equaltu,
                           int swaprs, int swaptu, int ptrans,
                           int *pused, int *psave,
                           int *ppair, double *pbatch,
                           double *work, double *cbatch)
{
    int n;    
    int npmin;
    int npmax;
    int wused;
    int inwork;

    if (ptrans && mijkl > 1 && nxyzt > 1)
    {
        erd__transpose_batch (mijkl, nxyzt, pbatch, work);
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
                           npr, nps, ccr, ccs,
                           ccbegr, ccbegs, ccendr, ccends,
                           primr, prims, equalrs, swaprs,
                           pused, psave, ppair,
                           work, &(work[wused]), pbatch);
        inwork = 0;
    }
    else
    {
        wused = n * nrs;
        erd__ctr_1st_half (n, npmax, npmin, mij, nrs, n, ncr, ncs,
                           npr, nps, ccr, ccs,
                           ccbegr, ccbegs, ccendr, ccends,
                           primr, prims, equalrs, swaprs,
                           pused, psave, ppair,
                           pbatch, &(work[wused]), work);
        inwork = 1;
    }


/*             ...reorder rs <-> kl , if needed. */
    if (mkl > 1 && nrs > 1)
    {
        if (inwork)
        {
            erd__map_ijkl_to_ikjl (nxyzt, mkl, nrs, 1,
                                   work, pbatch);
            inwork = 0;
        }
        else
        {
            erd__map_ijkl_to_ikjl (nxyzt, mkl, nrs, 1,
                                   pbatch, work);
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
                               ncu, npt, npu, cct, ccu,
                               ccbegt, ccbegu, ccendt, ccendu,
                               primt, primu, equaltu, swaptu,
                               pused, psave, ppair,
                               work, &(work[wused]), cbatch);
    }
    else
    {
        wused = 0;
        erd__ctr_2nd_half_new (n, npmax, npmin, mkl, ntu, n, nct,
                               ncu, npt, npu, cct, ccu,
                               ccbegt, ccbegu, ccendt, ccendu,
                               primt, primu, equaltu, swaptu,
                               pused, psave, ppair,
                               pbatch, work, cbatch);
    }

    return 0;
}