#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


#pragma offload_attribute(push, target(mic))

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__CTR_HALF */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation performs the first half contraction */
/*                step on the incomming integrals over primitives */
/*                in blocked form over invariant indices n: */
/*                    y (n,rs) = sum  ccr (r,i) * ccs (s,j) * x (n,ij) */
/*                                ij */
/*                where ccr and ccs are the arrays containing the */
/*                contraction coefficients. The sum is over the i and */
/*                j primitives which are transmitted in arrays PRIMR */
/*                and PRIMS, respectively, and may constitute only */
/*                a subset of the full range of primitives. */
/*                The contraction is split into two quarter steps, */
/*                the order of which is determined by the # of i and */
/*                j primitives: */
/*                   a) w (n,j/i) = sum ccr/s (r/s,i/j) * x (n,ij) */
/*                                  i/j */
/*                   b) y (n,rs)  = sum ccs/r (s/r,j/i) * w (n,j/i) */
/*                                  j/i */
/*                Size of the intermediate w (n,i) or w (n,j) array */
/*                has to be kept to a minimum and blocking over the */
/*                invariant indices n has to be performed such that */
/*                the intermediate w array does not get kicked out */
/*                from the cache lines after the first quarter */
/*                transformation. */
/*                In case of csh equality (EQUALRS = .true), we only */
/*                have to consider the lower triangle of the primitive */
/*                integrals, which, with the exception of the diagonals, */
/*                have to be used twice. */
/*                        --- SEGMENTED CONTRACTIONS --- */
/*                Segmented contractions are those defined to be */
/*                within a certain consecutive i- and j-range of */
/*                primitives. The segmented limits for each contraction */
/*                are sitting respectively in CCBEGR and CCBEGS (lowest */
/*                limit) and CCENDR and CCENDS (highest limit) and they */
/*                determine which of the actual i's and j's from the */
/*                PRIMR and PRIMS lists have to be considered for each */
/*                contraction index. */
/*                The code also allows efficient contractions in case */
/*                there is only one contraction coefficient present */
/*                in a certain contraction and its value is equal to 1. */
/*                In such cases we can save lots of multiplications by */
/*                1 and since these cases are quite common for certain */
/*                types of basis functions it is worth including some */
/*                IF's inside the contraction loops to gain speed. */
/*                  Input: */
/*                    N            =  # of invariant indices */
/*                    NPMAX(MIN)   =  the maximum (minimum) # of */
/*                                    primitives between both primitive */
/*                                    sets i,j */
/*                    MIJ          =  # of ij primitive products to */
/*                                    be transformed */
/*                    NRS          =  # of rs contractions to be done */
/*                    NBLOCK       =  blocking size for invariant */
/*                                    indices n */
/*                    NCR(S)       =  # of contractions for the i -> R */
/*                                    (j -> S) primitives */
/*                    NPR(S)       =  # of i(j) primitives */
/*                    CCR(S)       =  full set (including zeros) of */
/*                                    contraction coefficients for */
/*                                    R(S) contractions */
/*                    CCBEGR(S)    =  lowest nonzero primitive i(j) */
/*                                    index for R(S) contractions */
/*                    CCENDR(S)    =  highest nonzero primitive i(j) */
/*                                    index for R(S) contractions */
/*                    PRIMR(S)     =  primitive i(j) indices */
/*                    EQUALRS      =  is true, if only the lower */
/*                                    triangle of ij primitive indices */
/*                                    is present and consequently */
/*                                    only the lower triangle of rs */
/*                                    contractions needs to be evaluated */
/*                    SWAPRS       =  if this is true, the 1st quarter */
/*                                    transformation is over R followed */
/*                                    by the 2nd over S. If false, the */
/*                                    order is reversed: 1st over S then */
/*                                    2nd over R */
/*                    Pxxxx        =  intermediate storage arrays for */
/*                                    primitive labels to bundle */
/*                                    contraction steps in do loops */
/*                                    (xxxx = USED,SAVE,PAIR) */
/*                    X            =  array containing the primitive */
/*                                    integrals */
/*                    W            =  intermediate storage array */
/*                                    containing 1st quarter transformed */
/*                                    primitive integrals */
/*                  Output: */
/*                    Y            =  contains the final half transformed */
/*                                    integrals */
/* ------------------------------------------------------------------------ */
static int
erd__ctr_half (int n, int mij, double *ccr, double *ccs,
               int *primr, int *prims, double *x, double *y)
{
    int i;
    int j;
    int l;
    double c1;
    double c2;
    int ij;

    for (l = 0; l < n; l++)
    {
        y[l] = 0.0;
    }

/*             ...the full RS case. Check the order of the two quarter */
/*                transformations and proceed accordingly. */
/*                The case: # of I primitives R > # of J primitives S */
/*                The primitives I and J are ordered such that I varies */
/*                fastest. Outer contraction is over R, inner over S. */
    for (ij = 0; ij < mij; ij++)
    {
        j = prims[ij];
        i = primr[ij];
        c1 = ccr[i];
        c2 = ccs[j];
        for (l = 0; l < n; l++)
        {
            y[l] += c1 * c2 * x[l + ij * n];
        }
    }

    return 0;
}


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
int
erd__ctr_4index_block (int nxyzt, int mij, int mkl,
                       double *ccr, double *ccs,
                       double *cct, double *ccu,
                       int *primr, int *prims,
                       int *primt, int *primu,
                       double *pbatch, double *work, double *cbatch)
{
    int n;
    int mijkl;

    mijkl = mij * mkl;
    n = nxyzt * mkl;

    erd__transpose_batch (mijkl, nxyzt, pbatch, work);
    erd__ctr_half (n, mij, ccr, ccs, primr, prims, work, pbatch);
    erd__ctr_half (nxyzt, mkl, cct, ccu, primt, primu, pbatch, cbatch);

    return 0;
}

#pragma offload_attribute(pop)
