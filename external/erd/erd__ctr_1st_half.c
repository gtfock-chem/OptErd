#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__CTR_1ST_HALF */
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
int erd__ctr_1st_half (int n, int mij,
                       int ncr, int ncs, int npr, int nps,
                       double *ccr, double *ccs,
                       int *ccbegr, int *ccbegs, int *ccendr,
                       int *ccends, int *primr, int *prims,
                       int equalrs, int swaprs, int *pused,
                       int *psave, int *ppair,
                       double *x, double *w, double *y)
{
    /* System generated locals */
    int ccr_offset;
    int ccs_offset;
    int w_offset;
    int y_offset;
    int x_offset;

    int i;
    int j;
    int l;
    int m;
    int r;
    int s;
    
    double c1;
    int i1;
    int j1;
    int ij;
    int ni;
    int nj;
    int rs;
    int ij1;
    int pmin;
    int pmax;
    int inext;
    int jnext;
    int irange;
    int jrange;
    int cntrct;;

    --ppair;
    --psave;
    --pused;
    x_offset = 1 + n * 1;
    x -= x_offset;
    --prims;
    --primr;
    y_offset = 1 + n * 1;
    y -= y_offset;
    w_offset = 1 + n * 1;
    w -= w_offset;
    --ccendr;
    --ccbegr;
    --ccends;
    --ccbegs;
    ccr_offset = 1 + npr * 1;
    ccr -= ccr_offset;
    ccs_offset = 1 + nps * 1;
    ccs -= ccs_offset;

    if (equalrs)
    {
/*             ...the R >= S case. Here we always have equal # of */
/*                primitives I and J for both R and S. The primitives */
/*                I >= J are ordered such that J varies fastest. */
/*                Outer contraction is over S, inner over R. */
         #if 0
        for (ij = 1; ij <= mij; ij++)
        {
            j = prims[ij];
            i = primr[ij];
            for (l = 1; l <= n; l++)
            {
                y[l + rs * n] = 0.0;
            }
            for (l = 1; l <= n; l++)
            {
                y[l + rs * n] +=
                        ccr[i + r * npr] * 
                        ccs[j + s * nps] *
                        x[l + ij * n];
            }
        }
       #else
        {
            for (i = 1; i <= npr; i++)
            {
                pused[i] = 0;
            }
            if (nps == 1 && ccs[1 + 1 * nps] == 1.0)
            {
/*             ...we reach this point, if only one S contraction */
/*                coefficient equal to 1.0 is present. */
                for (ij = 1; ij <= mij; ij++)
                {
                    j = prims[ij];
                    if (j == 1)
                    {
                        i = primr[ij];
                        for (l = 1; l <= n; l++)
                        {
                            w[l + i * n] = x[l + ij * n];
                        }
                        pused[i] = 1;
                    }
                }
            }
            else
            {
/*             ...the general S contraction case with many coefficients. */
/*                For each I determine all allowed J's and perform the */
/*                present S contraction for each I over all J's */
/*                simultaneously. */
                nj = 0;
                for (ij = 1; ij <= mij; ij++)
                {
                    i = primr[ij];
                    j = prims[ij];
                    jrange = ((j >= 1) && (j <= nps));
                    if (jrange)
                    {
                        nj++;
                        psave[nj] = j;
                        ppair[nj] = ij;
                    }
                    if (ij == mij)
                    {
                        inext = 0;
                    }
                    else
                    {
                        inext = primr[ij + 1];
                    }
                    cntrct = ((inext != i) && (nj > 0));
                    if (cntrct)
                    {
                        for (l = 1; l <= n; l++)
                        {
                            w[l + i * n] = 0.0;
                        }
                        for (m = 1; m <= nj; ++m)
                        {
                            ij1 = ppair[m];
                            c1 = ccs[psave[m] + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] += c1 * x[l + ij1 * n];
                            }
                        }
                        pused[i] = 1;                        
/*             ...don't forget to add into W the contributions due */
/*                to the J indices due to triangularization. Observe */
/*                that there is no possibility of grouping the J's */
/*                together, unless we resort the I and J such that */
/*                I varies fastest. Also observe that only J's */
/*                will be addressed that are =< I, hence we do */
/*                not run into the danger of addressing W columns */
/*                prematurely which have J > I, which saves us */
/*                from checking addressing of W columns during the */
/*                above simultaneous J contractions for each I. */
                        if (i >= 1 && i <= nps)
                        {
                            c1 = ccs[i + 1 * nps];
                            for (m = 1; m <= nj; ++m)
                            {
                                j = psave[m];
                                if (j != i)
                                {
                                    ij1 = ppair[m];
                                    if (pused[j] == 1)
                                    {
                                        for (l = 1; l <= n; l++)
                                        {
                                            w[l + j * n] += c1 * x[l + ij1 * n];
                                        }
                                    }
                                    else
                                    {
                                        for (l = 1; l <= n; l++)
                                        {
                                            w[l + j * n] = c1 * x[l + ij1 * n];
                                        }
                                        pused[j] = 1;
                                    }
                                }
                            }
                        }
                        nj = 0;
                    }
                }
            } /* if (pmin == pmax && ccs[pmin + s * nps] == 1.0) */
            
/*             ...inner contraction over all R >= S. */
        //    for (r = 1; r <= 1; r++)
            {
                if (1 == npr &&
                    ccr[1 + 1 * npr] == 1.0 &&
                    pused[1] == 1)
                {
/*             ...we reach this point, if only one R contraction */
/*                coefficient equal to 1.0 is present. */
                    for (l = 1; l <= n; l++)
                    {
                        y[l + 1 * n] = w[l + 1 * n];
                    }
                }
                else
                {
/*             ...the general R contraction case with many coefficients. */
/*                Group all relevant I's together and perform the */
/*                present R contraction over all I's simultaneously. */
                    ni = 0;
                    for (i = 1; i <= npr; ++i)
                    {
                        if (pused[i] == 1)
                        {
                            ni++;
                            psave[ni] = i;
                        }
                    }

                    for (l = 1; l <= n; l++)
                    {
                        y[l + 1 * n] = 0.0;
                    }
                    for (m = 1; m <= ni; ++m)
                    {
                        i1 = psave[m];
                        c1 = ccr[i1 + 1 * npr];
                        for (l = 1; l <= n; l++)
                        {
                            y[l + 1 * n] += c1 * w[l + i1 * n];
                        }
                    }
                }
            }
        }
       #endif
    }
    else /* if (!equalrs) */
    {
/*             ...the full RS case. Check the order of the two quarter */
/*                transformations and proceed accordingly. */
/*                The case: # of I primitives R > # of J primitives S */
/*                The primitives I and J are ordered such that I varies */
/*                fastest. Outer contraction is over R, inner over S. */
        if (swaprs)
        {
            rs = 0;
            for (r = 1; r <= ncr; ++r)
            {
                pmin = ccbegr[r];
                pmax = ccendr[r];
                for (j = 1; j <= nps; ++j)
                {
                    pused[j] = 0;
                }
                if (pmin == pmax &&
                    ccr[pmin + r * npr] == 1.0)
                {
/*             ...we reach this point, if only one R contraction */
/*                coefficient equal to 1.0 is present. */
                    for (ij = 1; ij <= mij; ++ij)
                    {
                        i = primr[ij];
                        if (i == pmin)
                        {
                            j = prims[ij];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + j * n] = x[l + ij * n];
                            }
                            pused[j] = 1;
                        }
                    }
                }
                else
                {
/*             ...the general R contraction case with many coefficients. */
/*                For each J determine all allowed I's and perform the */
/*                present R contraction for each J over all I's */
/*                simultaneously. */
                    ni = 0;
                    for (ij = 1; ij <= mij; ++ij)
                    {
                        i = primr[ij];
                        j = prims[ij];
                        irange = ((i >= pmin) && (i <= pmax));
                        if (irange)
                        {
                            ++ni;
                            psave[ni] = i;
                            ppair[ni] = ij;
                        }
                        if (ij == mij)
                        {
                            jnext = 0;
                        }
                        else
                        {
                            jnext = prims[ij + 1];
                        }
                        cntrct = ((jnext != j) && (ni > 0));
                        if (cntrct)
                        {
                            for (l = 1; l <= n; l++)
                            {
                                w[l + j * n] = 0.0;
                            }
                            for (m = 1; m <= ni; ++m)
                            {
                                ij1 = ppair[m];
                                c1 = ccr[psave[m] + r * npr];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + j * n] += c1 * x[l + ij1 * n];
                                }
                            }
                            pused[j] = 1;
                            ni = 0;
                        }
                    }
                }

/*             ...inner contraction over all S. */
                for (s = 1; s <= ncs; ++s)
                {
                    ++rs;
                    pmin = ccbegs[s];
                    pmax = ccends[s];
                    if (pmin == pmax &&
                        ccs[pmin + s * nps] == 1.0 &&
                        pused[pmin] == 1)
                    {
/*             ...we reach this point, if only one S contraction */
/*                coefficient equal to 1.0 is present. */
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = w[l + pmin * n];
                        }
                    }
                    else
                    {
/*             ...the general S contraction case with many coefficients. */
/*                Group all relevant J's together and perform the */
/*                present S contraction over all J's simultaneously. */
                        nj = 0;
                        for (j = pmin; j <= pmax; ++j)
                        {
                            if (pused[j] == 1)
                            {
                                ++nj;
                                psave[nj] = j;
                            }
                        }
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = 0.0;
                        }
                        for (m = 1; m <= nj; ++m)
                        {
                            j1 = psave[m];
                            c1 = ccs[j1 + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] += c1 * w[l + j1 * n];
                            }
                        }
                    }
                }
            }
        }
        else /* if (!swaprs) */
        {
/*             ...the case: # of I primitives R =< # of J primitives S. */
/*                The primitives I and J are ordered such that J varies */
/*                fastest. Outer contraction is over S, inner over R. */
            rs = 0;
            for (s = 1; s <= ncs; ++s)
            {
                pmin = ccbegs[s];
                pmax = ccends[s];
                for (i = 1; i <= npr; ++i)
                {
                    pused[i] = 0;
                }
                if (pmin == pmax &&
                    ccs[pmin + s * nps] == 1.0)
                {
                    for (ij = 1; ij <= mij; ++ij)
                    {
                        j = prims[ij];
                        if (j == pmin)
                        {
                            i = primr[ij];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] = x[l + ij * n];
                            }
                            pused[i] = 1;
                        }
                    }
                }
                else
                {
/*             ...the general S contraction case with many coefficients. */
/*                For each I determine all allowed J's and perform the */
/*                present S contraction for each I over all J's */
/*                simultaneously. */
                    nj = 0;
                    for (ij = 1; ij <= mij; ij++)
                    {
                        i = primr[ij];
                        j = prims[ij];
                        jrange = j >= pmin && j <= pmax;
                        if (jrange)
                        {
                            nj++;
                            psave[nj] = j;
                            ppair[nj] = ij;
                        }
                        if (ij == mij)
                        {
                            inext = 0;
                        }
                        else
                        {
                            inext = primr[ij + 1];
                        }
                        cntrct = ((inext != i) && (nj > 0));
                        if (cntrct)
                        {
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] = 0.0;
                            }                        
                            for (m = 1; m <= nj; ++m)
                            {
                                ij1 = ppair[m];
                                c1 = ccs[psave[m] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] += c1 * x[l + ij1 * n];
                                }
                            }
                            pused[i] = 1;
                            nj = 0;
                        }
                    }
                }

/*             ...inner contraction over all R. */
                for (r = 1; r <= ncr; ++r)
                {
                    rs++;
                    pmin = ccbegr[r];
                    pmax = ccendr[r];
                    if (pmin == pmax &&
                        ccr[pmin + r * npr] == 1.0 &&
                        pused[pmin] == 1)
                    {
/*             ...we reach this point, if only one R contraction */
/*                coefficient equal to 1.0 is present. */
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = w[l + pmin * n];
                        }
                    }
                    else
                    {
/*             ...the general R contraction case with many coefficients. */
/*                Group all relevant I's together and perform the */
/*                present R contraction over all I's simultaneously. */
                        ni = 0;
                        for (i = pmin; i <= pmax; i++)
                        {
                            if (pused[i] == 1)
                            {
                                ni++;
                                psave[ni] = i;
                            }
                        }
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = 0.0;
                        }
                        for (m = 1; m <= ni; m++)
                        {
                            i1 = psave[m];
                            c1 = ccr[i1 + r * npr];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] += c1 * w[l + i1 * n];
                            }
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}