#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__CTR_2ND_HALF_NEW */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation performs the second half contraction */
/*                step on the incomming half transformed integrals */
/*                in blocked form over invariant indices n and generates */
/*                new fully transformed contracted batch of integrals: */
/*                   y (n,tu) = sum  cct (t,k) * ccu (u,l) * x (n,kl) */
/*                               kl */
/*                where cct and ccu are the arrays containing the */
/*                contraction coefficients. The sum is over the k and */
/*                l primitives which are transmitted in arrays PRIMT */
/*                and PRIMU, respectively, and may constitute only */
/*                a subset of the full range of primitives. */
/*                The contraction is split into two quarter steps, */
/*                the order of which is determined by the # of k and */
/*                l primitives: */
/*                   a) w (n,l/k) = sum cct/u (t/u,k/l) * x (n,kl) */
/*                                  k/l */
/*                   b) y (n,tu)  = sum ccu/t (u/t,l/k) * w (n,l/k) */
/*                                  l/k */
/*                Size of the intermediate w (n,k) or w (n,l) array */
/*                has to be kept to a minimum and blocking over the */
/*                invariant indices n has to be performed such that */
/*                the intermediate w array does not get kicked out */
/*                from the cache lines after the first quarter */
/*                transformation. */
/*                In case of csh equality (EQUALTU = .true), we only */
/*                have to consider the lower triangle of the primitive */
/*                integrals, which, with the exception of the diagonals, */
/*                have to be used twice. */
/*                        --- SEGMENTED CONTRACTIONS --- */
/*                Segmented contractions are those defined to be */
/*                within a certain consecutive k- and l-range of */
/*                primitives. The segmented limits for each contraction */
/*                are sitting respectively in CCBEGT and CCBEGU (lowest */
/*                limit) and CCENDT and CCENDU (highest limit) and they */
/*                determine which of the actual k's and l's from the */
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
/*                                    sets k,l */
/*                    MKL          =  # of kl primitive products to */
/*                                    be transformed */
/*                    NTU          =  # of tu contractions to be done */
/*                    NBLOCK       =  blocking size for invariant */
/*                                    indices n */
/*                    NCT(U)       =  # of contractions for the k -> T */
/*                                    (l -> U) primitives */
/*                    NPT(U)       =  # of k(l) primitives */
/*                    CCT(U)       =  full set (including zeros) of */
/*                                    contraction coefficients for */
/*                                    T(U) contractions */
/*                    CCBEGT(U)    =  lowest nonzero primitive k(l) */
/*                                    index for T(U) contractions */
/*                    CCENDT(U)    =  highest nonzero primitive k(l) */
/*                                    index for T(U) contractions */
/*                    PRIMT(U)     =  primitive k(l) indices */
/*                    EQUALTU      =  is true, if only the lower */
/*                                    triangle of kl primitive indices */
/*                                    is present and consequently */
/*                                    only the lower triangle of tu */
/*                                    contractions needs to be evaluated */
/*                    SWAPTU       =  if this is true, the 1st quarter */
/*                                    transformation is over T followed */
/*                                    by the 2nd over U. If false, the */
/*                                    order is reversed: 1st over U then */
/*                                    2nd over T */
/*                    Pxxxx        =  intermediate storage arrays for */
/*                                    primitive labels to bundle */
/*                                    contraction steps in do loops */
/*                                    (xxxx = USED,SAVE,PAIR) */
/*                    X            =  array containing the half */
/*                                    transformed integrals */
/*                    W            =  intermediate storage array */
/*                                    containing 3/4 transformed */
/*                                    integrals */
/*                    Y            =  original batch of fully contracted */
/*                                    integrals */
/*                  Output: */
/*                    Y            =  updated batch of fully contracted */
/*                                    integrals */
/* ------------------------------------------------------------------------ */
int erd__ctr_2nd_half_new (int n, int npmax, int npmin,
                           int mkl, int ntu, int nblock,
                           int nct, int ncu, int npt, int npu,
                           double *cct, double *ccu,
                           int *ccbegt, int *ccbegu,
                           int *ccendt, int *ccendu,
                           int *primt, int *primu,
                           int equaltu, int swaptu,
                           int *pused, int *psave, int *ppair,
                           double *x, double *w, double *y)
{
    /* System generated locals */
    int cct_offset;
    int ccu_offset;
    int w_offset;
    int x_offset;
    int y_offset;

    int i;
    int j;
    int k;
    int l;
    int m;
    int t;
    int u;
    double c1;
    double c2;
    double c3;
    double c4;
    double c5;
    double c6;
    double c7;
    double c8;
    int k1;
    int k2;
    int k3;
    int k4;
    int k5;
    int k6;
    int k7;
    int k8;
    int l1;
    int l2;
    int l3;
    int l4;
    int l5;
    int l6;
    int l7;
    int l8;
    
    int kl;
    int nk;
    int nl;
    int tu;
    int kl1;
    int kl2;
    int kl3;
    int kl4;
    int kl5;
    int kl6;
    int kl7;
    int kl8;
    int pmin;
    int pmax;
    int knext;
    int lnext;
    int nkbase;
    int nlbase;
    int krange;
    int lrange;
    int nkleft;
    int nlleft;
    int cntrct;
    int nkstep;
    int nlstep;
    int nkrest;
    int nlrest;

    /* Parameter adjustments */
    --ppair;
    --psave;
    --pused;
    x_offset = 1 + n * 1;
    x -= x_offset;
    --primu;
    --primt;
    y_offset = 1 + n * 1;
    y -= y_offset;
    w_offset = 1 + n * 1;
    w -= w_offset;
    --ccendt;
    --ccbegt;
    --ccendu;
    --ccbegu;
    cct_offset = 1 + npt * 1;
    cct -= cct_offset;
    ccu_offset = 1 + npu * 1;
    ccu -= ccu_offset;

    if (equaltu)
    {
/*             ...the T >= U case. Here we always have equal # of */
/*                primitives K and L for both T and U. The primitives */
/*                K >= L are ordered such that L varies fastest. */
/*                Outer contraction is over U, inner over T. */
        tu = 0;
        for (u = 1; u <= ncu; ++u)
        {
            pmin = ccbegu[u];
            pmax = ccendu[u];
            for (k = 1; k <= npt; ++k)
            {
                pused[k] = 0;
            }
            if (pmin == pmax &&
                ccu[pmin + u * npu] == 1.0)
            {
/*             ...we reach this point, if only one U contraction */
/*                coefficient equal to 1.0 is present. */
                for (kl = 1; kl <= mkl; ++kl)
                {
                    l = primu[kl];
                    if (l == pmin)
                    {
                        k = primt[kl];
                        for (j = 1; j <= n; ++j)
                        {
                            w[j + k * n] = x[j + kl * n];
                        }
                        pused[k] = 1;
                    }
                }
            }
            else
            {
/*             ...the general U contraction case with many coefficients. */
/*                For each K determine all allowed L's and perform the */
/*                present U contraction for each K over all L's */
/*                simultaneously. */
                nl = 0;
                for (kl = 1; kl <= mkl; ++kl)
                {
                    k = primt[kl];
                    l = primu[kl];
                    lrange = ((l >= pmin) && (l <= pmax));
                    if (lrange)
                    {
                        ++nl;
                        psave[nl] = l;
                        ppair[nl] = kl;
                    }
                    if (kl == mkl)
                    {
                        knext = 0;
                    }
                    else
                    {
                        knext = primt[kl + 1];
                    }
                    cntrct = knext != k && nl > 0;
                    if (cntrct)
                    {
                        for (j = 1; j <= n; ++j)
                        {
                            w[j + k * n] = 0.0;
                        }                        
                        for (m = 1; m <= nl; ++m)
                        {
                            kl1 = ppair[m];
                            c1 = ccu[psave[m] + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] += c1 * x[j + kl1 * n];
                            }
                        }
                        pused[k] = 1;
/*             ...don't forget to add into W the contributions due */
/*                to the L indices due to triangularization. Observe */
/*                that there is no possibility of grouping the L's */
/*                together, unless we resort the K and L such that */
/*                K varies fastest. Also observe that only L's */
/*                will be addressed that are =< K, hence we do */
/*                not run into the danger of addressing W columns */
/*                prematurely which have L > K, which saves us */
/*                from checking addressing of W columns during the */
/*                above simultaneous L contractions for each K. */
                        if (k >= pmin && k <= pmax)
                        {
                            c1 = ccu[k + u * npu];
                            for (m = 1; m <= nl; ++m)
                            {
                                l = psave[m];
                                if (l != k)
                                {
                                    kl1 = ppair[m];
                                    if (pused[l] == 1)
                                    {
                                        for (j = 1; j <= n; ++j)
                                        {
                                            w[j + l * n] += c1 * x[j +
                                                                        kl1 *
                                                                        n];
                                        }
                                    }
                                    else
                                    {
                                        for (j = 1; j <= n; ++j)
                                        {
                                            w[j + l * n] = c1 * x[j +
                                                                       kl1 *
                                                                       n];
                                        }
                                        pused[l] = 1;
                                    }
                                }
                            }
                        }
                        nl = 0;
                    }
                }
            }

/*             ...inner contraction over all T >= U. */
            for (t = u; t <= nct; ++t)
            {
                ++tu;
                pmin = ccbegt[t];
                pmax = ccendt[t];
                if (pmin == pmax &&
                    cct[pmin + t * npt] == 1.0
                    && pused[pmin] == 1)
                {
/*             ...we reach this point, if only one T contraction */
/*                coefficient equal to 1.0 is present. */
                    for (j = 1; j <= n; ++j)
                    {
                        y[j + tu * n] = w[j + pmin * n];
                    }
                }
                else
                {
/*             ...the general T contraction case with many coefficients. */
/*                Group all relevant K's together and perform the */
/*                present T contraction over all K's simultaneously. */
                    nk = 0;
                    for (k = pmin; k <= pmax; ++k)
                    {
                        if (pused[k] == 1)
                        {
                            ++nk;
                            psave[nk] = k;
                        }
                    }
                    for (j = 1; j <= n; ++j)
                    {
                        y[j + tu * n] = 0.0;
                    }
                    for (m = 1; m <= nk; ++m)
                    {
                        k1 = psave[m];
                        c1 = cct[k1 + t * npt];
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] += c1 * w[j + k1 * n];
                        }
                    }
                }
            }
        }
    }
    else
    {
/*             ...the full TU case. Check the order of the two quarter */
/*                transformations and proceed accordingly. */
/*                The case: # of K primitives T > # of L primitives U */
/*                The primitives K and L are ordered such that K varies */
/*                fastest. Outer contraction is over T, inner over U. */
        if (swaptu)
        {
            tu = 0;
            for (t = 1; t <= nct; ++t)
            {
                pmin = ccbegt[t];
                pmax = ccendt[t];
                for (l = 1; l <= npu; ++l)
                {
                    pused[l] = 0;
                }
                if (pmin == pmax && cct[pmin + t * npt] == 1.)
                {
/*             ...we reach this point, if only one T contraction */
/*                coefficient equal to 1.0 is present. */
                    for (kl = 1; kl <= mkl; ++kl)
                    {
                        k = primt[kl];
                        if (k == pmin)
                        {
                            l = primu[kl];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + l * n] = x[j + kl * n];
                            }
                            pused[l] = 1;
                        }
                    }
                }
                else
                {
/*             ...the general T contraction case with many coefficients. */
/*                For each L determine all allowed K's and perform the */
/*                present T contraction for each L over all K's */
/*                simultaneously. */
                    nk = 0;
                    for (kl = 1; kl <= mkl; ++kl)
                    {
                        k = primt[kl];
                        l = primu[kl];
                        krange = ((k >= pmin) && (k <= pmax));
                        if (krange)
                        {
                            ++nk;
                            psave[nk] = k;
                            ppair[nk] = kl;
                        }
                        if (kl == mkl)
                        {
                            lnext = 0;
                        }
                        else
                        {
                            lnext = primu[kl + 1];
                        }
                        cntrct = lnext != l && nk > 0;
                        if (cntrct)
                        {
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + l * n] = 0.0;
                            }
                            for (m = 1; m <= nk; ++m)
                            {
                                kl1 = ppair[m];
                                c1 = cct[psave[m] + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + l * n] += c1 * x[j + kl1 * n];
                                }
                            }
                            pused[l] = 1;
                            nk = 0;
                        }
                    }
                }
/*             ...inner contraction over all U. */
                for (u = 1; u <= ncu; ++u)
                {
                    ++tu;
                    pmin = ccbegu[u];
                    pmax = ccendu[u];
                    if (pmin == pmax &&
                        ccu[pmin + u * npu] == 1.0 &&
                        pused[pmin] == 1)
                    {
/*             ...we reach this point, if only one U contraction */
/*                coefficient equal to 1.0 is present. */
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = w[j + pmin * n];
                        }
                    }
                    else
                    {
/*             ...the general U contraction case with many coefficients. */
/*                Group all relevant L's together and perform the */
/*                present U contraction over all L's simultaneously. */
                        nl = 0;
                        for (l = pmin; l <= pmax; ++l)
                        {
                            if (pused[l] == 1)
                            {
                                ++nl;
                                psave[nl] = l;
                            }
                        }
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = 0.0;
                        }
                        for (m = 1; m <= nl; ++m)
                        {
                            l1 = psave[m];
                            c1 = ccu[l1 + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] += c1 * w[j + l1 * n];
                            }
                        }
                    }
                }
            }
        }
        else
        {

/*             ...the case: # of K primitives T =< # of L primitives U. */
/*                The primitives K and L are ordered such that L varies */
/*                fastest. Outer contraction is over U, inner over T. */
            tu = 0;
            for (u = 1; u <= ncu; ++u)
            {
                pmin = ccbegu[u];
                pmax = ccendu[u];
                for (k = 1; k <= npt; ++k)
                {
                    pused[k] = 0;
                }
                if (pmin == pmax && ccu[pmin + u * npu] == 1.)
                {
/*             ...we reach this point, if only one U contraction */
/*                coefficient equal to 1.0 is present. */
                    for (kl = 1; kl <= mkl; ++kl)
                    {
                        l = primu[kl];
                        if (l == pmin)
                        {
                            k = primt[kl];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] = x[j + kl * n];
                            }
                            pused[k] = 1;
                        }
                    }
                }
                else
                {
/*             ...the general U contraction case with many coefficients. */
/*                For each K determine all allowed L's and perform the */
/*                present U contraction for each K over all L's */
/*                simultaneously. */
                    nl = 0;
                    for (kl = 1; kl <= mkl; ++kl)
                    {
                        k = primt[kl];
                        l = primu[kl];
                        lrange = l >= pmin && l <= pmax;
                        if (lrange)
                        {
                            ++nl;
                            psave[nl] = l;
                            ppair[nl] = kl;
                        }
                        if (kl == mkl)
                        {
                            knext = 0;
                        }
                        else
                        {
                            knext = primt[kl + 1];
                        }
                        cntrct = knext != k && nl > 0;
                        if (cntrct)
                        {
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] = 0.0;
                            }
                            for (m = 1; m <= nl; ++m)
                            {
                                kl1 = ppair[m];
                                c1 = ccu[psave[m] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] += c1 * x[j + kl1 * n];
                                }
                            }
                            pused[k] = 1;
                            nl = 0;
                        }
                    }
                }
/*             ...inner contraction over all T. */
                for (t = 1; t <= nct; ++t)
                {
                    ++tu;
                    pmin = ccbegt[t];
                    pmax = ccendt[t];
                    if (pmin == pmax && cct[pmin + t * npt] == 1. &&
                        pused[pmin] == 1)
                    {

/*             ...we reach this point, if only one T contraction */
/*                coefficient equal to 1.0 is present. */
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = w[j + pmin * n];
                        }
                    }
                    else
                    {
/*             ...the general T contraction case with many coefficients. */
/*                Group all relevant K's together and perform the */
/*                present T contraction over all K's simultaneously. */
                        nk = 0;
                        for (k = pmin; k <= pmax; ++k)
                        {
                            if (pused[k] == 1)
                            {
                                ++nk;
                                psave[nk] = k;
                            }
                        }
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = 0.0;
                        }
                        for (m = 1; m <= nk; ++m)
                        {
                            k1 = psave[m];
                            c1 = cct[k1 + t * npt];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] += c1 * w[j + k1 * n];
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}
