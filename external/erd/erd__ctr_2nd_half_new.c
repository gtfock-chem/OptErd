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
                        if (nl == 1)
                        {
                            kl1 = ppair[1];
                            c1 = ccu[psave[1] + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] = c1 * x[j + kl1 * n];
                            }
                        }
                        else if (nl == 2)
                        {
                            kl1 = ppair[1];
                            kl2 = ppair[2];
                            c1 = ccu[psave[1] + u * npu];
                            c2 = ccu[psave[2] + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] = c1 * x[j + kl1 * n] +
                                    c2 * x[j + kl2 * n];
                            }
                        }
                        else if (nl == 3)
                        {
                            kl1 = ppair[1];
                            kl2 = ppair[2];
                            kl3 = ppair[3];
                            c1 = ccu[psave[1] + u * npu];
                            c2 = ccu[psave[2] + u * npu];
                            c3 = ccu[psave[3] + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] = c1 * x[j + kl1 * n] +
                                    c2 * x[j + kl2 * n] +
                                    c3 * x[j + kl3 * n];
                            }
                        }
                        else if (nl == 4)
                        {
                            kl1 = ppair[1];
                            kl2 = ppair[2];
                            kl3 = ppair[3];
                            kl4 = ppair[4];
                            c1 = ccu[psave[1] + u * npu];
                            c2 = ccu[psave[2] + u * npu];
                            c3 = ccu[psave[3] + u * npu];
                            c4 = ccu[psave[4] + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] = c1 * x[j + kl1 *
                                                           n] +
                                    c2 * x[j + kl2 * n] + c3 * x[j +
                                                                      kl3 *
                                                                      n]
                                    + c4 * x[j + kl4 * n];
                            }
                        }
                        else if (nl == 5)
                        {
                            kl1 = ppair[1];
                            kl2 = ppair[2];
                            kl3 = ppair[3];
                            kl4 = ppair[4];
                            kl5 = ppair[5];
                            c1 = ccu[psave[1] + u * npu];
                            c2 = ccu[psave[2] + u * npu];
                            c3 = ccu[psave[3] + u * npu];
                            c4 = ccu[psave[4] + u * npu];
                            c5 = ccu[psave[5] + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] = c1 * x[j + kl1 *
                                                           n] +
                                    c2 * x[j + kl2 * n] + c3 * x[j +
                                                                      kl3 *
                                                                      n]
                                    + c4 * x[j + kl4 * n] + c5 * x[j +
                                                                        kl5 *
                                                                        n];
                            }
                        }
                        else if (nl == 6)
                        {
                            kl1 = ppair[1];
                            kl2 = ppair[2];
                            kl3 = ppair[3];
                            kl4 = ppair[4];
                            kl5 = ppair[5];
                            kl6 = ppair[6];
                            c1 = ccu[psave[1] + u * npu];
                            c2 = ccu[psave[2] + u * npu];
                            c3 = ccu[psave[3] + u * npu];
                            c4 = ccu[psave[4] + u * npu];
                            c5 = ccu[psave[5] + u * npu];
                            c6 = ccu[psave[6] + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] = c1 * x[j + kl1 *
                                                           n] +
                                    c2 * x[j + kl2 * n] + c3 * x[j +
                                                                      kl3 *
                                                                      n]
                                    + c4 * x[j + kl4 * n] + c5 * x[j +
                                                                        kl5 *
                                                                        n]
                                    + c6 * x[j + kl6 * n];
                            }
                        }
                        else if (nl == 7)
                        {
                            kl1 = ppair[1];
                            kl2 = ppair[2];
                            kl3 = ppair[3];
                            kl4 = ppair[4];
                            kl5 = ppair[5];
                            kl6 = ppair[6];
                            kl7 = ppair[7];
                            c1 = ccu[psave[1] + u * npu];
                            c2 = ccu[psave[2] + u * npu];
                            c3 = ccu[psave[3] + u * npu];
                            c4 = ccu[psave[4] + u * npu];
                            c5 = ccu[psave[5] + u * npu];
                            c6 = ccu[psave[6] + u * npu];
                            c7 = ccu[psave[7] + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] = c1 * x[j + kl1 *
                                                           n] +
                                    c2 * x[j + kl2 * n] + c3 * x[j +
                                                                      kl3 *
                                                                      n]
                                    + c4 * x[j + kl4 * n] + c5 * x[j +
                                                                        kl5 *
                                                                        n]
                                    + c6 * x[j + kl6 * n] + c7 * x[j +
                                                                        kl7 *
                                                                        n];
                            }
                        }
                        else if (nl == 8)
                        {
                            kl1 = ppair[1];
                            kl2 = ppair[2];
                            kl3 = ppair[3];
                            kl4 = ppair[4];
                            kl5 = ppair[5];
                            kl6 = ppair[6];
                            kl7 = ppair[7];
                            kl8 = ppair[8];
                            c1 = ccu[psave[1] + u * npu];
                            c2 = ccu[psave[2] + u * npu];
                            c3 = ccu[psave[3] + u * npu];
                            c4 = ccu[psave[4] + u * npu];
                            c5 = ccu[psave[5] + u * npu];
                            c6 = ccu[psave[6] + u * npu];
                            c7 = ccu[psave[7] + u * npu];
                            c8 = ccu[psave[8] + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] = c1 * x[j + kl1 *
                                                           n] +
                                    c2 * x[j + kl2 * n] + c3 * x[j +
                                                                      kl3 *
                                                                      n]
                                    + c4 * x[j + kl4 * n] + c5 * x[j +
                                                                        kl5 *
                                                                        n]
                                    + c6 * x[j + kl6 * n] + c7 * x[j +
                                                                        kl7 *
                                                                        n]
                                    + c8 * x[j + kl8 * n];
                            }
                        }
                        else
                        {
                            kl1 = ppair[1];
                            kl2 = ppair[2];
                            kl3 = ppair[3];
                            kl4 = ppair[4];
                            kl5 = ppair[5];
                            kl6 = ppair[6];
                            kl7 = ppair[7];
                            kl8 = ppair[8];
                            c1 = ccu[psave[1] + u * npu];
                            c2 = ccu[psave[2] + u * npu];
                            c3 = ccu[psave[3] + u * npu];
                            c4 = ccu[psave[4] + u * npu];
                            c5 = ccu[psave[5] + u * npu];
                            c6 = ccu[psave[6] + u * npu];
                            c7 = ccu[psave[7] + u * npu];
                            c8 = ccu[psave[8] + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                w[j + k * n] = c1 * x[j + kl1 *
                                                           n] +
                                    c2 * x[j + kl2 * n] + c3 * x[j +
                                                                      kl3 *
                                                                      n]
                                    + c4 * x[j + kl4 * n] + c5 * x[j +
                                                                        kl5 *
                                                                        n]
                                    + c6 * x[j + kl6 * n] + c7 * x[j +
                                                                        kl7 *
                                                                        n]
                                    + c8 * x[j + kl8 * n];
                            }
                            nlbase = 8;
                            nlleft = nl - 8;
                            nlstep = nlleft / 8;
                            nlrest = nlleft % 8;
                            for (m = 1; m <= nlstep; ++m)
                            {
                                kl1 = ppair[nlbase + 1];
                                kl2 = ppair[nlbase + 2];
                                kl3 = ppair[nlbase + 3];
                                kl4 = ppair[nlbase + 4];
                                kl5 = ppair[nlbase + 5];
                                kl6 = ppair[nlbase + 6];
                                kl7 = ppair[nlbase + 7];
                                kl8 = ppair[nlbase + 8];
                                c1 = ccu[psave[nlbase + 1] + u * npu];
                                c2 = ccu[psave[nlbase + 2] + u * npu];
                                c3 = ccu[psave[nlbase + 3] + u * npu];
                                c4 = ccu[psave[nlbase + 4] + u * npu];
                                c5 = ccu[psave[nlbase + 5] + u * npu];
                                c6 = ccu[psave[nlbase + 6] + u * npu];
                                c7 = ccu[psave[nlbase + 7] + u * npu];
                                c8 = ccu[psave[nlbase + 8] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] = w[j + k * n] +
                                        c1 * x[j + kl1 * n] +
                                        c2 * x[j + kl2 * n] +
                                        c3 * x[j + kl3 * n] +
                                        c4 * x[j + kl4 * n] +
                                        c5 * x[j + kl5 * n] +
                                        c6 * x[j + kl6 * n] +
                                        c7 * x[j + kl7 * n] +
                                        c8 * x[j + kl8 * n];
                                }
                                nlbase += 8;
                            }
                            for (m = 1; m <= nlrest; ++m)
                            {
                                kl1 = ppair[nlbase + m];
                                c1 = ccu[psave[nlbase + m] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] += c1 * x[j + kl1
                                                                * n];
                                }
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
                    if (nk == 0)
                    {
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = 0.;
                        }
                    }
                    else if (nk == 1)
                    {
                        k1 = psave[1];
                        c1 = cct[k1 + t * npt];
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = c1 * w[j + k1 * n];
                        }
                    }
                    else if (nk == 2)
                    {
                        k1 = psave[1];
                        k2 = psave[2];
                        c1 = cct[k1 + t * npt];
                        c2 = cct[k2 + t * npt];
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = c1 * w[j + k1 * n]
                                + c2 * w[j + k2 * n];
                        }
                    }
                    else if (nk == 3)
                    {
                        k1 = psave[1];
                        k2 = psave[2];
                        k3 = psave[3];
                        c1 = cct[k1 + t * npt];
                        c2 = cct[k2 + t * npt];
                        c3 = cct[k3 + t * npt];
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = c1 * w[j + k1 * n] +
                                c2 * w[j + k2 * n] +
                                c3 * w[j + k3 * n];
                        }
                    }
                    else if (nk == 4)
                    {
                        k1 = psave[1];
                        k2 = psave[2];
                        k3 = psave[3];
                        k4 = psave[4];
                        c1 = cct[k1 + t * npt];
                        c2 = cct[k2 + t * npt];
                        c3 = cct[k3 + t * npt];
                        c4 = cct[k4 + t * npt];
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = c1 * w[j + k1 * n]
                                + c2 * w[j + k2 * n] + c3 * w[j + k3 * n] +
                                c4 * w[j + k4 * n];
                        }
                    }
                    else if (nk == 5)
                    {
                        k1 = psave[1];
                        k2 = psave[2];
                        k3 = psave[3];
                        k4 = psave[4];
                        k5 = psave[5];
                        c1 = cct[k1 + t * npt];
                        c2 = cct[k2 + t * npt];
                        c3 = cct[k3 + t * npt];
                        c4 = cct[k4 + t * npt];
                        c5 = cct[k5 + t * npt];
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = c1 * w[j + k1 * n]
                                + c2 * w[j + k2 * n] + c3 * w[j +
                                                                   k3 *
                                                                   n] +
                                c4 * w[j + k4 * n] + c5 * w[j +
                                                                 k5 * n];
                        }
                    }
                    else if (nk == 6)
                    {
                        k1 = psave[1];
                        k2 = psave[2];
                        k3 = psave[3];
                        k4 = psave[4];
                        k5 = psave[5];
                        k6 = psave[6];
                        c1 = cct[k1 + t * npt];
                        c2 = cct[k2 + t * npt];
                        c3 = cct[k3 + t * npt];
                        c4 = cct[k4 + t * npt];
                        c5 = cct[k5 + t * npt];
                        c6 = cct[k6 + t * npt];
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = c1 * w[j + k1 * n]
                                + c2 * w[j + k2 * n] + c3 * w[j +
                                                                   k3 *
                                                                   n] +
                                c4 * w[j + k4 * n] + c5 * w[j +
                                                                 k5 *
                                                                 n] +
                                c6 * w[j + k6 * n];
                        }
                    }
                    else if (nk == 7)
                    {
                        k1 = psave[1];
                        k2 = psave[2];
                        k3 = psave[3];
                        k4 = psave[4];
                        k5 = psave[5];
                        k6 = psave[6];
                        k7 = psave[7];
                        c1 = cct[k1 + t * npt];
                        c2 = cct[k2 + t * npt];
                        c3 = cct[k3 + t * npt];
                        c4 = cct[k4 + t * npt];
                        c5 = cct[k5 + t * npt];
                        c6 = cct[k6 + t * npt];
                        c7 = cct[k7 + t * npt];
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = c1 * w[j + k1 * n]
                                + c2 * w[j + k2 * n] + c3 * w[j +
                                                                   k3 *
                                                                   n] +
                                c4 * w[j + k4 * n] + c5 * w[j +
                                                                 k5 *
                                                                 n] +
                                c6 * w[j + k6 * n] + c7 * w[j +
                                                                 k7 * n];
                        }
                    }
                    else if (nk == 8)
                    {
                        k1 = psave[1];
                        k2 = psave[2];
                        k3 = psave[3];
                        k4 = psave[4];
                        k5 = psave[5];
                        k6 = psave[6];
                        k7 = psave[7];
                        k8 = psave[8];
                        c1 = cct[k1 + t * npt];
                        c2 = cct[k2 + t * npt];
                        c3 = cct[k3 + t * npt];
                        c4 = cct[k4 + t * npt];
                        c5 = cct[k5 + t * npt];
                        c6 = cct[k6 + t * npt];
                        c7 = cct[k7 + t * npt];
                        c8 = cct[k8 + t * npt];
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = c1 * w[j + k1 * n]
                                + c2 * w[j + k2 * n] + c3 * w[j +
                                                                   k3 *
                                                                   n] +
                                c4 * w[j + k4 * n] + c5 * w[j +
                                                                 k5 *
                                                                 n] +
                                c6 * w[j + k6 * n] + c7 * w[j +
                                                                 k7 *
                                                                 n] +
                                c8 * w[j + k8 * n];
                        }
                    }
                    else
                    {
                        k1 = psave[1];
                        k2 = psave[2];
                        k3 = psave[3];
                        k4 = psave[4];
                        k5 = psave[5];
                        k6 = psave[6];
                        k7 = psave[7];
                        k8 = psave[8];
                        c1 = cct[k1 + t * npt];
                        c2 = cct[k2 + t * npt];
                        c3 = cct[k3 + t * npt];
                        c4 = cct[k4 + t * npt];
                        c5 = cct[k5 + t * npt];
                        c6 = cct[k6 + t * npt];
                        c7 = cct[k7 + t * npt];
                        c8 = cct[k8 + t * npt];
                        for (j = 1; j <= n; ++j)
                        {
                            y[j + tu * n] = c1 * w[j + k1 * n]
                                + c2 * w[j + k2 * n] + c3 * w[j +
                                                                   k3 *
                                                                   n] +
                                c4 * w[j + k4 * n] + c5 * w[j +
                                                                 k5 *
                                                                 n] +
                                c6 * w[j + k6 * n] + c7 * w[j +
                                                                 k7 *
                                                                 n] +
                                c8 * w[j + k8 * n];
                        }
                        nkbase = 8;
                        nkleft = nk - 8;
                        nkstep = nkleft / 8;
                        nkrest = nkleft % 8;
                        for (m = 1; m <= nkstep; ++m)
                        {
                            k1 = psave[nkbase + 1];
                            k2 = psave[nkbase + 2];
                            k3 = psave[nkbase + 3];
                            k4 = psave[nkbase + 4];
                            k5 = psave[nkbase + 5];
                            k6 = psave[nkbase + 6];
                            k7 = psave[nkbase + 7];
                            k8 = psave[nkbase + 8];
                            c1 = cct[k1 + t * npt];
                            c2 = cct[k2 + t * npt];
                            c3 = cct[k3 + t * npt];
                            c4 = cct[k4 + t * npt];
                            c5 = cct[k5 + t * npt];
                            c6 = cct[k6 + t * npt];
                            c7 = cct[k7 + t * npt];
                            c8 = cct[k8 + t * npt];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = y[j + tu *
                                                       n] + c1 * w[j +
                                                                        k1 *
                                                                        n]
                                    + c2 * w[j + k2 * n] + c3 * w[j +
                                                                       k3 *
                                                                       n]
                                    + c4 * w[j + k4 * n] + c5 * w[j +
                                                                       k5 *
                                                                       n]
                                    + c6 * w[j + k6 * n] + c7 * w[j +
                                                                       k7 *
                                                                       n]
                                    + c8 * w[j + k8 * n];
                            }
                            nkbase += 8;
                        }
                        for (m = 1; m <= nkrest; ++m)
                        {
                            k1 = psave[nkbase + m];
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
                            if (nk == 1)
                            {
                                kl1 = ppair[1];
                                c1 = cct[psave[1] + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + l * n] = c1 * x[j + kl1 * n];
                                }
                            }
                            else if (nk == 2)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                c1 = cct[psave[1] + t * npt];
                                c2 = cct[psave[2] + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + l * n] = c1 * x[j + kl1 * n] +
                                        c2 * x[j + kl2 * n];
                                }
                            }
                            else if (nk == 3)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                c1 = cct[psave[1] + t * npt];
                                c2 = cct[psave[2] + t * npt];
                                c3 = cct[psave[3] + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + l * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n];
                                }
                            }
                            else if (nk == 4)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                c1 = cct[psave[1] + t * npt];
                                c2 = cct[psave[2] + t * npt];
                                c3 = cct[psave[3] + t * npt];
                                c4 = cct[psave[4] + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + l * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n];
                                }
                            }
                            else if (nk == 5)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                kl5 = ppair[5];
                                c1 = cct[psave[1] + t * npt];
                                c2 = cct[psave[2] + t * npt];
                                c3 = cct[psave[3] + t * npt];
                                c4 = cct[psave[4] + t * npt];
                                c5 = cct[psave[5] + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + l * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n] +
                                        c5 * x[j + kl5 * n];
                                }
                            }
                            else if (nk == 6)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                kl5 = ppair[5];
                                kl6 = ppair[6];
                                c1 = cct[psave[1] + t * npt];
                                c2 = cct[psave[2] + t * npt];
                                c3 = cct[psave[3] + t * npt];
                                c4 = cct[psave[4] + t * npt];
                                c5 = cct[psave[5] + t * npt];
                                c6 = cct[psave[6] + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + l * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n] +
                                        c5 * x[j + kl5 * n] + c6 * x[j +
                                                                          kl6
                                                                          *
                                                                          n];
                                }
                            }
                            else if (nk == 7)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                kl5 = ppair[5];
                                kl6 = ppair[6];
                                kl7 = ppair[7];
                                c1 = cct[psave[1] + t * npt];
                                c2 = cct[psave[2] + t * npt];
                                c3 = cct[psave[3] + t * npt];
                                c4 = cct[psave[4] + t * npt];
                                c5 = cct[psave[5] + t * npt];
                                c6 = cct[psave[6] + t * npt];
                                c7 = cct[psave[7] + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + l * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n] +
                                        c5 * x[j + kl5 * n] + c6 * x[j +
                                                                          kl6
                                                                          *
                                                                          n]
                                        + c7 * x[j + kl7 * n];
                                }
                            }
                            else if (nk == 8)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                kl5 = ppair[5];
                                kl6 = ppair[6];
                                kl7 = ppair[7];
                                kl8 = ppair[8];
                                c1 = cct[psave[1] + t * npt];
                                c2 = cct[psave[2] + t * npt];
                                c3 = cct[psave[3] + t * npt];
                                c4 = cct[psave[4] + t * npt];
                                c5 = cct[psave[5] + t * npt];
                                c6 = cct[psave[6] + t * npt];
                                c7 = cct[psave[7] + t * npt];
                                c8 = cct[psave[8] + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + l * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n] +
                                        c5 * x[j + kl5 * n] + c6 * x[j +
                                                                          kl6
                                                                          *
                                                                          n]
                                        + c7 * x[j + kl7 * n] +
                                        c8 * x[j + kl8 * n];
                                }
                            }
                            else
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                kl5 = ppair[5];
                                kl6 = ppair[6];
                                kl7 = ppair[7];
                                kl8 = ppair[8];
                                c1 = cct[psave[1] + t * npt];
                                c2 = cct[psave[2] + t * npt];
                                c3 = cct[psave[3] + t * npt];
                                c4 = cct[psave[4] + t * npt];
                                c5 = cct[psave[5] + t * npt];
                                c6 = cct[psave[6] + t * npt];
                                c7 = cct[psave[7] + t * npt];
                                c8 = cct[psave[8] + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + l * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n] +
                                        c5 * x[j + kl5 * n] + c6 * x[j +
                                                                          kl6
                                                                          *
                                                                          n]
                                        + c7 * x[j + kl7 * n] +
                                        c8 * x[j + kl8 * n];
                                }
                                nkbase = 8;
                                nkleft = nk - 8;
                                nkstep = nkleft / 8;
                                nkrest = nkleft % 8;
                                for (m = 1; m <= nkstep; ++m)
                                {
                                    kl1 = ppair[nkbase + 1];
                                    kl2 = ppair[nkbase + 2];
                                    kl3 = ppair[nkbase + 3];
                                    kl4 = ppair[nkbase + 4];
                                    kl5 = ppair[nkbase + 5];
                                    kl6 = ppair[nkbase + 6];
                                    kl7 = ppair[nkbase + 7];
                                    kl8 = ppair[nkbase + 8];
                                    c1 = cct[psave[nkbase + 1] +
                                             t * npt];
                                    c2 = cct[psave[nkbase + 2] +
                                             t * npt];
                                    c3 = cct[psave[nkbase + 3] +
                                             t * npt];
                                    c4 = cct[psave[nkbase + 4] +
                                             t * npt];
                                    c5 = cct[psave[nkbase + 5] +
                                             t * npt];
                                    c6 = cct[psave[nkbase + 6] +
                                             t * npt];
                                    c7 = cct[psave[nkbase + 7] +
                                             t * npt];
                                    c8 = cct[psave[nkbase + 8] +
                                             t * npt];
                                    for (j = 1; j <= n; ++j)
                                    {
                                        w[j + l * n] = w[j + l * n]
                                            + c1 * x[j + kl1 *
                                                     n] + c2 * x[j +
                                                                      kl2 *
                                                                      n]
                                            + c3 * x[j + kl3 * n] +
                                            c4 * x[j + kl4 * n] +
                                            c5 * x[j + kl5 * n] +
                                            c6 * x[j + kl6 * n] +
                                            c7 * x[j + kl7 * n] +
                                            c8 * x[j + kl8 * n];
                                    }
                                    nkbase += 8;
                                }
                                for (m = 1; m <= nkrest; ++m)
                                {
                                    kl1 = ppair[nkbase + m];
                                    c1 = cct[psave[nkbase + m] +
                                             t * npt];
                                    for (j = 1; j <= n; ++j)
                                    {
                                        w[j + l * n] += c1 * x[j +
                                                                    kl1 *
                                                                    n];
                                    }
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
                        if (nl == 0)
                        {
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = 0.;
                            }
                        }
                        else if (nl == 1)
                        {
                            l1 = psave[1];
                            c1 = ccu[l1 + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + l1 * n];
                            }
                        }
                        else if (nl == 2)
                        {
                            l1 = psave[1];
                            l2 = psave[2];
                            c1 = ccu[l1 + u * npu];
                            c2 = ccu[l2 + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + l1 *
                                                            n] +
                                    c2 * w[j + l2 * n];
                            }
                        }
                        else if (nl == 3)
                        {
                            l1 = psave[1];
                            l2 = psave[2];
                            l3 = psave[3];
                            c1 = ccu[l1 + u * npu];
                            c2 = ccu[l2 + u * npu];
                            c3 = ccu[l3 + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + l1 *
                                                            n] +
                                    c2 * w[j + l2 * n] + c3 * w[j +
                                                                     l3 *
                                                                     n];
                            }
                        }
                        else if (nl == 4)
                        {
                            l1 = psave[1];
                            l2 = psave[2];
                            l3 = psave[3];
                            l4 = psave[4];
                            c1 = ccu[l1 + u * npu];
                            c2 = ccu[l2 + u * npu];
                            c3 = ccu[l3 + u * npu];
                            c4 = ccu[l4 + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + l1 *
                                                            n] +
                                    c2 * w[j + l2 * n] + c3 * w[j +
                                                                     l3 *
                                                                     n] +
                                    c4 * w[j + l4 * n];
                            }
                        }
                        else if (nl == 5)
                        {
                            l1 = psave[1];
                            l2 = psave[2];
                            l3 = psave[3];
                            l4 = psave[4];
                            l5 = psave[5];
                            c1 = ccu[l1 + u * npu];
                            c2 = ccu[l2 + u * npu];
                            c3 = ccu[l3 + u * npu];
                            c4 = ccu[l4 + u * npu];
                            c5 = ccu[l5 + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + l1 *
                                                            n] +
                                    c2 * w[j + l2 * n] + c3 * w[j +
                                                                     l3 *
                                                                     n] +
                                    c4 * w[j + l4 * n] + c5 * w[j +
                                                                     l5 *
                                                                     n];
                            }
                        }
                        else if (nl == 6)
                        {
                            l1 = psave[1];
                            l2 = psave[2];
                            l3 = psave[3];
                            l4 = psave[4];
                            l5 = psave[5];
                            l6 = psave[6];
                            c1 = ccu[l1 + u * npu];
                            c2 = ccu[l2 + u * npu];
                            c3 = ccu[l3 + u * npu];
                            c4 = ccu[l4 + u * npu];
                            c5 = ccu[l5 + u * npu];
                            c6 = ccu[l6 + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + l1 *
                                                            n] +
                                    c2 * w[j + l2 * n] + c3 * w[j +
                                                                     l3 *
                                                                     n] +
                                    c4 * w[j + l4 * n] + c5 * w[j +
                                                                     l5 *
                                                                     n] +
                                    c6 * w[j + l6 * n];
                            }
                        }
                        else if (nl == 7)
                        {
                            l1 = psave[1];
                            l2 = psave[2];
                            l3 = psave[3];
                            l4 = psave[4];
                            l5 = psave[5];
                            l6 = psave[6];
                            l7 = psave[7];
                            c1 = ccu[l1 + u * npu];
                            c2 = ccu[l2 + u * npu];
                            c3 = ccu[l3 + u * npu];
                            c4 = ccu[l4 + u * npu];
                            c5 = ccu[l5 + u * npu];
                            c6 = ccu[l6 + u * npu];
                            c7 = ccu[l7 + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + l1 *
                                                            n] +
                                    c2 * w[j + l2 * n] + c3 * w[j +
                                                                     l3 *
                                                                     n] +
                                    c4 * w[j + l4 * n] + c5 * w[j +
                                                                     l5 *
                                                                     n] +
                                    c6 * w[j + l6 * n] + c7 * w[j +
                                                                     l7 *
                                                                     n];
                            }
                        }
                        else if (nl == 8)
                        {
                            l1 = psave[1];
                            l2 = psave[2];
                            l3 = psave[3];
                            l4 = psave[4];
                            l5 = psave[5];
                            l6 = psave[6];
                            l7 = psave[7];
                            l8 = psave[8];
                            c1 = ccu[l1 + u * npu];
                            c2 = ccu[l2 + u * npu];
                            c3 = ccu[l3 + u * npu];
                            c4 = ccu[l4 + u * npu];
                            c5 = ccu[l5 + u * npu];
                            c6 = ccu[l6 + u * npu];
                            c7 = ccu[l7 + u * npu];
                            c8 = ccu[l8 + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + l1 *
                                                            n] +
                                    c2 * w[j + l2 * n] + c3 * w[j +
                                                                     l3 *
                                                                     n] +
                                    c4 * w[j + l4 * n] + c5 * w[j +
                                                                     l5 *
                                                                     n] +
                                    c6 * w[j + l6 * n] + c7 * w[j +
                                                                     l7 *
                                                                     n] +
                                    c8 * w[j + l8 * n];
                            }
                        }
                        else
                        {
                            l1 = psave[1];
                            l2 = psave[2];
                            l3 = psave[3];
                            l4 = psave[4];
                            l5 = psave[5];
                            l6 = psave[6];
                            l7 = psave[7];
                            l8 = psave[8];
                            c1 = ccu[l1 + u * npu];
                            c2 = ccu[l2 + u * npu];
                            c3 = ccu[l3 + u * npu];
                            c4 = ccu[l4 + u * npu];
                            c5 = ccu[l5 + u * npu];
                            c6 = ccu[l6 + u * npu];
                            c7 = ccu[l7 + u * npu];
                            c8 = ccu[l8 + u * npu];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + l1 *
                                                            n] +
                                    c2 * w[j + l2 * n] + c3 * w[j +
                                                                     l3 *
                                                                     n] +
                                    c4 * w[j + l4 * n] + c5 * w[j +
                                                                     l5 *
                                                                     n] +
                                    c6 * w[j + l6 * n] + c7 * w[j +
                                                                     l7 *
                                                                     n] +
                                    c8 * w[j + l8 * n];
                            }
                            nlbase = 8;
                            nlleft = nl - 8;
                            nlstep = nlleft / 8;
                            nlrest = nlleft % 8;
                            for (m = 1; m <= nlstep; ++m)
                            {
                                l1 = psave[nlbase + 1];
                                l2 = psave[nlbase + 2];
                                l3 = psave[nlbase + 3];
                                l4 = psave[nlbase + 4];
                                l5 = psave[nlbase + 5];
                                l6 = psave[nlbase + 6];
                                l7 = psave[nlbase + 7];
                                l8 = psave[nlbase + 8];
                                c1 = ccu[l1 + u * npu];
                                c2 = ccu[l2 + u * npu];
                                c3 = ccu[l3 + u * npu];
                                c4 = ccu[l4 + u * npu];
                                c5 = ccu[l5 + u * npu];
                                c6 = ccu[l6 + u * npu];
                                c7 = ccu[l7 + u * npu];
                                c8 = ccu[l8 + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    y[j + tu * n] = y[j + tu
                                                           * n] +
                                        c1 * w[j + l1 * n] + c2 * w[j +
                                                                         l2 *
                                                                         n]
                                        + c3 * w[j + l3 * n] + c4 * w[j +
                                                                           l4
                                                                           *
                                                                           n]
                                        + c5 * w[j + l5 * n] + c6 * w[j +
                                                                           l6
                                                                           *
                                                                           n]
                                        + c7 * w[j + l7 * n] + c8 * w[j +
                                                                           l8
                                                                           *
                                                                           n];
                                }
                                nlbase += 8;
                            }
                            for (m = 1; m <= nlrest; ++m)
                            {
                                l1 = psave[nlbase + m];
                                c1 = ccu[l1 + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    y[j + tu * n] += c1 * w[j + l1
                                                                 * n];
                                }
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
                            if (nl == 1)
                            {
                                kl1 = ppair[1];
                                c1 = ccu[psave[1] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] = c1 * x[j + kl1 *
                                                               n];
                                }
                            }
                            else if (nl == 2)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                c1 = ccu[psave[1] + u * npu];
                                c2 = ccu[psave[2] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n];
                                }
                            }
                            else if (nl == 3)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                c1 = ccu[psave[1] + u * npu];
                                c2 = ccu[psave[2] + u * npu];
                                c3 = ccu[psave[3] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n];
                                }
                            }
                            else if (nl == 4)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                c1 = ccu[psave[1] + u * npu];
                                c2 = ccu[psave[2] + u * npu];
                                c3 = ccu[psave[3] + u * npu];
                                c4 = ccu[psave[4] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n];
                                }
                            }
                            else if (nl == 5)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                kl5 = ppair[5];
                                c1 = ccu[psave[1] + u * npu];
                                c2 = ccu[psave[2] + u * npu];
                                c3 = ccu[psave[3] + u * npu];
                                c4 = ccu[psave[4] + u * npu];
                                c5 = ccu[psave[5] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n] +
                                        c5 * x[j + kl5 * n];
                                }
                            }
                            else if (nl == 6)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                kl5 = ppair[5];
                                kl6 = ppair[6];
                                c1 = ccu[psave[1] + u * npu];
                                c2 = ccu[psave[2] + u * npu];
                                c3 = ccu[psave[3] + u * npu];
                                c4 = ccu[psave[4] + u * npu];
                                c5 = ccu[psave[5] + u * npu];
                                c6 = ccu[psave[6] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n] +
                                        c5 * x[j + kl5 * n] + c6 * x[j +
                                                                          kl6
                                                                          *
                                                                          n];
                                }
                            }
                            else if (nl == 7)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                kl5 = ppair[5];
                                kl6 = ppair[6];
                                kl7 = ppair[7];
                                c1 = ccu[psave[1] + u * npu];
                                c2 = ccu[psave[2] + u * npu];
                                c3 = ccu[psave[3] + u * npu];
                                c4 = ccu[psave[4] + u * npu];
                                c5 = ccu[psave[5] + u * npu];
                                c6 = ccu[psave[6] + u * npu];
                                c7 = ccu[psave[7] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n] +
                                        c5 * x[j + kl5 * n] + c6 * x[j +
                                                                          kl6
                                                                          *
                                                                          n]
                                        + c7 * x[j + kl7 * n];
                                }
                            }
                            else if (nl == 8)
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                kl5 = ppair[5];
                                kl6 = ppair[6];
                                kl7 = ppair[7];
                                kl8 = ppair[8];
                                c1 = ccu[psave[1] + u * npu];
                                c2 = ccu[psave[2] + u * npu];
                                c3 = ccu[psave[3] + u * npu];
                                c4 = ccu[psave[4] + u * npu];
                                c5 = ccu[psave[5] + u * npu];
                                c6 = ccu[psave[6] + u * npu];
                                c7 = ccu[psave[7] + u * npu];
                                c8 = ccu[psave[8] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n] +
                                        c5 * x[j + kl5 * n] + c6 * x[j +
                                                                          kl6
                                                                          *
                                                                          n]
                                        + c7 * x[j + kl7 * n] +
                                        c8 * x[j + kl8 * n];
                                }
                            }
                            else
                            {
                                kl1 = ppair[1];
                                kl2 = ppair[2];
                                kl3 = ppair[3];
                                kl4 = ppair[4];
                                kl5 = ppair[5];
                                kl6 = ppair[6];
                                kl7 = ppair[7];
                                kl8 = ppair[8];
                                c1 = ccu[psave[1] + u * npu];
                                c2 = ccu[psave[2] + u * npu];
                                c3 = ccu[psave[3] + u * npu];
                                c4 = ccu[psave[4] + u * npu];
                                c5 = ccu[psave[5] + u * npu];
                                c6 = ccu[psave[6] + u * npu];
                                c7 = ccu[psave[7] + u * npu];
                                c8 = ccu[psave[8] + u * npu];
                                for (j = 1; j <= n; ++j)
                                {
                                    w[j + k * n] = c1 * x[j + kl1 *
                                                               n] +
                                        c2 * x[j + kl2 * n] + c3 * x[j +
                                                                          kl3
                                                                          *
                                                                          n]
                                        + c4 * x[j + kl4 * n] +
                                        c5 * x[j + kl5 * n] + c6 * x[j +
                                                                          kl6
                                                                          *
                                                                          n]
                                        + c7 * x[j + kl7 * n] +
                                        c8 * x[j + kl8 * n];
                                }
                                nlbase = 8;
                                nlleft = nl - 8;
                                nlstep = nlleft / 8;
                                nlrest = nlleft % 8;
                                for (m = 1; m <= nlstep; ++m)
                                {
                                    kl1 = ppair[nlbase + 1];
                                    kl2 = ppair[nlbase + 2];
                                    kl3 = ppair[nlbase + 3];
                                    kl4 = ppair[nlbase + 4];
                                    kl5 = ppair[nlbase + 5];
                                    kl6 = ppair[nlbase + 6];
                                    kl7 = ppair[nlbase + 7];
                                    kl8 = ppair[nlbase + 8];
                                    c1 = ccu[psave[nlbase + 1] +
                                             u * npu];
                                    c2 = ccu[psave[nlbase + 2] +
                                             u * npu];
                                    c3 = ccu[psave[nlbase + 3] +
                                             u * npu];
                                    c4 = ccu[psave[nlbase + 4] +
                                             u * npu];
                                    c5 = ccu[psave[nlbase + 5] +
                                             u * npu];
                                    c6 = ccu[psave[nlbase + 6] +
                                             u * npu];
                                    c7 = ccu[psave[nlbase + 7] +
                                             u * npu];
                                    c8 = ccu[psave[nlbase + 8] +
                                             u * npu];
                                    for (j = 1; j <= n; ++j)
                                    {
                                        w[j + k * n] = w[j + k * n]
                                            + c1 * x[j + kl1 *
                                                     n] + c2 * x[j +
                                                                      kl2 *
                                                                      n]
                                            + c3 * x[j + kl3 * n] +
                                            c4 * x[j + kl4 * n] +
                                            c5 * x[j + kl5 * n] +
                                            c6 * x[j + kl6 * n] +
                                            c7 * x[j + kl7 * n] +
                                            c8 * x[j + kl8 * n];
                                    }
                                    nlbase += 8;
                                }
                                for (m = 1; m <= nlrest; ++m)
                                {
                                    kl1 = ppair[nlbase + m];
                                    c1 = ccu[psave[nlbase + m] +
                                             u * npu];
                                    for (j = 1; j <= n; ++j)
                                    {
                                        w[j + k * n] += c1 * x[j +
                                                                    kl1 *
                                                                    n];
                                    }
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
                        if (nk == 0)
                        {
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = 0.;
                            }
                        }
                        else if (nk == 1)
                        {
                            k1 = psave[1];
                            c1 = cct[k1 + t * npt];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + k1 * n];
                            }
                        }
                        else if (nk == 2)
                        {
                            k1 = psave[1];
                            k2 = psave[2];
                            c1 = cct[k1 + t * npt];
                            c2 = cct[k2 + t * npt];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + k1 *
                                                            n] +
                                    c2 * w[j + k2 * n];
                            }
                        }
                        else if (nk == 3)
                        {
                            k1 = psave[1];
                            k2 = psave[2];
                            k3 = psave[3];
                            c1 = cct[k1 + t * npt];
                            c2 = cct[k2 + t * npt];
                            c3 = cct[k3 + t * npt];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + k1 *
                                                            n] +
                                    c2 * w[j + k2 * n] + c3 * w[j +
                                                                     k3 *
                                                                     n];
                            }
                        }
                        else if (nk == 4)
                        {
                            k1 = psave[1];
                            k2 = psave[2];
                            k3 = psave[3];
                            k4 = psave[4];
                            c1 = cct[k1 + t * npt];
                            c2 = cct[k2 + t * npt];
                            c3 = cct[k3 + t * npt];
                            c4 = cct[k4 + t * npt];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + k1 *
                                                            n] +
                                    c2 * w[j + k2 * n] + c3 * w[j +
                                                                     k3 *
                                                                     n] +
                                    c4 * w[j + k4 * n];
                            }
                        }
                        else if (nk == 5)
                        {
                            k1 = psave[1];
                            k2 = psave[2];
                            k3 = psave[3];
                            k4 = psave[4];
                            k5 = psave[5];
                            c1 = cct[k1 + t * npt];
                            c2 = cct[k2 + t * npt];
                            c3 = cct[k3 + t * npt];
                            c4 = cct[k4 + t * npt];
                            c5 = cct[k5 + t * npt];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + k1 *
                                                            n] +
                                    c2 * w[j + k2 * n] + c3 * w[j +
                                                                     k3 *
                                                                     n] +
                                    c4 * w[j + k4 * n] + c5 * w[j +
                                                                     k5 *
                                                                     n];
                            }
                        }
                        else if (nk == 6)
                        {
                            k1 = psave[1];
                            k2 = psave[2];
                            k3 = psave[3];
                            k4 = psave[4];
                            k5 = psave[5];
                            k6 = psave[6];
                            c1 = cct[k1 + t * npt];
                            c2 = cct[k2 + t * npt];
                            c3 = cct[k3 + t * npt];
                            c4 = cct[k4 + t * npt];
                            c5 = cct[k5 + t * npt];
                            c6 = cct[k6 + t * npt];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + k1 *
                                                            n] +
                                    c2 * w[j + k2 * n] + c3 * w[j +
                                                                     k3 *
                                                                     n] +
                                    c4 * w[j + k4 * n] + c5 * w[j +
                                                                     k5 *
                                                                     n] +
                                    c6 * w[j + k6 * n];
                            }
                        }
                        else if (nk == 7)
                        {
                            k1 = psave[1];
                            k2 = psave[2];
                            k3 = psave[3];
                            k4 = psave[4];
                            k5 = psave[5];
                            k6 = psave[6];
                            k7 = psave[7];
                            c1 = cct[k1 + t * npt];
                            c2 = cct[k2 + t * npt];
                            c3 = cct[k3 + t * npt];
                            c4 = cct[k4 + t * npt];
                            c5 = cct[k5 + t * npt];
                            c6 = cct[k6 + t * npt];
                            c7 = cct[k7 + t * npt];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + k1 *
                                                            n] +
                                    c2 * w[j + k2 * n] + c3 * w[j +
                                                                     k3 *
                                                                     n] +
                                    c4 * w[j + k4 * n] + c5 * w[j +
                                                                     k5 *
                                                                     n] +
                                    c6 * w[j + k6 * n] + c7 * w[j +
                                                                     k7 *
                                                                     n];
                            }
                        }
                        else if (nk == 8)
                        {
                            k1 = psave[1];
                            k2 = psave[2];
                            k3 = psave[3];
                            k4 = psave[4];
                            k5 = psave[5];
                            k6 = psave[6];
                            k7 = psave[7];
                            k8 = psave[8];
                            c1 = cct[k1 + t * npt];
                            c2 = cct[k2 + t * npt];
                            c3 = cct[k3 + t * npt];
                            c4 = cct[k4 + t * npt];
                            c5 = cct[k5 + t * npt];
                            c6 = cct[k6 + t * npt];
                            c7 = cct[k7 + t * npt];
                            c8 = cct[k8 + t * npt];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + k1 *
                                                            n] +
                                    c2 * w[j + k2 * n] + c3 * w[j +
                                                                     k3 *
                                                                     n] +
                                    c4 * w[j + k4 * n] + c5 * w[j +
                                                                     k5 *
                                                                     n] +
                                    c6 * w[j + k6 * n] + c7 * w[j +
                                                                     k7 *
                                                                     n] +
                                    c8 * w[j + k8 * n];
                            }
                        }
                        else
                        {
                            k1 = psave[1];
                            k2 = psave[2];
                            k3 = psave[3];
                            k4 = psave[4];
                            k5 = psave[5];
                            k6 = psave[6];
                            k7 = psave[7];
                            k8 = psave[8];
                            c1 = cct[k1 + t * npt];
                            c2 = cct[k2 + t * npt];
                            c3 = cct[k3 + t * npt];
                            c4 = cct[k4 + t * npt];
                            c5 = cct[k5 + t * npt];
                            c6 = cct[k6 + t * npt];
                            c7 = cct[k7 + t * npt];
                            c8 = cct[k8 + t * npt];
                            for (j = 1; j <= n; ++j)
                            {
                                y[j + tu * n] = c1 * w[j + k1 *
                                                            n] +
                                    c2 * w[j + k2 * n] + c3 * w[j +
                                                                     k3 *
                                                                     n] +
                                    c4 * w[j + k4 * n] + c5 * w[j +
                                                                     k5 *
                                                                     n] +
                                    c6 * w[j + k6 * n] + c7 * w[j +
                                                                     k7 *
                                                                     n] +
                                    c8 * w[j + k8 * n];
                            }
                            nkbase = 8;
                            nkleft = nk - 8;
                            nkstep = nkleft / 8;
                            nkrest = nkleft % 8;
                            for (m = 1; m <= nkstep; ++m)
                            {
                                k1 = psave[nkbase + 1];
                                k2 = psave[nkbase + 2];
                                k3 = psave[nkbase + 3];
                                k4 = psave[nkbase + 4];
                                k5 = psave[nkbase + 5];
                                k6 = psave[nkbase + 6];
                                k7 = psave[nkbase + 7];
                                k8 = psave[nkbase + 8];
                                c1 = cct[k1 + t * npt];
                                c2 = cct[k2 + t * npt];
                                c3 = cct[k3 + t * npt];
                                c4 = cct[k4 + t * npt];
                                c5 = cct[k5 + t * npt];
                                c6 = cct[k6 + t * npt];
                                c7 = cct[k7 + t * npt];
                                c8 = cct[k8 + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    y[j + tu * n] = y[j + tu
                                                           * n] +
                                        c1 * w[j + k1 * n] + c2 * w[j +
                                                                         k2 *
                                                                         n]
                                        + c3 * w[j + k3 * n] + c4 * w[j +
                                                                           k4
                                                                           *
                                                                           n]
                                        + c5 * w[j + k5 * n] + c6 * w[j +
                                                                           k6
                                                                           *
                                                                           n]
                                        + c7 * w[j + k7 * n] + c8 * w[j +
                                                                           k8
                                                                           *
                                                                           n];
                                }
                                nkbase += 8;
                            }
                            for (m = 1; m <= nkrest; ++m)
                            {
                                k1 = psave[nkbase + m];
                                c1 = cct[k1 + t * npt];
                                for (j = 1; j <= n; ++j)
                                {
                                    y[j + tu * n] += c1 * w[j + k1
                                                                 * n];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}
