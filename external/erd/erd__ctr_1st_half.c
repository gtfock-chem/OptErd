#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


int erd__ctr_1st_half (int n, int npmax, int npmin,
                       int mij, int nrs, int nblock,
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
    double c2;
    double c3;
    double c4;
    double c5;
    double c6;
    double c7;
    double c8;
    int i1;
    int i2;
    int i3;
    int i4;
    int i5;
    int i6;
    int i7;
    int i8;
    int j1;
    int j2;
    int j3;
    int j4;
    int j5;
    int j6;
    int j7;
    int j8;
    int ij;
    int ni;
    int nj;
    int rs;
    int ij1;
    int ij2;
    int ij3;
    int ij4;
    int ij5;
    int ij6;
    int ij7;
    int ij8;
    int pmin;
    int pmax;
    int inext;
    int jnext;
    int nibase;
    int njbase;
    int irange;
    int jrange;
    int nileft;
    int njleft;
    int cntrct;
    int nistep;
    int njstep;
    int nirest;
    int njrest;

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
        rs = 0;
        for (s = 1; s <= ncs; s++)
        {
            pmin = ccbegs[s];
            pmax = ccends[s];
            for (i = 1; i <= npr; i++)
            {
                pused[i] = 0;
            }
            if (pmin == pmax && ccs[pmin + s * nps] == 1.0)
            {
/*             ...we reach this point, if only one S contraction */
/*                coefficient equal to 1.0 is present. */
                for (ij = 1; ij <= mij; ij++)
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
                    jrange = ((j >= pmin) && (j <= pmax));
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
                        if (nj == 1)
                        {
                            ij1 = ppair[1];
                            c1 = ccs[psave[1] + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] = c1 * x[l + ij1 * n];
                            }
                        }
                        else if (nj == 2)
                        {
                            ij1 = ppair[1];
                            ij2 = ppair[2];
                            c1 = ccs[psave[1] + s * nps];
                            c2 = ccs[psave[2] + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] = c1 * x[l + ij1 * n] +
                                    c2 * x[l + ij2 * n];
                            }
                        }
                        else if (nj == 3)
                        {
                            ij1 = ppair[1];
                            ij2 = ppair[2];
                            ij3 = ppair[3];
                            c1 = ccs[psave[1] + s * nps];
                            c2 = ccs[psave[2] + s * nps];
                            c3 = ccs[psave[3] + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] = c1 * x[l + ij1 * n] +
                                    c2 * x[l + ij2 * n] +
                                    c3 * x[l + ij3 * n];
                            }
                        }
                        else if (nj == 4)
                        {
                            ij1 = ppair[1];
                            ij2 = ppair[2];
                            ij3 = ppair[3];
                            ij4 = ppair[4];
                            c1 = ccs[psave[1] + s * nps];
                            c2 = ccs[psave[2] + s * nps];
                            c3 = ccs[psave[3] + s * nps];
                            c4 = ccs[psave[4] + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] = c1 * x[l + ij1 * n] +
                                    c2 * x[l + ij2 * n] +
                                    c3 * x[l + ij3 * n] +
                                    c4 * x[l + ij4 * n];
                            }
                        }
                        else if (nj == 5)
                        {
                            ij1 = ppair[1];
                            ij2 = ppair[2];
                            ij3 = ppair[3];
                            ij4 = ppair[4];
                            ij5 = ppair[5];
                            c1 = ccs[psave[1] + s * nps];
                            c2 = ccs[psave[2] + s * nps];
                            c3 = ccs[psave[3] + s * nps];
                            c4 = ccs[psave[4] + s * nps];
                            c5 = ccs[psave[5] + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] = c1 * x[l + ij1 * n] +
                                    c2 * x[l + ij2 * n] +
                                    c3 * x[l + ij3 * n] +
                                    c4 * x[l + ij4 * n] +
                                    c5 * x[l + ij5 * n];
                            }
                        }
                        else if (nj == 6)
                        {
                            ij1 = ppair[1];
                            ij2 = ppair[2];
                            ij3 = ppair[3];
                            ij4 = ppair[4];
                            ij5 = ppair[5];
                            ij6 = ppair[6];
                            c1 = ccs[psave[1] + s * nps];
                            c2 = ccs[psave[2] + s * nps];
                            c3 = ccs[psave[3] + s * nps];
                            c4 = ccs[psave[4] + s * nps];
                            c5 = ccs[psave[5] + s * nps];
                            c6 = ccs[psave[6] + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] = c1 * x[l + ij1 * n] +
                                    c2 * x[l + ij2 * n] +
                                    c3 * x[l + ij3 * n] +
                                    c4 * x[l + ij4 * n] +
                                    c5 * x[l + ij5 * n] +
                                    c6 * x[l + ij6 * n];
                            }
                        }
                        else if (nj == 7)
                        {
                            ij1 = ppair[1];
                            ij2 = ppair[2];
                            ij3 = ppair[3];
                            ij4 = ppair[4];
                            ij5 = ppair[5];
                            ij6 = ppair[6];
                            ij7 = ppair[7];
                            c1 = ccs[psave[1] + s * nps];
                            c2 = ccs[psave[2] + s * nps];
                            c3 = ccs[psave[3] + s * nps];
                            c4 = ccs[psave[4] + s * nps];
                            c5 = ccs[psave[5] + s * nps];
                            c6 = ccs[psave[6] + s * nps];
                            c7 = ccs[psave[7] + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] = c1 * x[l + ij1 * n] +
                                    c2 * x[l + ij2 * n] +
                                    c3 * x[l + ij3 * n] +
                                    c4 * x[l + ij4 * n] +
                                    c5 * x[l + ij5 * n] + 
                                    c6 * x[l + ij6 * n] +
                                    c7 * x[l + ij7 * n];
                            }
                        }
                        else if (nj == 8)
                        {
                            ij1 = ppair[1];
                            ij2 = ppair[2];
                            ij3 = ppair[3];
                            ij4 = ppair[4];
                            ij5 = ppair[5];
                            ij6 = ppair[6];
                            ij7 = ppair[7];
                            ij8 = ppair[8];
                            c1 = ccs[psave[1] + s * nps];
                            c2 = ccs[psave[2] + s * nps];
                            c3 = ccs[psave[3] + s * nps];
                            c4 = ccs[psave[4] + s * nps];
                            c5 = ccs[psave[5] + s * nps];
                            c6 = ccs[psave[6] + s * nps];
                            c7 = ccs[psave[7] + s * nps];
                            c8 = ccs[psave[8] + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] = c1 * x[l + ij1 * n] +
                                    c2 * x[l + ij2 * n] +
                                    c3 * x[l + ij3 * n] +
                                    c4 * x[l + ij4 * n] +
                                    c5 * x[l + ij5 * n] +
                                    c6 * x[l + ij6 * n] +
                                    c7 * x[l + ij7 * n] +
                                    c8 * x[l + ij8 * n];
                            }
                        }
                        else
                        {
                            ij1 = ppair[1];
                            ij2 = ppair[2];
                            ij3 = ppair[3];
                            ij4 = ppair[4];
                            ij5 = ppair[5];
                            ij6 = ppair[6];
                            ij7 = ppair[7];
                            ij8 = ppair[8];
                            c1 = ccs[psave[1] + s * nps];
                            c2 = ccs[psave[2] + s * nps];
                            c3 = ccs[psave[3] + s * nps];
                            c4 = ccs[psave[4] + s * nps];
                            c5 = ccs[psave[5] + s * nps];
                            c6 = ccs[psave[6] + s * nps];
                            c7 = ccs[psave[7] + s * nps];
                            c8 = ccs[psave[8] + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                w[l + i * n] = c1 * x[l + ij1 * n] +
                                    c2 * x[l + ij2 * n] +
                                    c3 * x[l + ij3 * n] +
                                    c4 * x[l + ij4 * n] +
                                    c5 * x[l + ij5 * n] +
                                    c6 * x[l + ij6 * n] +
                                    c7 * x[l + ij7 * n] +
                                    c8 * x[l + ij8 * n];
                            }
                            njbase = 8;
                            njleft = nj - 8;
                            njstep = njleft / 8;
                            njrest = njleft % 8;
                            for (m = 1; m <= njstep; ++m)
                            {
                                ij1 = ppair[njbase + 1];
                                ij2 = ppair[njbase + 2];
                                ij3 = ppair[njbase + 3];
                                ij4 = ppair[njbase + 4];
                                ij5 = ppair[njbase + 5];
                                ij6 = ppair[njbase + 6];
                                ij7 = ppair[njbase + 7];
                                ij8 = ppair[njbase + 8];
                                c1 = ccs[psave[njbase + 1] + s * nps];
                                c2 = ccs[psave[njbase + 2] + s * nps];
                                c3 = ccs[psave[njbase + 3] + s * nps];
                                c4 = ccs[psave[njbase + 4] + s * nps];
                                c5 = ccs[psave[njbase + 5] + s * nps];
                                c6 = ccs[psave[njbase + 6] + s * nps];
                                c7 = ccs[psave[njbase + 7] + s * nps];
                                c8 = ccs[psave[njbase + 8] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] = w[l + i * n] +
                                        c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n] +
                                        c5 * x[l + ij5 * n] +
                                        c6 * x[l + ij6 * n] +
                                        c7 * x[l + ij7 * n] +
                                        c8 * x[l + ij8 * n];
                                }
                                njbase += 8;
                            }
                            for (m = 1; m <= njrest; ++m)
                            {
                                ij1 = ppair[njbase + m];
                                c1 = ccs[psave[njbase + m] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] += c1 * x[l + ij1 * n];
                                }
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
                        if (i >= pmin && i <= pmax)
                        {
                            c1 = ccs[i + s * nps];
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
            for (r = s; r <= ncr; r++)
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
                    for (i = pmin; i <= pmax; ++i)
                    {
                        if (pused[i] == 1)
                        {
                            ni++;
                            psave[ni] = i;
                        }
                    }
                    if (ni == 0)
                    {
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = 0.;
                        }
                    }
                    else if (ni == 1)
                    {
                        i1 = psave[1];
                        c1 = ccr[i1 + r * npr];
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = c1 * w[l + i1 * n];
                        }
                    }
                    else if (ni == 2)
                    {
                        i1 = psave[1];
                        i2 = psave[2];
                        c1 = ccr[i1 + r * npr];
                        c2 = ccr[i2 + r * npr];
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = c1 * w[l + i1 * n]
                                + c2 * w[l + i2 * n];
                        }
                    }
                    else if (ni == 3)
                    {
                        i1 = psave[1];
                        i2 = psave[2];
                        i3 = psave[3];
                        c1 = ccr[i1 + r * npr];
                        c2 = ccr[i2 + r * npr];
                        c3 = ccr[i3 + r * npr];
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = c1 * w[l + i1 * n] +
                                c2 * w[l + i2 * n] +
                                c3 * w[l + i3 * n];
                        }
                    }
                    else if (ni == 4)
                    {
                        i1 = psave[1];
                        i2 = psave[2];
                        i3 = psave[3];
                        i4 = psave[4];
                        c1 = ccr[i1 + r * npr];
                        c2 = ccr[i2 + r * npr];
                        c3 = ccr[i3 + r * npr];
                        c4 = ccr[i4 + r * npr];
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = c1 * w[l + i1 * n] +
                                c2 * w[l + i2 * n] +
                                c3 * w[l + i3 * n] +
                                c4 * w[l + i4 * n];
                        }
                    }
                    else if (ni == 5)
                    {
                        i1 = psave[1];
                        i2 = psave[2];
                        i3 = psave[3];
                        i4 = psave[4];
                        i5 = psave[5];
                        c1 = ccr[i1 + r * npr];
                        c2 = ccr[i2 + r * npr];
                        c3 = ccr[i3 + r * npr];
                        c4 = ccr[i4 + r * npr];
                        c5 = ccr[i5 + r * npr];
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = c1 * w[l + i1 * n] +
                                c2 * w[l + i2 * n] +
                                c3 * w[l + i3 * n] +
                                c4 * w[l + i4 * n] +
                                c5 * w[l + i5 * n];
                        }
                    }
                    else if (ni == 6)
                    {
                        i1 = psave[1];
                        i2 = psave[2];
                        i3 = psave[3];
                        i4 = psave[4];
                        i5 = psave[5];
                        i6 = psave[6];
                        c1 = ccr[i1 + r * npr];
                        c2 = ccr[i2 + r * npr];
                        c3 = ccr[i3 + r * npr];
                        c4 = ccr[i4 + r * npr];
                        c5 = ccr[i5 + r * npr];
                        c6 = ccr[i6 + r * npr];
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = c1 * w[l + i1 * n] +
                                c2 * w[l + i2 * n] +
                                c3 * w[l + i3 * n] +
                                c4 * w[l + i4 * n] +
                                c5 * w[l + i5 * n] +
                                c6 * w[l + i6 * n];
                        }
                    }
                    else if (ni == 7)
                    {
                        i1 = psave[1];
                        i2 = psave[2];
                        i3 = psave[3];
                        i4 = psave[4];
                        i5 = psave[5];
                        i6 = psave[6];
                        i7 = psave[7];
                        c1 = ccr[i1 + r * npr];
                        c2 = ccr[i2 + r * npr];
                        c3 = ccr[i3 + r * npr];
                        c4 = ccr[i4 + r * npr];
                        c5 = ccr[i5 + r * npr];
                        c6 = ccr[i6 + r * npr];
                        c7 = ccr[i7 + r * npr];
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = c1 * w[l + i1 * n] +
                                c2 * w[l + i2 * n] +
                                c3 * w[l + i3 * n] +
                                c4 * w[l + i4 * n] +
                                c5 * w[l + i5 * n] +
                                c6 * w[l + i6 * n] +
                                c7 * w[l + i7 * n];
                        }
                    }
                    else if (ni == 8)
                    {
                        i1 = psave[1];
                        i2 = psave[2];
                        i3 = psave[3];
                        i4 = psave[4];
                        i5 = psave[5];
                        i6 = psave[6];
                        i7 = psave[7];
                        i8 = psave[8];
                        c1 = ccr[i1 + r * npr];
                        c2 = ccr[i2 + r * npr];
                        c3 = ccr[i3 + r * npr];
                        c4 = ccr[i4 + r * npr];
                        c5 = ccr[i5 + r * npr];
                        c6 = ccr[i6 + r * npr];
                        c7 = ccr[i7 + r * npr];
                        c8 = ccr[i8 + r * npr];
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = c1 * w[l + i1 * n] +
                                c2 * w[l + i2 * n] +
                                c3 * w[l + i3 * n] +
                                c4 * w[l + i4 * n] +
                                c5 * w[l + i5 * n] +
                                c6 * w[l + i6 * n] +
                                c7 * w[l + i7 * n] +
                                c8 * w[l + i8 * n];
                        }
                    }
                    else if (ni > 8)
                    {
                        i1 = psave[1];
                        i2 = psave[2];
                        i3 = psave[3];
                        i4 = psave[4];
                        i5 = psave[5];
                        i6 = psave[6];
                        i7 = psave[7];
                        i8 = psave[8];
                        c1 = ccr[i1 + r * npr];
                        c2 = ccr[i2 + r * npr];
                        c3 = ccr[i3 + r * npr];
                        c4 = ccr[i4 + r * npr];
                        c5 = ccr[i5 + r * npr];
                        c6 = ccr[i6 + r * npr];
                        c7 = ccr[i7 + r * npr];
                        c8 = ccr[i8 + r * npr];
                        for (l = 1; l <= n; l++)
                        {
                            y[l + rs * n] = c1 * w[l + i1 * n] +
                                c2 * w[l + i2 * n] +
                                c3 * w[l + i3 * n] +
                                c4 * w[l + i4 * n] +
                                c5 * w[l + i5 * n] +
                                c6 * w[l + i6 * n] +
                                c7 * w[l + i7 * n] +
                                c8 * w[l + i8 * n];
                        }
                        nibase = 8;
                        nileft = ni - 8;
                        nistep = nileft / 8;
                        nirest = nileft % 8;
                        for (m = 1; m <= nistep; ++m)
                        {
                            i1 = psave[nibase + 1];
                            i2 = psave[nibase + 2];
                            i3 = psave[nibase + 3];
                            i4 = psave[nibase + 4];
                            i5 = psave[nibase + 5];
                            i6 = psave[nibase + 6];
                            i7 = psave[nibase + 7];
                            i8 = psave[nibase + 8];
                            c1 = ccr[i1 + r * npr];
                            c2 = ccr[i2 + r * npr];
                            c3 = ccr[i3 + r * npr];
                            c4 = ccr[i4 + r * npr];
                            c5 = ccr[i5 + r * npr];
                            c6 = ccr[i6 + r * npr];
                            c7 = ccr[i7 + r * npr];
                            c8 = ccr[i8 + r * npr];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = y[l + rs * n] +
                                    c1 * w[l + i1 * n] +
                                    c2 * w[l + i2 * n] +
                                    c3 * w[l + i3 * n] +
                                    c4 * w[l + i4 * n] +
                                    c5 * w[l + i5 * n] + 
                                    c6 * w[l + i6 * n] +
                                    c7 * w[l + i7 * n] +
                                    c8 * w[l + i8 * n];
                            }
                            nibase += 8;
                        }
                        i3 = nirest;
                        for (m = 1; m <= i3; ++m)
                        {
                            i1 = psave[nibase + m];
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
                            if (ni == 1)
                            {
                                ij1 = ppair[1];
                                c1 = ccr[psave[1] + r * npr];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + j * n] = c1 * x[l + ij1 * n];
                                }
                            }
                            else if (ni == 2)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                c1 = ccr[psave[1] + r * npr];
                                c2 = ccr[psave[2] + r * npr];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + j * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n];
                                }
                            }
                            else if (ni == 3)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                c1 = ccr[psave[1] + r * npr];
                                c2 = ccr[psave[2] + r * npr];
                                c3 = ccr[psave[3] + r * npr];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + j * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n];
                                }
                            }
                            else if (ni == 4)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                c1 = ccr[psave[1] + r * npr];
                                c2 = ccr[psave[2] + r * npr];
                                c3 = ccr[psave[3] + r * npr];
                                c4 = ccr[psave[4] + r * npr];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + j * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n];
                                }
                            }
                            else if (ni == 5)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                ij5 = ppair[5];
                                c1 = ccr[psave[1] + r * npr];
                                c2 = ccr[psave[2] + r * npr];
                                c3 = ccr[psave[3] + r * npr];
                                c4 = ccr[psave[4] + r * npr];
                                c5 = ccr[psave[5] + r * npr];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + j * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n] +
                                        c5 * x[l + ij5 * n];
                                }
                            }
                            else if (ni == 6)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                ij5 = ppair[5];
                                ij6 = ppair[6];
                                c1 = ccr[psave[1] + r * npr];
                                c2 = ccr[psave[2] + r * npr];
                                c3 = ccr[psave[3] + r * npr];
                                c4 = ccr[psave[4] + r * npr];
                                c5 = ccr[psave[5] + r * npr];
                                c6 = ccr[psave[6] + r * npr];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + j * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n] +
                                        c5 * x[l + ij5 * n] +
                                        c6 * x[l + ij6 * n];
                                }
                            }
                            else if (ni == 7)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                ij5 = ppair[5];
                                ij6 = ppair[6];
                                ij7 = ppair[7];
                                c1 = ccr[psave[1] + r * npr];
                                c2 = ccr[psave[2] + r * npr];
                                c3 = ccr[psave[3] + r * npr];
                                c4 = ccr[psave[4] + r * npr];
                                c5 = ccr[psave[5] + r * npr];
                                c6 = ccr[psave[6] + r * npr];
                                c7 = ccr[psave[7] + r * npr];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + j * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n] +
                                        c5 * x[l + ij5 * n] +
                                        c6 * x[l + ij6 * n] +
                                        c7 * x[l + ij7 * n];
                                }
                            }
                            else if (ni == 8)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                ij5 = ppair[5];
                                ij6 = ppair[6];
                                ij7 = ppair[7];
                                ij8 = ppair[8];
                                c1 = ccr[psave[1] + r * npr];
                                c2 = ccr[psave[2] + r * npr];
                                c3 = ccr[psave[3] + r * npr];
                                c4 = ccr[psave[4] + r * npr];
                                c5 = ccr[psave[5] + r * npr];
                                c6 = ccr[psave[6] + r * npr];
                                c7 = ccr[psave[7] + r * npr];
                                c8 = ccr[psave[8] + r * npr];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + j * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n] +
                                        c5 * x[l + ij5 * n] +
                                        c6 * x[l + ij6 * n] +
                                        c7 * x[l + ij7 * n] +
                                        c8 * x[l + ij8 * n];
                                }
                            }
                            else
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                ij5 = ppair[5];
                                ij6 = ppair[6];
                                ij7 = ppair[7];
                                ij8 = ppair[8];
                                c1 = ccr[psave[1] + r * npr];
                                c2 = ccr[psave[2] + r * npr];
                                c3 = ccr[psave[3] + r * npr];
                                c4 = ccr[psave[4] + r * npr];
                                c5 = ccr[psave[5] + r * npr];
                                c6 = ccr[psave[6] + r * npr];
                                c7 = ccr[psave[7] + r * npr];
                                c8 = ccr[psave[8] + r * npr];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + j * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n] +
                                        c5 * x[l + ij5 * n] +
                                        c6 * x[l + ij6 * n] +
                                        c7 * x[l + ij7 * n] +
                                        c8 * x[l + ij8 * n];
                                }
                                nibase = 8;
                                nileft = ni - 8;
                                nistep = nileft / 8;
                                nirest = nileft % 8;
                                i3 = nistep;
                                for (m = 1; m <= nistep; ++m)
                                {
                                    ij1 = ppair[nibase + 1];
                                    ij2 = ppair[nibase + 2];
                                    ij3 = ppair[nibase + 3];
                                    ij4 = ppair[nibase + 4];
                                    ij5 = ppair[nibase + 5];
                                    ij6 = ppair[nibase + 6];
                                    ij7 = ppair[nibase + 7];
                                    ij8 = ppair[nibase + 8];
                                    c1 = ccr[psave[nibase + 1] + r * npr];
                                    c2 = ccr[psave[nibase + 2] + r * npr];
                                    c3 = ccr[psave[nibase + 3] + r * npr];
                                    c4 = ccr[psave[nibase + 4] + r * npr];
                                    c5 = ccr[psave[nibase + 5] + r * npr];
                                    c6 = ccr[psave[nibase + 6] + r * npr];
                                    c7 = ccr[psave[nibase + 7] + r * npr];
                                    c8 = ccr[psave[nibase + 8] + r * npr];
                                    for (l = 1; l <= n; l++)
                                    {
                                        w[l + j * n] = w[l + j * n] +
                                            c1 * x[l + ij1 * n] +
                                            c2 * x[l + ij2 * n] +
                                            c3 * x[l + ij3 * n] +
                                            c4 * x[l + ij4 * n] +
                                            c5 * x[l + ij5 * n] +
                                            c6 * x[l + ij6 * n] +
                                            c7 * x[l + ij7 * n] +
                                            c8 * x[l + ij8 * n];
                                    }
                                    nibase += 8;
                                }
                                for (m = 1; m <= nirest; ++m)
                                {
                                    ij1 = ppair[nibase + m];
                                    c1 = ccr[psave[nibase + m] + r * npr];
                                    for (l = 1; l <= n; l++)
                                    {
                                        w[l + j * n] += c1 * x[l + ij1 * n];
                                    }
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
                        if (nj == 0)
                        {
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = 0.;
                            }
                        }
                        else if (nj == 1)
                        {
                            j1 = psave[1];
                            c1 = ccs[j1 + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + j1 * n];
                            }
                        }
                        else if (nj == 2)
                        {
                            j1 = psave[1];
                            j2 = psave[2];
                            c1 = ccs[j1 + s * nps];
                            c2 = ccs[j2 + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + j1 * n] +
                                    c2 * w[l + j2 * n];
                            }
                        }
                        else if (nj == 3)
                        {
                            j1 = psave[1];
                            j2 = psave[2];
                            j3 = psave[3];
                            c1 = ccs[j1 + s * nps];
                            c2 = ccs[j2 + s * nps];
                            c3 = ccs[j3 + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + j1 * n] +
                                    c2 * w[l + j2 * n] +
                                    c3 * w[l + j3 * n];
                            }
                        }
                        else if (nj == 4)
                        {
                            j1 = psave[1];
                            j2 = psave[2];
                            j3 = psave[3];
                            j4 = psave[4];
                            c1 = ccs[j1 + s * nps];
                            c2 = ccs[j2 + s * nps];
                            c3 = ccs[j3 + s * nps];
                            c4 = ccs[j4 + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + j1 * n] +
                                    c2 * w[l + j2 * n] +
                                    c3 * w[l + j3 * n] +
                                    c4 * w[l + j4 * n];
                            }
                        }
                        else if (nj == 5)
                        {
                            j1 = psave[1];
                            j2 = psave[2];
                            j3 = psave[3];
                            j4 = psave[4];
                            j5 = psave[5];
                            c1 = ccs[j1 + s * nps];
                            c2 = ccs[j2 + s * nps];
                            c3 = ccs[j3 + s * nps];
                            c4 = ccs[j4 + s * nps];
                            c5 = ccs[j5 + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + j1 * n] +
                                    c2 * w[l + j2 * n] +
                                    c3 * w[l + j3 * n] +
                                    c4 * w[l + j4 * n] +
                                    c5 * w[l + j5 * n];
                            }
                        }
                        else if (nj == 6)
                        {
                            j1 = psave[1];
                            j2 = psave[2];
                            j3 = psave[3];
                            j4 = psave[4];
                            j5 = psave[5];
                            j6 = psave[6];
                            c1 = ccs[j1 + s * nps];
                            c2 = ccs[j2 + s * nps];
                            c3 = ccs[j3 + s * nps];
                            c4 = ccs[j4 + s * nps];
                            c5 = ccs[j5 + s * nps];
                            c6 = ccs[j6 + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + j1 * n] +
                                    c2 * w[l + j2 * n] +
                                    c3 * w[l + j3 * n] +
                                    c4 * w[l + j4 * n] +
                                    c5 * w[l + j5 * n] +
                                    c6 * w[l + j6 * n];
                            }
                        }
                        else if (nj == 7)
                        {
                            j1 = psave[1];
                            j2 = psave[2];
                            j3 = psave[3];
                            j4 = psave[4];
                            j5 = psave[5];
                            j6 = psave[6];
                            j7 = psave[7];
                            c1 = ccs[j1 + s * nps];
                            c2 = ccs[j2 + s * nps];
                            c3 = ccs[j3 + s * nps];
                            c4 = ccs[j4 + s * nps];
                            c5 = ccs[j5 + s * nps];
                            c6 = ccs[j6 + s * nps];
                            c7 = ccs[j7 + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + j1 * n] +
                                    c2 * w[l + j2 * n] +
                                    c3 * w[l + j3 * n] +
                                    c4 * w[l + j4 * n] +
                                    c5 * w[l + j5 * n] +
                                    c6 * w[l + j6 * n] +
                                    c7 * w[l + j7 * n];
                            }
                        }
                        else if (nj == 8)
                        {
                            j1 = psave[1];
                            j2 = psave[2];
                            j3 = psave[3];
                            j4 = psave[4];
                            j5 = psave[5];
                            j6 = psave[6];
                            j7 = psave[7];
                            j8 = psave[8];
                            c1 = ccs[j1 + s * nps];
                            c2 = ccs[j2 + s * nps];
                            c3 = ccs[j3 + s * nps];
                            c4 = ccs[j4 + s * nps];
                            c5 = ccs[j5 + s * nps];
                            c6 = ccs[j6 + s * nps];
                            c7 = ccs[j7 + s * nps];
                            c8 = ccs[j8 + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + j1 * n] +
                                    c2 * w[l + j2 * n] +
                                    c3 * w[l + j3 * n] +
                                    c4 * w[l + j4 * n] +
                                    c5 * w[l + j5 * n] +
                                    c6 * w[l + j6 * n] +
                                    c7 * w[l + j7 * n] +
                                    c8 * w[l + j8 * n];
                            }
                        }
                        else if (nj > 8)
                        {
                            j1 = psave[1];
                            j2 = psave[2];
                            j3 = psave[3];
                            j4 = psave[4];
                            j5 = psave[5];
                            j6 = psave[6];
                            j7 = psave[7];
                            j8 = psave[8];
                            c1 = ccs[j1 + s * nps];
                            c2 = ccs[j2 + s * nps];
                            c3 = ccs[j3 + s * nps];
                            c4 = ccs[j4 + s * nps];
                            c5 = ccs[j5 + s * nps];
                            c6 = ccs[j6 + s * nps];
                            c7 = ccs[j7 + s * nps];
                            c8 = ccs[j8 + s * nps];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + j1 * n] +
                                    c2 * w[l + j2 * n] +
                                    c3 * w[l + j3 * n] +
                                    c4 * w[l + j4 * n] +
                                    c5 * w[l + j5 * n] +
                                    c6 * w[l + j6 * n] +
                                    c7 * w[l + j7 * n] +
                                    c8 * w[l + j8 * n];
                            }
                            njbase = 8;
                            njleft = nj - 8;
                            njstep = njleft / 8;
                            njrest = njleft % 8;
                            for (m = 1; m <= njstep; ++m)
                            {
                                j1 = psave[njbase + 1];
                                j2 = psave[njbase + 2];
                                j3 = psave[njbase + 3];
                                j4 = psave[njbase + 4];
                                j5 = psave[njbase + 5];
                                j6 = psave[njbase + 6];
                                j7 = psave[njbase + 7];
                                j8 = psave[njbase + 8];
                                c1 = ccs[j1 + s * nps];
                                c2 = ccs[j2 + s * nps];
                                c3 = ccs[j3 + s * nps];
                                c4 = ccs[j4 + s * nps];
                                c5 = ccs[j5 + s * nps];
                                c6 = ccs[j6 + s * nps];
                                c7 = ccs[j7 + s * nps];
                                c8 = ccs[j8 + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    y[l + rs * n] = y[l + rs * n] +
                                        c1 * w[l + j1 * n] +
                                        c2 * w[l + j2 * n] +
                                        c3 * w[l + j3 * n] +
                                        c4 * w[l + j4 * n] +
                                        c5 * w[l + j5 * n] +
                                        c6 * w[l + j6 * n] +
                                        c7 * w[l + j7 * n] +
                                        c8 * w[l + j8 * n];
                                }
                                njbase += 8;
                            }
                            for (m = 1; m <= njrest; ++m)
                            {
                                j1 = psave[njbase + m];
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
                            if (nj == 1)
                            {
                                ij1 = ppair[1];
                                c1 = ccs[psave[1] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] = c1 * x[l + ij1 * n];
                                }
                            }
                            else if (nj == 2)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                c1 = ccs[psave[1] + s * nps];
                                c2 = ccs[psave[2] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n];
                                }
                            }
                            else if (nj == 3)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                c1 = ccs[psave[1] + s * nps];
                                c2 = ccs[psave[2] + s * nps];
                                c3 = ccs[psave[3] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n];
                                }
                            }
                            else if (nj == 4)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                c1 = ccs[psave[1] + s * nps];
                                c2 = ccs[psave[2] + s * nps];
                                c3 = ccs[psave[3] + s * nps];
                                c4 = ccs[psave[4] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n];
                                }
                            }
                            else if (nj == 5)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                ij5 = ppair[5];
                                c1 = ccs[psave[1] + s * nps];
                                c2 = ccs[psave[2] + s * nps];
                                c3 = ccs[psave[3] + s * nps];
                                c4 = ccs[psave[4] + s * nps];
                                c5 = ccs[psave[5] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n] +
                                        c5 * x[l + ij5 * n];
                                }
                            }
                            else if (nj == 6)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                ij5 = ppair[5];
                                ij6 = ppair[6];
                                c1 = ccs[psave[1] + s * nps];
                                c2 = ccs[psave[2] + s * nps];
                                c3 = ccs[psave[3] + s * nps];
                                c4 = ccs[psave[4] + s * nps];
                                c5 = ccs[psave[5] + s * nps];
                                c6 = ccs[psave[6] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n] +
                                        c5 * x[l + ij5 * n] +
                                        c6 * x[l + ij6 * n];
                                }
                            }
                            else if (nj == 7)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                ij5 = ppair[5];
                                ij6 = ppair[6];
                                ij7 = ppair[7];
                                c1 = ccs[psave[1] + s * nps];
                                c2 = ccs[psave[2] + s * nps];
                                c3 = ccs[psave[3] + s * nps];
                                c4 = ccs[psave[4] + s * nps];
                                c5 = ccs[psave[5] + s * nps];
                                c6 = ccs[psave[6] + s * nps];
                                c7 = ccs[psave[7] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n] +
                                        c5 * x[l + ij5 * n] +
                                        c6 * x[l + ij6 * n] +
                                        c7 * x[l + ij7 * n];
                                }
                            }
                            else if (nj == 8)
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                ij5 = ppair[5];
                                ij6 = ppair[6];
                                ij7 = ppair[7];
                                ij8 = ppair[8];
                                c1 = ccs[psave[1] + s * nps];
                                c2 = ccs[psave[2] + s * nps];
                                c3 = ccs[psave[3] + s * nps];
                                c4 = ccs[psave[4] + s * nps];
                                c5 = ccs[psave[5] + s * nps];
                                c6 = ccs[psave[6] + s * nps];
                                c7 = ccs[psave[7] + s * nps];
                                c8 = ccs[psave[8] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n] +
                                        c5 * x[l + ij5 * n] +
                                        c6 * x[l + ij6 * n] +
                                        c7 * x[l + ij7 * n] +
                                        c8 * x[l + ij8 * n];
                                }
                            }
                            else
                            {
                                ij1 = ppair[1];
                                ij2 = ppair[2];
                                ij3 = ppair[3];
                                ij4 = ppair[4];
                                ij5 = ppair[5];
                                ij6 = ppair[6];
                                ij7 = ppair[7];
                                ij8 = ppair[8];
                                c1 = ccs[psave[1] + s * nps];
                                c2 = ccs[psave[2] + s * nps];
                                c3 = ccs[psave[3] + s * nps];
                                c4 = ccs[psave[4] + s * nps];
                                c5 = ccs[psave[5] + s * nps];
                                c6 = ccs[psave[6] + s * nps];
                                c7 = ccs[psave[7] + s * nps];
                                c8 = ccs[psave[8] + s * nps];
                                for (l = 1; l <= n; l++)
                                {
                                    w[l + i * n] = c1 * x[l + ij1 * n] +
                                        c2 * x[l + ij2 * n] +
                                        c3 * x[l + ij3 * n] +
                                        c4 * x[l + ij4 * n] +
                                        c5 * x[l + ij5 * n] +
                                        c6 * x[l + ij6 * n] +
                                        c7 * x[l + ij7 * n] +
                                        c8 * x[l + ij8 * n];
                                }
                                njbase = 8;
                                njleft = nj - 8;
                                njstep = njleft / 8;
                                njrest = njleft % 8;
                                for (m = 1; m <= njstep; ++m)
                                {
                                    ij1 = ppair[njbase + 1];
                                    ij2 = ppair[njbase + 2];
                                    ij3 = ppair[njbase + 3];
                                    ij4 = ppair[njbase + 4];
                                    ij5 = ppair[njbase + 5];
                                    ij6 = ppair[njbase + 6];
                                    ij7 = ppair[njbase + 7];
                                    ij8 = ppair[njbase + 8];
                                    c1 = ccs[psave[njbase + 1] + s * nps];
                                    c2 = ccs[psave[njbase + 2] + s * nps];
                                    c3 = ccs[psave[njbase + 3] + s * nps];
                                    c4 = ccs[psave[njbase + 4] + s * nps];
                                    c5 = ccs[psave[njbase + 5] + s * nps];
                                    c6 = ccs[psave[njbase + 6] + s * nps];
                                    c7 = ccs[psave[njbase + 7] + s * nps];
                                    c8 = ccs[psave[njbase + 8] + s * nps];
                                    for (l = 1; l <= n; l++)
                                    {
                                        w[l + i * n] = w[l + i * n] +
                                            c1 * x[l + ij1 * n] +
                                            c2 * x[l + ij2 * n] +
                                            c3 * x[l + ij3 * n] +
                                            c4 * x[l + ij4 * n] +
                                            c5 * x[l + ij5 * n] +
                                            c6 * x[l + ij6 * n] +
                                            c7 * x[l + ij7 * n] +
                                            c8 * x[l + ij8 * n];
                                    }
                                    njbase += 8;
                                }
                                for (m = 1; m <= njrest; ++m)
                                {
                                    ij1 = ppair[njbase + m];
                                    c1 = ccs[psave[njbase + m] + s * nps];
                                    for (l = 1; l <= n; l++)
                                    {
                                        w[l + i * n] += c1 * x[l + ij1 * n];
                                    }
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
                        if (ni == 0)
                        {
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = 0.;
                            }
                        }
                        else if (ni == 1)
                        {
                            i1 = psave[1];
                            c1 = ccr[i1 + r * npr];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + i1 * n];
                            }
                        }
                        else if (ni == 2)
                        {
                            i1 = psave[1];
                            i2 = psave[2];
                            c1 = ccr[i1 + r * npr];
                            c2 = ccr[i2 + r * npr];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + i1 * n] +
                                    c2 * w[l + i2 * n];
                            }
                        }
                        else if (ni == 3)
                        {
                            i1 = psave[1];
                            i2 = psave[2];
                            i3 = psave[3];
                            c1 = ccr[i1 + r * npr];
                            c2 = ccr[i2 + r * npr];
                            c3 = ccr[i3 + r * npr];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + i1 * n] +
                                    c2 * w[l + i2 * n] +
                                    c3 * w[l + i3 * n];
                            }
                        }
                        else if (ni == 4)
                        {
                            i1 = psave[1];
                            i2 = psave[2];
                            i3 = psave[3];
                            i4 = psave[4];
                            c1 = ccr[i1 + r * npr];
                            c2 = ccr[i2 + r * npr];
                            c3 = ccr[i3 + r * npr];
                            c4 = ccr[i4 + r * npr];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + i1 * n] +
                                    c2 * w[l + i2 * n] +
                                    c3 * w[l + i3 * n] +
                                    c4 * w[l + i4 * n];
                            }
                        }
                        else if (ni == 5)
                        {
                            i1 = psave[1];
                            i2 = psave[2];
                            i3 = psave[3];
                            i4 = psave[4];
                            i5 = psave[5];
                            c1 = ccr[i1 + r * npr];
                            c2 = ccr[i2 + r * npr];
                            c3 = ccr[i3 + r * npr];
                            c4 = ccr[i4 + r * npr];
                            c5 = ccr[i5 + r * npr];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + i1 * n] +
                                    c2 * w[l + i2 * n] +
                                    c3 * w[l + i3 * n] +
                                    c4 * w[l + i4 * n] +
                                    c5 * w[l + i5 * n];
                            }
                        }
                        else if (ni == 6)
                        {
                            i1 = psave[1];
                            i2 = psave[2];
                            i3 = psave[3];
                            i4 = psave[4];
                            i5 = psave[5];
                            i6 = psave[6];
                            c1 = ccr[i1 + r * npr];
                            c2 = ccr[i2 + r * npr];
                            c3 = ccr[i3 + r * npr];
                            c4 = ccr[i4 + r * npr];
                            c5 = ccr[i5 + r * npr];
                            c6 = ccr[i6 + r * npr];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + i1 * n] +
                                    c2 * w[l + i2 * n] +
                                    c3 * w[l + i3 * n] +
                                    c4 * w[l + i4 * n] +
                                    c5 * w[l + i5 * n] +
                                    c6 * w[l + i6 * n];
                            }
                        }
                        else if (ni == 7)
                        {
                            i1 = psave[1];
                            i2 = psave[2];
                            i3 = psave[3];
                            i4 = psave[4];
                            i5 = psave[5];
                            i6 = psave[6];
                            i7 = psave[7];
                            c1 = ccr[i1 + r * npr];
                            c2 = ccr[i2 + r * npr];
                            c3 = ccr[i3 + r * npr];
                            c4 = ccr[i4 + r * npr];
                            c5 = ccr[i5 + r * npr];
                            c6 = ccr[i6 + r * npr];
                            c7 = ccr[i7 + r * npr];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + i1 * n] +
                                    c2 * w[l + i2 * n] +
                                    c3 * w[l + i3 * n] +
                                    c4 * w[l + i4 * n] +
                                    c5 * w[l + i5 * n] +
                                    c6 * w[l + i6 * n] +
                                    c7 * w[l + i7 * n];
                            }
                        }
                        else if (ni == 8)
                        {
                            i1 = psave[1];
                            i2 = psave[2];
                            i3 = psave[3];
                            i4 = psave[4];
                            i5 = psave[5];
                            i6 = psave[6];
                            i7 = psave[7];
                            i8 = psave[8];
                            c1 = ccr[i1 + r * npr];
                            c2 = ccr[i2 + r * npr];
                            c3 = ccr[i3 + r * npr];
                            c4 = ccr[i4 + r * npr];
                            c5 = ccr[i5 + r * npr];
                            c6 = ccr[i6 + r * npr];
                            c7 = ccr[i7 + r * npr];
                            c8 = ccr[i8 + r * npr];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + i1 * n] +
                                    c2 * w[l + i2 * n] +
                                    c3 * w[l + i3 * n] +
                                    c4 * w[l + i4 * n] +
                                    c5 * w[l + i5 * n] +
                                    c6 * w[l + i6 * n] +
                                    c7 * w[l + i7 * n] +
                                    c8 * w[l + i8 * n];
                            }
                        }
                        else if (ni > 8)
                        {
                            i1 = psave[1];
                            i2 = psave[2];
                            i3 = psave[3];
                            i4 = psave[4];
                            i5 = psave[5];
                            i6 = psave[6];
                            i7 = psave[7];
                            i8 = psave[8];
                            c1 = ccr[i1 + r * npr];
                            c2 = ccr[i2 + r * npr];
                            c3 = ccr[i3 + r * npr];
                            c4 = ccr[i4 + r * npr];
                            c5 = ccr[i5 + r * npr];
                            c6 = ccr[i6 + r * npr];
                            c7 = ccr[i7 + r * npr];
                            c8 = ccr[i8 + r * npr];
                            for (l = 1; l <= n; l++)
                            {
                                y[l + rs * n] = c1 * w[l + i1 * n] +
                                    c2 * w[l + i2 * n] +
                                    c3 * w[l + i3 * n] +
                                    c4 * w[l + i4 * n] +
                                    c5 * w[l + i5 * n] +
                                    c6 * w[l + i6 * n] +
                                    c7 * w[l + i7 * n] +
                                    c8 * w[l + i8 * n];
                            }
                            nibase = 8;
                            nileft = ni - 8;
                            nistep = nileft / 8;
                            nirest = nileft % 8;
                            for (m = 1; m <= nistep; ++m)
                            {
                                i1 = psave[nibase + 1];
                                i2 = psave[nibase + 2];
                                i3 = psave[nibase + 3];
                                i4 = psave[nibase + 4];
                                i5 = psave[nibase + 5];
                                i6 = psave[nibase + 6];
                                i7 = psave[nibase + 7];
                                i8 = psave[nibase + 8];
                                c1 = ccr[i1 + r * npr];
                                c2 = ccr[i2 + r * npr];
                                c3 = ccr[i3 + r * npr];
                                c4 = ccr[i4 + r * npr];
                                c5 = ccr[i5 + r * npr];
                                c6 = ccr[i6 + r * npr];
                                c7 = ccr[i7 + r * npr];
                                c8 = ccr[i8 + r * npr];
                                for (l = 1; l <= n; l++)
                                {
                                    y[l + rs * n] = y[l + rs * n] +
                                        c1 * w[l + i1 * n] +
                                        c2 * w[l + i2 * n] +
                                        c3 * w[l + i3 * n] +
                                        c4 * w[l + i4 * n] +
                                        c5 * w[l + i5 * n] +
                                        c6 * w[l + i6 * n] +
                                        c7 * w[l + i7 * n] +
                                        c8 * w[l + i8 * n];
                                }
                                nibase += 8;
                            }
                            for (m = 1; m <= nirest; m++)
                            {
                                i1 = psave[nibase + m];
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
    }
    
    return 0;
}


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


/*  AUTHOR      : Norbert Flocke */
/* ------------------------------------------------------------------------ */
int erd__ctr_1st_half_ (int * n, int * npmax, int * npmin,
                        int * mij, int * nrs, int * nblock,
                        int * ncr, int * ncs, int * npr,
                        int * nps, double * ccr, double * ccs,
                        int * ccbegr, int * ccbegs, int * ccendr,
                        int * ccends, int * primr, int * prims,
                        int * equalrs, int * swaprs, int * pused,
                        int * psave, int * ppair, double * x,
                        double * w, double * y)
{
    erd__ctr_1st_half (*n, *npmax, *npmin, *mij, *nrs, *nblock,
                       *ncr, *ncs, *npr, *nps,
                       ccr, ccs, ccbegr, ccbegs, ccendr, ccends,
                       primr, prims, *equalrs, *swaprs,
                       pused, psave, ppair,
                       x, w, y);
    
    return 0;
}