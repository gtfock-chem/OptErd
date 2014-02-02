#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "erd.h"


#pragma offload_attribute(push, target(mic))

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__HRR_STEP */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation performs a single horizontal recurrence */
/*                relation (HRR) step on the input transformation matrix */
/*                WIN of formal dimension NXYZET x NAB, where NXYZET is */
/*                the dimension of the starting shell combination (e0|. */
/*                The columns of the matrix WIN are combined such that: */
/*                   (a(b+1i)| = ((a+1i)b| + (Ai-Bi)(ab| ; i=x,y,z */
/*                Note, that the shell symbol a stands for a range of */
/*                shells, starting from shell x and ending at shell p-1 */
/*                for (a(b+1i)| and (ab| and ending at shell p for */
/*                ((a+1i)b|. */
/*                Since the input transformation matrix WIN contains */
/*                already all the HRR info from previous steps, the */
/*                output transformation matrix will have the complete */
/*                information for performing a HRR of the following */
/*                kind: */
/*                           (e0| --> (a(b+1i)| */
/*                Strategy used to perform the HRR step: */
/*                ------------------------------------- */
/*                The HRR is split into two parts: part I) deals with */
/*                the (Ai-Bi)(ab| part and part II) with the ((a+1i)b| */
/*                part. */
/*                Part I):  In this case we have the following scheme, */
/*                          with a-shell range a = x to p-1: */
/*                           W = WIN (xa,ya,za,xb,yb,zb) -> */
/*                               + ABX*W to WOUT (xa,ya,za,xb+1,yb,zb) */
/*                               + ABY*W to WOUT (xa,ya,za,xb,yb+1,zb) */
/*                               + ABZ*W to WOUT (xa,ya,za,xb,yb,zb+1) */
/*                Part II): Here we have the following scheme, with */
/*                          a-shell range a = x+1 to p: */
/*                           W = WIN (xa,ya,za,xb,yb,zb) -> */
/*                               + W to WOUT (xa-1,ya,za,xb+1,yb,zb) */
/*                               + W to WOUT (xa,ya-1,za,xb,yb+1,zb) */
/*                               + W to WOUT (xa,ya,za-1,xb,yb,zb+1) */
/*                How to get from the b- to the (b+1)-shell monomials: */
/*                --------------------------------------------------- */
/*                To perform parts I) and II) of the HRR we need an */
/*                algorithm which creates a unique set of (b+1)-shell */
/*                monomials from those of the b-shell. The following */
/*                strategy is adopted: */
/*                     1) if x-exponent in b-shell monomial is > 0, */
/*                        add +1 to x-exponent. */
/*                     2) if x-exponent in b-shell monomial is = 0 */
/*                        and y-exponent is > 0, add +1 to x-exponent */
/*                        and y-exponent. */
/*                     3) if x-exponent and y-exponent in b-shell */
/*                        monomial is = 0, add +1 to all exponents. */
/*                Example for b=2 --> b+1=3 case: */
/*                            200 --> 300 */
/*                            110 --> 210 */
/*                            101 --> 201 */
/*                            020 --> 120,030 */
/*                            011 --> 111,021 */
/*                            002 --> 102,012,003 */
/*                  Input: */
/*                    NAB          =  total # of monomials of the */
/*                                    input ((a+1i)b| part, i.e. # */
/*                                    of columns of input transformation */
/*                                    matrix */
/*                    NABO         =  total # of monomials of the */
/*                                    output (a(b+1i)| part, i.e. # */
/*                                    of columns of output transformation */
/*                                    matrix */
/*                    MROWIN       =  maximum # of nonzero row elements */
/*                                    per column in input transformation */
/*                                    matrix */
/*                    MROWOUT      =  maximum # of nonzero row elements */
/*                                    per column in output transformation */
/*                                    matrix */
/*                    NXYZX        =  total # of monomials for shell x */
/*                    NXYZP        =  total # of monomials for shell p */
/*                    NXYZA        =  sum of total # of monomials for */
/*                                    shell range a = x to p */
/*                    NXYZB        =  total # of monomials for shell b */
/*                    NXYZAO       =  sum of total # of monomials for */
/*                                    shell range a = x to p-1 */
/*                    SHELLx       =  shell type for shells x = X,P,B */
/*                    ABm          =  the m=x,y,z-coordinate differences */
/*                                    between centers A and B */
/*                    CPAIR        =  int array that will hold the */
/*                                    pair of column indices of the input */
/*                                    transformation matrix that have */
/*                                    to be combined to form each output */
/*                                    transformation matrix column */
/*                    NROWIN (J)   =  # of nonzero row elements of J-th */
/*                                    input transformation matrix column */
/*                    ROWIN (I,J)  =  nonzero row indices of J-th input */
/*                                    transformation matrix column */
/*                    WIN (I,J)    =  nonzero elements of the input */
/*                                    transformation matrix corresponding */
/*                                    to nonzero row index I and column J */
/*                  Output: */
/*                    NROWOUT (J)  =  # of nonzero row elements of J-th */
/*                                    output transformation matrix column */
/*                    ROWOUT (I,J) =  nonzero row indices of J-th output */
/*                                    transformation matrix column */
/*                    WOUT (I,J)   =  nonzero elements of the output */
/*                                    transformation matrix corresponding */
/*                                    to nonzero row index I and column J */
/* ------------------------------------------------------------------------ */
int erd__hrr_step (int nab, int nabo, int mrowin,
                   int mrowout, int nxyzx, int nxyzp,
                   int nxyza, int nxyzb, int nxyzao,
                   int shellx, int shellp, int shellb,
                   double abx, double aby, double abz,
                   int *cpair, int *nrowin, int *rowin, double *win,
                   int *nrowout, int *rowout, double *wout)
{
    int cpair_offset, rowin_offset, rowout_offset,
        win_offset, wout_offset;

    int b, i, j, m, n, x, y, z, c1, c2, i1, i2, xa, ya, xo, yo,
        zo, row1, row2, nout, nrow1, nrow2;
    double coeff;
    int rdiff, offya, nxbgt0, shella, offybo;

    --nrowin;
    --nrowout;
    cpair_offset = 1 + nabo * 1;
    cpair -= cpair_offset;
    win_offset = 1 + mrowin * 1;
    win -= win_offset;
    rowin_offset = 1 + mrowin * 1;
    rowin -= rowin_offset;
    wout_offset = 1 + mrowout * 1;
    wout -= wout_offset;
    rowout_offset = 1 + mrowout * 1;
    rowout -= rowout_offset;

    nxbgt0 = nxyzb - shellb - 1;
    offybo = nxyzao * (shellb + 1);
/*             ...determine the column pairs to be 'added' to define */
/*                the output tranformation matrix. */
/*                      ------- Part I) section --------- */
/*                Outer loop over b- to (b+1)-shell monomials with */
/*                xb-exponent > 0. Inner loop over a-part contractions */
/*                for part I) of HRR for a-shell range a = x to p-1. */
    x = -nxyza;
    xo = -nxyzao;
    for (b = 1; b <= nxbgt0; b++)
    {
        x += nxyza;
        xo += nxyzao;
        if (fabs(abx) > 0.0)
        {
            for (j = 1; j <= nxyzao; ++j)
            {
                cpair[xo + j + nabo] = x + j;
            }
        }
    }
/*             ...outer loop over b- to (b+1)-shell monomials with */
/*                xb-exponent = 0 and yb-exponent > 0. Inner loop over */
/*                a-part contractions for part I) of HRR for a-shell */
/*                range a = x to p-1. */
    yo = xo + offybo;
    for (b = 1; b <= shellb; ++b)
    {
        x += nxyza;
        xo += nxyzao;
        yo += nxyzao;
        if (fabs(abx) > 0.)
        {
            for (j = 1; j <= nxyzao; ++j)
            {
                cpair[xo + j + nabo] = x + j;
            }
        }
        if (fabs (aby) > 0.)
        {
            for (j = 1; j <= nxyzao; ++j)
            {
                cpair[yo + j + nabo] = x + j;
            }
        }
    }

/*             ...last b- to (b+1)-shell monomial with xb-exponent = 0 */
/*                and yb-exponent = 0. Inner loop over a-part */
/*                contractions for part I) of HRR for a-shell range */
/*                a = x to p-1. */
    x += nxyza;
    xo += nxyzao;
    yo += nxyzao;
    zo = yo + nxyzao;
    if (fabs(abx) > 0.0)
    {
        for (j = 1; j <= nxyzao; ++j)
        {
            cpair[xo + j + nabo] = x + j;
        }
    }
    if (fabs (aby) > 0.)
    {
        for (j = 1; j <= nxyzao; ++j)
        {
            cpair[yo + j + nabo] = x + j;
        }
    }
    if (fabs (abz) > 0.)
    {
        for (j = 1; j <= nxyzao; ++j)
        {
            cpair[zo + j + nabo] = x + j;
        }
    }


/*                      ------- Part II) section --------- */
/*             ...outer loop over b- to (b+1)-shell monomials with */
/*                xb-exponent > 0. Inner loop over a-part contractions */
/*                for part II) of HRR for a-shell range a = x+1 to p. */
    x = 0;
    xo = 0;
    for (b = 1; b <= nxbgt0; ++b)
    {
        x += nxyzx;
        for (shella = shellx; shella <= shellp - 1; ++shella)
        {
            for (xa = shella; xa >= 0; --xa)
            {
                for (ya = shella - xa; ya >= 0; --ya)
                {
                    ++x;
                    ++xo;
                    cpair[xo + nabo * 2] = x;
                }
            }
            x = x + shella + 2;
        }
    }


/*             ...outer loop over b- to (b+1)-shell monomials with */
/*                xb-exponent = 0 and yb-exponent > 0. Inner loop over */
/*                a-part contractions for part II) of HRR for a-shell */
/*                range a = x+1 to p. */
    yo = xo + offybo;
    for (b = 1; b <= shellb; ++b)
    {
        x += nxyzx;
        for (shella = shellx; shella <= shellp - 1; ++shella)
        {
            offya = 0;
            for (xa = shella; xa >= 0; --xa)
            {
                ++offya;
                for (ya = offya - 1; ya >= 0; --ya)
                {
                    ++x;
                    y = x + offya;
                    ++xo;
                    ++yo;
                    cpair[xo + nabo * 2] = x;
                    cpair[yo + nabo * 2] = y;
                }
            }
            x = y + 1;
        }
    }

/*             ...last b- to (b+1)-shell monomial with xb-exponent = 0 */
/*                and yb-exponent = 0. Inner loop over a-part */
/*                contractions for part II) of HRR for a-shell range */
/*                a = x+1 to p. */
    x += nxyzx;
    zo = yo + nxyzao;
    for (shella = shellx; shella <= shellp - 1; ++shella)
    {
        offya = 0;
        for (xa = shella; xa >= 0; --xa)
        {
            ++offya;
            for (ya = offya - 1; ya >= 0; --ya)
            {
                ++x;
                y = x + offya;
                z = y + 1;
                ++xo;
                ++yo;
                ++zo;
                cpair[xo + nabo * 2] = x;
                cpair[yo + nabo * 2] = y;
                cpair[zo + nabo * 2] = z;
            }
        }
        x = z;
    }


/*             ...the column pairs are ready. Construct the new */
/*                transformation matrix. */
    for (n = 1; n <= nabo; ++n)
    {
        if (n <= xo)
        {
            coeff = abx;
        }
        else if (n <= yo)
        {
            coeff = aby;
        }
        else
        {
            coeff = abz;
        }
        if (fabs(coeff) > 0.)
        {
            c1 = cpair[n + nabo];
            c2 = cpair[n + nabo * 2];
            nrow1 = nrowin[c1];
            nrow2 = nrowin[c2];
            m = MIN(nrow1, nrow2);
            i1 = 1;
            i2 = 1;
            nout = 0;
            for (i = 2; i <= m; ++i)
            {
                row1 = rowin[i1 + c1 * mrowin];
                row2 = rowin[i2 + c2 * mrowin];
                rdiff = row1 - row2;
                if (rdiff == 0)
                {
                    ++nout;
                    wout[nout + n * mrowout] =
                        coeff * win[i1 + c1 * mrowin] +
                        win[i2 + c2 * mrowin];
                    rowout[nout + n * mrowout] = row1;
                    ++i1;
                    ++i2;
                }
                else if (rdiff < 0)
                {
                    ++nout;
                    wout[nout + n * mrowout] = 
                        coeff * win[i1 + c1 * mrowin];
                    rowout[nout + n * mrowout] = row1;
                    ++i1;
                }
                else
                {
                    ++nout;
                    wout[nout + n * mrowout] = win[i2 + c2 * mrowin];
                    rowout[nout + n * mrowout] = row2;
                    ++i2;
                }
            }
          L3000:
            row1 = rowin[i1 + c1 * mrowin];
            row2 = rowin[i2 + c2 * mrowin];
            rdiff = row1 - row2;
            if (rdiff == 0)
            {
                ++nout;
                wout[nout + n * mrowout] =
                    coeff * win[i1 + c1 * mrowin]
                    + win[i2 + c2 * mrowin];
                rowout[nout + n * mrowout] = row1;
                ++i1;
                ++i2;
            }
            else if (rdiff < 0)
            {
                ++nout;
                wout[nout + n * mrowout] =
                    coeff * win[i1 + c1 * mrowin];
                rowout[nout + n * mrowout] = row1;
                ++i1;
            }
            else
            {
                ++nout;
                wout[nout + n * mrowout] = win[i2 + c2 * mrowin];
                rowout[nout + n * mrowout] = row2;
                ++i2;
            }
            if (i1 > nrow1)
            {
                for (i = i2; i <= nrow2; ++i)
                {
                    row2 = rowin[i + c2 * mrowin];
                    ++nout;
                    wout[nout + n * mrowout] = win[i + c2 * mrowin];
                    rowout[nout + n * mrowout] = row2;
                }
            }
            else if (i2 > nrow2)
            {
                for (i = i1; i <= nrow1; ++i)
                {
                    row1 = rowin[i + c1 * mrowin];
                    ++nout;
                    wout[nout + n * mrowout] =
                        coeff * win[i + c1 * mrowin];
                    rowout[nout + n * mrowout] = row1;
                }
            }
            else
            {
                goto L3000;
            }
            nrowout[n] = nout;
        }
        else
        {
            c2 = cpair[n + nabo * 2];
            nrow2 = nrowin[c2];
            for (i2 = 1; i2 <= nrow2; ++i2)
            {
                wout[i2 + n * mrowout] = win[i2 + c2 * mrowin];
                rowout[i2 + n * mrowout] = rowin[i2 + c2 * mrowin];
            }
            nrowout[n] = nrow2;
        }
    }


    return 0;
}

#pragma offload_attribute(pop)
