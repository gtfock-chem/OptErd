#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "erd.h"
#include "jacobi.h"


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

static double d_sign(double *a, double *b)
{
    double x;
    x = (*a >= 0 ? *a : - *a);
    return( *b >= 0 ? x : -x);
}


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__RYS_X_ROOTS_WEIGHTS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : Master routine to evaluate roots and weights in the */
/*                interval [0,1] over the Rys weight function: */
/*                                       exp(-T*x) */
/*                             W   (x) = --------- */
/*                              Rys      2*sqrt(x) */
/*                using the general Gaussian Quadrature technique */
/*                consisting of the following basic steps: */
/*                   1) Calculate the auxillary polynomial moments */
/*                   2) Calculate the auxillary polynomial coefficients */
/*                   3) Set up the tridiagonal symmetric terminal */
/*                      matrix */
/*                   4) Solve the tridiagonal symmetric terminal */
/*                      matrix for roots and weights. */
/*                In this routine all T-values are treated all at once. */
/*                The value of the T-parameter dictates which type of */
/*                auxillary set of polynomials is to be used for the */
/*                modified Chebyshev algorithm. The auxillary polynomials */
/*                are often chosen (as is the case here) to be orthogonal */
/*                relative to some classical weight function. */
/*                The classical weights used to establish the auxillary */
/*                polynomials are the following (notation according to */
/*                M.Abramowitz and I.A.Stegun, Handbook of Mathematical */
/*                Functions, 1964): */
/*                i) Range of validity: 0 =< T =< 30 */
/*                   Here we use shifted Jacobi weights and polynomials: */
/*                        W      (p,q,x)  =  (1-x)^(p-q) * x^(q-1) */
/*                         Jacobi */
/*                        i-th Jacobi polynomial  =  G (p,q,x) */
/*                                                    i */
/*                   with conditions: x in interval [0,1] */
/*                                    p-q > -1 , q > 0 */
/*                                    i = 0,1,2,... */
/*                ii) Range of validity: 1 < T =< oo (infinity) */
/*                   Here we use generalized Laguerre weights and */
/*                   polynomials: */
/*                        W        (a,x)  =  exp^(-x) * x^a */
/*                         Laguerre */
/*                        i-th Laguerre polynomial  =  L (a,x) */
/*                                                      i */
/*                   with conditions: x in interval [0,inf] */
/*                                    a > -1 */
/*                                    i = 0,1,2,... */
/*                Range of validity means that for all the specified */
/*                T's within the range the resulting moment integrals */
/*                are accurate to within 1.D-16. */
/*                  Input: */
/*                    NT           =  # of T-values */
/*                    NTGQP        =  # of roots times # of T-values */
/*                    NGQP         =  # of gaussian quadrature points */
/*                                    (roots) */
/*                    NMOM         =  # of necessary moment integrals */
/*                                    to calculate the quadrature roots */
/*                    TVAL         =  the T-values */
/*                    RYSZERO      =  the zeroth Rys moments for all */
/*                                    T-values */
/*                    A,B          =  will contain the recurrence */
/*                                    coefficients for the auxillary */
/*                                    polynomials */
/*                    MOM          =  will contain the normed auxillary */
/*                                    polynomial modified moments */
/*                    DIA,OFF      =  will contain the diagonal and */
/*                                    offdiagonal elements of the */
/*                                    tridiagonal symmetric terminal */
/*                                    matrix */
/*                    ROW1,ROW2    =  will be used to evaluate the */
/*                                    tridiagonal elements of the */
/*                                    symmetric terminal matrix in an */
/*                                    efficient way using Sack and */
/*                                    Donovan's method */
/*                  Output: */
/*                    RTS          =  the roots array */
/*                    WTS          =  the weights array */
/* ------------------------------------------------------------------------ */
int erd__rys_x_roots_weights_ (int * nt, int * ngqp, int * nmom, double * tval,
                               double * ryszero, double * a,
                               double * b, double * mom,
                               double * dia, double * off,
                               double * row1, double * row2,
                               double * rts, double * wts)
{
    int i__1, i__2, i__3;
    double d__1, d__2;
    double c__, d__, f, g;
    int i__, j, m, n;
    double p, r__, s, t, r1;
    int ip1, jp1;
    double lim1, lim2, lim3, binc, sinc;
    int imax, jmax;
    double momi;
    int iter;
    double tinv, texp, zmom, root;
    int nrts;
    double zinv, test1, test2, tinv2, scale, sigma, theta, momim1,
        momip1, tinvhf, tpower, tinvsq, momzero;

    --ryszero;
    --tval;
    --wts;
    --rts;
    --off;
    --dia;
    --row2;
    --row1;
    --mom;
    b -= 2;
    --a;

    nrts = 0;
    i__1 = *nt;
    for (n = 1; n <= i__1; ++n)
    {
        t = tval[n];
        momzero = ryszero[n];
        if (t <= 15.)
        {
/*             ...The Jacobi section. Check first, if the number of */
/*                moments wanted is within limits for the Jacobi case. */
/*                If ok, we proceed by calculating two steps at the */
/*                same time: */
/*                1) normed shifted Jacobi modified moments defined as */
/*                   integrals on [0,1] over the Rys weight function: */
/*                                        exp(-Tx) */
/*                             W   (x) = --------- */
/*                              Rys      2*sqrt(x) */
/*                   and shifted Jacobi polynomials: */
/*                                G (0.5,0.5,x) */
/*                                 i */
/*                   A 3-term recursion relation can be given for these */
/*                   moments as follows: */
/*                       MOM (i+1) = R * MOM (i)  +  S * MOM (i-1) */
/*                   where */
/*                      R  =  (2i+1)/2T  +  (2i+1)/(4i-1)(4i+3) */
/*                      S  =  2i(2i+1)(2i-1)**2 / (4i-3)(4i+1)(4i-1)**2 */
/*                   However, in order to evaluate the moments to the */
/*                   required accuracy, a downward recursion must be */
/*                   applied starting with some tiny seed values for */
/*                   MOM (i+1) and MOM (i): */
/*                       MOM (i-1) = 1/S * ( MOM (i+1) - R * MOM (i)) */
/*                   The sequence of moments thus obtained is not */
/*                   normalized. Hence the downward recursion is */
/*                   performed down to the zeroth moment ZMOM and the */
/*                   relevant moments are normalized by subsequent */
/*                   division with ZMOM. */
/*                   The above downward recursion technique is extremely */
/*                   unstable if T is very small, in which case even the */
/*                   small seed values are not small enough to prevent */
/*                   overflow of the lowest not normalized moments. */
/*                   In such a case it is better to build the normed */
/*                   sequence starting with the first moment after the */
/*                   following considerations. An expansion of MOM (i) */
/*                   in terms of powers of T shows that the leading */
/*                   term is in the i-th power of T, i.e. */
/*                                    inf           k */
/*                       MOM (i)  =   sum  c (i,k) T */
/*                                   k = i */
/*                   with leading coefficient c (i,i) getting smaller */
/*                   with increasing i. Hence for very small T it is */
/*                   sufficient to set: */
/*                                             i */
/*                       MOM (i)  =   c (i,i) T     (for T very small) */
/*                   an equation which is valid in terms of computer */
/*                   accuray if T is less or equal to the minimum */
/*                   possible nonzero number within the precision of */
/*                   the mantissa allowed (here in the present code */
/*                   this is double precision, hence T less or equal */
/*                   1.D-16). Rather than evaluating each leading */
/*                   coefficient c (i,i) individualy (they are given */
/*                   by very! complicated expressions in terms of sums */
/*                   of double factorials), the procedure adopted here */
/*                   was to predetermine each c (i,i) numerically by */
/*                   running the routine with T = 1.D-16 and presetting */
/*                   the resulting coefficients in a data array CSMALL */
/*                   to be used whenever T =< 1.D-16. */
/*                   If T is in the range 1.D-16 to TMAX = 30.D0, then */
/*                   the routine was calibrated such that the mantissa */
/*                   of the resulting moments are accurate to the 16th */
/*                   decimal place. Calibration in this context means the */
/*                   predetermination of the maximum moment (controled */
/*                   in the code by its index IMAX) that has to be */
/*                   generated to perform the downward recurrence */
/*                   relation starting with the tiny seeds. Hence the */
/*                   calibration depends on four things: */
/*                       i) The total maximum number MOMMAX of normed */
/*                          moments wanted. */
/*                      ii) The mantissa accuray wanted for the normed */
/*                          moments. */
/*                     iii) The tiny values of the two seeds (the nonzero */
/*                          seed should correpond to the smallest */
/*                          possible nonzero number representable on */
/*                          the present computer). */
/*                      iv) T range. Each T range requires a different */
/*                          maximum moment (i.e. a different IMAX value) */
/*                          to start with. */
/*                   In order to perform a calibration for a different */
/*                   machine, one has to write a separate small program */
/*                   using this routine with different seeds, T- and IMAX */
/*                   values. Once the seeds have been set, this little */
/*                   program should recalculate the MOMMAX normed */
/*                   moments for increasing IMAX values starting with */
/*                   IMAX= MOMMAX. The IMAX value finally taken for a */
/*                   particular T-value should then obey the following */
/*                   inequality: */
/*                   | MOM (i,IMAX) - MOM (i,IMAX-1)| */
/*                   | -----------------------------| < mantissa accuracy */
/*                   |         MOM (i,IMAX-1)       | */
/*                       COMMENTS FOR FAST EVALUATION OF THE MOMENTS */
/*                      --------------------------------------------- */
/*                   As can be seen from the downward recursion formula */
/*                   for the moments, all we need are the values of the R */
/*                   and 1/S recursion parameters, which are conveniently */
/*                   precomputed (with R being split as R = R1/2T + R2) */
/*                   and supplied via an include table 'erd__jacobi.inc'. */
/*                   Hence 'erd__jacobi.inc' will have the values of the */
/*                   complicated R2 and SINV = 1/S expressions in the */
/*                   range from i = 1 to 100. The R1 values are simply */
/*                   calculated by using the appropriate decrement value */
/*                   of 2. The include file 'erd__jacobi.inc' will also */
/*                   contain the CSMALL array for the moments of very */
/*                   small T's. */
/*                2) recurrence coefficients for the shifted Jacobi */
/*                   polynomials G (0.5,0.5,x), denoted simply by G (x): */
/*                         G (x)  =  1 */
/*                          0 */
/*                         G (x)  =  (x - A ) */
/*                          1              1 */
/*                         G   (x)  =  (x - A   ) * G (x) - B   G   (x) */
/*                          i+1              i+1     i       i+1 i-1 */
/*                   The result consists of the recurrence coefficients: */
/*                               A (I) , I = 1,NMOM */
/*                               B (I) , I = 2,NMOM */
/*                   whose values are given by the following expressions: */
/*                        A (i+1) = 4i*(2i+1)-1 / (4i+3)(4i-1) */
/*                        B (i+1) = (4i**2)((4i**2)-4i+1) / */
/*                                  (4i-3)(4i+1)((4i-1)**2) */
/*                   Since these are complicated and time consuming */
/*                   expressions, they are precalculated and included */
/*                   into the include file 'erd__jacobi.inc'. */
            assert (*nmom <= 30);


/*             ...the very small T case. */
            if (t <= 1e-16)
            {
                imax = MIN (*nmom, 16);
                a[1] = ajac[0];
                mom[1] = csmall[0] * t;
                tpower = t;
                i__2 = imax;
                for (i__ = 2; i__ <= i__2; ++i__)
                {
                    tpower *= t;
                    a[i__] = ajac[i__ - 1];
                    b[i__] = bjac[i__ - 2];
                    mom[i__] = csmall[i__ - 1] * tpower;
                }
                i__2 = *nmom;
                for (i__ = imax + 1; i__ <= i__2; ++i__)
                {
                    a[i__] = ajac[i__ - 1];
                    b[i__] = bjac[i__ - 2];
                    mom[i__] = 0.;
                }
            }
            else
            {
/*             ...the general Jacobi case. Set maximum number of moments */
/*                necessary to get required moments to an accuracy of */
/*                at least 1.D-16. See above in the routine description */
/*                for details and calibration of the setting. */
                if (*nmom <= 5)
                {
                    if (t < 1e-6)
                    {
                        imax = *nmom + 1;
                    }
                    else if (t < .1)
                    {
                        imax = *nmom + 3;
                    }
                    else if (t < 2.)
                    {
                        imax = *nmom + 7;
                    }
                    else if (t < 10.)
                    {
                        imax = *nmom + 13;
                    }
                    else
                    {
                        imax = *nmom + 22;
                    }
                }
                else
                {
                    if (t < 1e-6)
                    {
                        imax = *nmom;
                    }
                    else if (t < .1)
                    {
                        imax = *nmom + 2;
                    }
                    else if (t < 2.)
                    {
                        imax = *nmom + 4;
                    }
                    else if (t < 10.)
                    {
                        imax = *nmom + 8;
                    }
                    else
                    {
                        imax = *nmom + 16;
                    }
                }

/*             ...proceed by setting seed values for downward recursion */
/*                to obtain minimal solution and start recursive */
/*                evaluation for all Jacobi moments down to first moment */
/*                (two loops are necessary here due to calculation of */
/*                higher moments than actually needed later on). */
                momi = 1e-300;
                momip1 = 0.;
                tinvhf = .5 / t;
                r1 = (double) ((imax << 1) + 5);
                i__2 = *nmom + 2;
                for (i__ = imax + 1; i__ >= i__2; --i__)
                {
                    r1 += -2.;
                    r__ = r1 * tinvhf + r2[i__ - 1];
                    momim1 = sinv[i__ - 1] * (momip1 - r__ * momi);
                    momip1 = momi;
                    momi = momim1;
                }
                for (i__ = *nmom + 1; i__ >= 2; --i__)
                {
                    r1 += -2.;
                    r__ = r1 * tinvhf + r2[i__ - 1];
                    momim1 = sinv[i__ - 1] * (momip1 - r__ * momi);
                    mom[i__ - 1] = momim1;
                    momip1 = momi;
                    momi = momim1;
                }

/*             ...evaluate zeroth moment and normalize sequence. */
/*                If the absolute zeroth moment is less than the */
/*                approximate absolute nonzero minimum (set here */
/*                equal to 1.D-300), the normalization looses its */
/*                meaning and the calculation must be stopped. */
/*                Set also here the recurrence relation coefficients */
/*                A and B for the shifted Jacobi polynomials. */
                r__ = tinvhf * 3. + r2[0];
                zmom = sinv[0] * (momip1 - r__ * momi);
                assert (fabs (zmom) >= 1e-300);                
                a[1] = ajac[0];
                zinv = 1. / zmom;
                mom[1] *= zinv;
                i__2 = *nmom;
                for (i__ = 2; i__ <= i__2; ++i__)
                {
                    a[i__] = ajac[i__ - 1];
                    b[i__] = bjac[i__ - 2];
                    mom[i__] *= zinv;
                }
            }
        }
        else
        {
/*             ...The Laguerre section. As for the Jacobi case, */
/*                we calculate two steps at the same time: */
/*                1) normed Laguerre modified moments defined as */
/*                   integrals on [0,1] over the Rys weight function: */
/*                                        exp(-Tx) */
/*                             W   (x) = --------- */
/*                              Rys      2*sqrt(x) */
/*                   and monic generalized 'scaled' Laguerre polynomials: */
/*                                L (-0.5,Tx) */
/*                                 i */
/*                   A closed formula can be given for these moments in */
/*                   terms of monic generalized 'scaled' Laguerre */
/*                   polynomials with generalization +0.5: */
/*                       MOM (i) = SCALE * L   (+0.5,T) */
/*                                          i-1 */
/*                   where */
/*                               - exp(-T) */
/*                       SCALE = ---------   ;  F0(T) = Rys zero moment */
/*                                2T*F0(T) */
/*                   The recursion relation for the +0.5 polynomials is: */
/*                     L   (+0.5,T) = R * L   (+0.5,T) - S * L   (+0.5,T) */
/*                      i-1                i-2                i-3 */
/*                   where */
/*                          R  =  (T - 2i + 5/2) / T */
/*                          S  =  (i - 2)(i - 3/2) / T*T */
/*                   All moments MOM (i);i=1,2,...NMOM using the above */
/*                   outlined algorithm are evaluated. */
/*                2) recurrence coefficients for T-scaled generalized */
/*                   monic Laguerre polynomials L (-0.5,Tx), denoted */
/*                   simply by L (Tx): */
/*                       L (Tx)    =   1 */
/*                        0 */
/*                       L (Tx)    =   (x - A ) */
/*                        1                  1 */
/*                       L   (Tx)  =   (x - A   ) * L (Tx) - B   L   (Tx) */
/*                        i+1                i+1     i        i+1 i-1 */
/*                   The result consists of the recurrence coefficients: */
/*                             A (I) , I = 1,NMOM */
/*                             B (I) , I = 2,NMOM */
/*                   whose values are given by the following expressions: */
/*                           A (1) = 1 / 2T */
/*                         A (i+1) = (2i+1/2) / T */
/*                         B (i+1) = i*(i-1/2) / (T*T) */
            texp = exp (-t);
            tinv = 1. / t;
            tinv2 = tinv * 2.;
            tinvhf = tinv * .5;
            tinvsq = tinv * tinv;
            scale = -tinvhf * texp / momzero;
            if (*nmom == 1)
            {
                a[1] = tinvhf;
                mom[1] = scale;
            }
            else
            {
                a[1] = tinvhf;
                a[2] = tinvhf + tinv2;
                b[2] = tinvsq * .5;
                mom[1] = scale;
                r__ = 1. - tinv * 1.5;
                mom[2] = scale * r__;
                s = 0.;
                binc = .5;
                sinc = -.5;
                lim2 = r__;
                lim3 = 1.;
                i__2 = *nmom;
                for (i__ = 3; i__ <= i__2; ++i__)
                {
                    binc += 2.;
                    a[i__] = a[i__ - 1] + tinv2;
                    b[i__] = b[i__ - 1] + binc * tinvsq;
                    sinc += 2.;
                    r__ -= tinv2;
                    s += sinc * tinvsq;
                    lim1 = r__ * lim2 - s * lim3;
                    mom[i__] = scale * lim1;
                    lim3 = lim2;
                    lim2 = lim1;
                }
            }
        }

/*             ...This section calculates a symmetric terminal matrix */
/*                using the normed modified moments MOM and related */
/*                recurrence coefficients A and B of monic polynomials */
/*                established in the previous section. The algorithm is */
/*                a transcription of the LQMD algorithm described by Sack */
/*                and Donovan in Numer. Mathematik 18, 465-478 (1972). */
/*                The needed data is as follows (NMOM = 2*NGQP-1): */
/*                   1) Normed modified moments: MOM (I),I=1,NMOM */
/*                   2) Recurrence coefficients of the corresponding */
/*                      monic polynomials: A (I),I=1,NMOM and B (I), */
/*                      I=2,NMOM. The recurrence relation for the monic */
/*                      polynomials is defined as: */
/*                         P (x)  =  1 */
/*                          0 */
/*                         P (x)  =  (x - A ) */
/*                          1              1 */
/*                         P   (x)  =  (x - A   ) * P (x) - B   P   (x) */
/*                          i+1              i+1     i       i+1 i-1 */
/*                The result will consist in the diagonal and offdiagonal */
/*                parts of the symmetric terminal matrix: */
/*                            DIA (I) , I = 1,NGQP */
/*                            OFF (I) , I = 1,NGQP-1 */
/*                Proceed now with Sack and Donovan's algorithm. */
/*                Handle low number (1 or 2) of quadrature points first. */
        if (*ngqp == 1)
        {
            dia[1] = mom[1] + a[1];
        }
        else if (*ngqp == 2)
        {
            sigma = mom[1] + a[1];
            dia[1] = sigma;
            theta = (a[2] - sigma) * mom[1] + mom[2] + b[2];
            off[1] = sqrt (theta);
            dia[2] = ((a[3] - sigma) * mom[2] + mom[3] + b[3] * mom[1]) /
                theta - mom[1] + a[2];
        }
        else
        {
/*             ...Handle case for number of quadrature points > 2. */
/*                Set maximum values for I and J and evaluate first */
/*                diagonal element. */
            imax = *ngqp - 1;
            jmax = *ngqp + imax;
            i__2 = jmax;
            for (j = 1; j <= i__2; ++j)
            {
                row1[j] = mom[j];
            }
            sigma = row1[1] + a[1];
            dia[1] = sigma;

/*             ...evaluate 2nd row of terminal matrix. */
            row2[1] = (a[2] - sigma) * row1[1] + row1[2] + b[2];
            theta = row2[1];
            off[1] = sqrt (theta);
            --jmax;
            i__2 = jmax;
            for (j = 2; j <= i__2; ++j)
            {
                jp1 = j + 1;
                row2[j] = (a[jp1] - sigma) * row1[j] + row1[jp1] + b[jp1] *
                    row1[j - 1];
            }
            sigma = row2[2] / theta - row1[1] + a[2];
            dia[2] = sigma;

/*             ...proceed with higher rows. */
            i__2 = imax;
            for (i__ = 2; i__ <= i__2; ++i__)
            {
                ip1 = i__ + 1;
                --jmax;
                if (i__ % 2 == 0)
                {
                    i__3 = jmax;
                    for (j = i__; j <= i__3; ++j)
                    {
                        jp1 = j + 1;
                        row1[j] =
                            (a[jp1] - sigma) * row2[j] + row2[jp1] +
                            b[jp1] * row2[j - 1] - theta * row1[j];
                    }
                    sigma = a[ip1] - row2[i__] / row2[i__ - 1] + row1[ip1] /
                        row1[i__];
                    theta = row1[i__] / row2[i__ - 1];
                }
                else
                {
                    i__3 = jmax;
                    for (j = i__; j <= i__3; ++j)
                    {
                        jp1 = j + 1;
                        row2[j] =
                            (a[jp1] - sigma) * row1[j] + row1[jp1] +
                            b[jp1] * row1[j - 1] - theta * row2[j];
                    }
                    sigma = a[ip1] - row1[i__] / row1[i__ - 1] + row2[ip1] /
                        row2[i__];
                    theta = row2[i__] / row1[i__ - 1];
                }
                dia[ip1] = sigma;
                off[i__] = sqrt (theta);
            }
        }

/*             ...The last section computes the gaussian quadrature roots */
/*                and weights from a previously established symmetric */
/*                terminal matrix by the Golub-Welsch algorithm (see */
/*                G.H. Golub and J.H. Welsch, Math. of Computation 23, */
/*                p. 221-230 and A1-A10, 1969), which is based on */
/*                a result of Wilf (see H.S. Wilf, Mathematics for the */
/*                Physical Sciences, New York: Wiley, Problem 9, p. 80). */
/*                Wilf has shown that if Z (K,I) is the K-th element of */
/*                the I-th normalized eigenvector of the terminal matrix */
/*                corresponding to the I-th eigenvalue D (I), then the */
/*                roots RTS (i.e. the zeros of the N-th orthogonal */
/*                monic polynomial) and weights WTS to be used for the */
/*                gaussian quadrature are given by: */
/*                           RTS (I) = D (I) */
/*                           WTS (I) = MOMZERO * (Z (1,I)**2) */
/*                where MOMZERO is the value of the definite integral */
/*                over the weight function W (x) alone. In our case */
/*                it is equal to the value of the zeroth Rys moment. */
/*                The present section performs hence a diagonalization */
/*                of the tridiagonal symmetric terminal matrix keeping */
/*                only the first component of the eigenvectors and sets */
/*                the roots and weights equal to the above relations. */
/*                The diagonalization code was derived from the routine */
/*                IMTQL2 in the EISPACK collection and uses the implicit */
/*                QL method. */
/*                Note, that the original diagonals DIA and offdiagonals */
/*                OFF of the terminal matrix are destroyed during the */
/*                diagonalization process !!! */
/*                Handle special case, if order of terminal matrix is 1. */
        if (*ngqp == 1)
        {
            ++nrts;
            rts[nrts] = dia[1];
            wts[nrts] = momzero;
        }
        else
        {
/*             ...initialize vector for collecting first component */
/*                of eigenvectors. To save space, array A is used, */
/*                which can be done safely, since its dimension NMOM */
/*                (# of moments) is always >= # of roots NGQP. */
            a[1] = 1.;
            i__2 = *ngqp;
            for (j = 2; j <= i__2; ++j)
            {
                a[j] = 0.;
            }

/*             ...QL iterations. */
            off[*ngqp] = 0.;
            i__2 = *ngqp;
            for (j = 1; j <= i__2; ++j)
            {
                iter = 0;
              L3000:
                i__3 = *ngqp;
                for (m = j; m <= i__3; ++m)
                {
                    if (m == *ngqp)
                    {
                        goto L3300;
                    }
                    test1 = (d__1 = dia[m], fabs (d__1)) +
                        (d__2 = dia[m + 1], fabs (d__2));
                    test2 = test1 + (d__1 = off[m], fabs (d__1));
                    if (test2 == test1)
                    {
                        goto L3300;
                    }
                }
L3300:
                p = dia[j];
                if (m == j)
                {
                    goto L300;
                }
                assert (iter != 30);                
                ++iter;
                g = (dia[j + 1] - p) / (off[j] * 2.);
                r__ = sqrt (g * g + 1.);
                g = dia[m] - p + off[j] / (g + d_sign (&r__, &g));
                s = 1.;
                c__ = 1.;
                p = 0.;
                i__3 = j;
                for (i__ = m - 1; i__ >= i__3; --i__)
                {
                    f = s * off[i__];
                    d__ = c__ * off[i__];
                    r__ = sqrt (f * f + g * g);
                    off[i__ + 1] = r__;
                    if (r__ == 0.)
                    {
                        dia[i__ + 1] -= p;
                        off[m] = 0.;
                        goto L3000;
                    }
                    s = f / r__;
                    c__ = g / r__;
                    g = dia[i__ + 1] - p;
                    r__ = (dia[i__] - g) * s + c__ * 2. * d__;
                    p = s * r__;
                    dia[i__ + 1] = g + p;
                    g = c__ * r__ - d__;
                    f = a[i__ + 1];
                    a[i__ + 1] = s * a[i__] + c__ * f;
                    a[i__] = c__ * a[i__] - s * f;
                }
                dia[j] -= p;
                off[j] = g;
                off[m] = 0.;
                goto L3000;
L300:
                ;
            }

/*             ...calculate roots and weights. Since it is known that */
/*                the roots must lay between 0 and 1, a check is */
/*                made on them to see if they are actually within this */
/*                range. */
            i__2 = *ngqp;
            for (i__ = 1; i__ <= i__2; ++i__)
            {
                root = dia[i__];
                assert (root >= 0. && root < 1.);
                rts[nrts + i__] = root;
                d__1 = a[i__];
                wts[nrts + i__] = momzero * (d__1 * d__1);
            }
            nrts += *ngqp;
        }
    }

    return 0;
}


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
