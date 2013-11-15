#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "erd.h"


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__SET_ABCD */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine handles the logistics on how to evaluate */
/*                the (12|34) integral batch in the most efficient way. */
/*                It performs the label map: */
/*                             (12|34) --> (AB|CD) */
/*                according to certain criteria that have to be met for */
/*                efficiency. The freedom we have in making the internal */
/*                association 1,2,3,4 -> A,B,C,D follows from the 8-fold */
/*                permutational symmetry of the integrals in (12|34): */
/*                     (12|34) = (21|34) = (12|43) = (21|43) = */
/*                     (34|12) = (43|12) = (34|21) = (43|21) */
/*                where the first line has only 12 -> 21 and 34 -> 43 */
/*                switches and the second line is obtained from the */
/*                first by bra <-> ket transpositions. */
/*                The type of switch to be applied is simply governed */
/*                by the demand that the final A,B,C,D shell labels */
/*                obey the relations A>=B and C>=D, since this means */
/*                the least amount of work (# of steps) for the HRR */
/*                procedure. The decision to perform a bra <-> ket */
/*                transposition comes from handling memory allocations */
/*                during both HRR on the bra and ket sides. */
/*                  Input (x = 1,2,3 and 4): */
/*                   NCGTOx        =  # of contractions for csh x */
/*                   NPGTOx        =  # of primitives per contraction */
/*                                    for csh x */
/*                   SHELLx        =  the shell type for csh x */
/*                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers */
/*                                    y = 1,2,3 and 4 */
/*                   EXPx          =  primitive exponents for csh x */
/*                   CCx           =  contraction coeffs for csh x */
/*                   SPHERIC       =  is true, if spherical integrals */
/*                                    are wanted, false if cartesian */
/*                                    ones are wanted */
/*                  Output (x = A,B,C and D): */
/*                   NCGTOx        =  # of contractions for csh x */
/*                   NPGTOx        =  # of primitives per contraction */
/*                                    for csh x */
/*                   SHELLx        =  the shell type for csh x */
/*                   SHELLy        =  the shell sums: y = P,Q,T = */
/*                                    A+B,C+D,P+Q */
/*                   MXSHELL       =  the largest (maximum) shell type */
/*                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers */
/*                                    y = A,B,C and D */
/*                   ATOMIC        =  indicates, if purely atomic */
/*                                    integrals will be evaluated */
/*                   ATOMxy        =  indicates, if atoms x and y are */
/*                                    considered to be equal for the */
/*                                    pairs xy = AB and CD */
/*                   EQUALxy       =  indicates, if csh x and csh y are */
/*                                    considered to be equal for the */
/*                                    pairs xy = AB and CD */
/*                   xxX,xxY,xxZ   =  the x,y,z-coordinate differences */
/*                                    for centers xx = AB and CD */
/*                   NxxCOOR       =  # of non-zero x,y,z-coordinate */
/*                                    differences for centers xx = AB */
/*                                    and CD */
/*                   RNxxSQ        =  square of the magnitude of the */
/*                                    distance between centers xx = AB */
/*                                    and CD */
/*                   SPNORM        =  normalization factor due to */
/*                                    presence of s- and p-type shells. */
/*                                    For each s-type shell there is a */
/*                                    factor of 1 and for each p-type */
/*                                    shell a factor of 2 */
/*                   NXYZx         =  # of cartesian monomials for csh x */
/*                   NXYZE(F)T     =  sum of # of cartesian monomials */
/*                                    for all shells in the range */
/*                                    E = A,...,A+B and in the range */
/*                                    F = C,...,C+D */
/*                   NXYZy         =  # of cartesian monomials for */
/*                                    y = P,Q shells */
/*                   NRYx          =  # of spherical functions for csh x */
/*                   INDEXx        =  index A,B,C,D -> 1,2,3,4 map */
/*                   SWAPxy        =  is .true. for xy = 12 and 34, if */
/*                                    a swap 1 <-> 2 and 3 <-> 4 has */
/*                                    been performed */
/*                   SWAPRS(TU)    =  is set .true. if the contraction */
/*                                    order of the primitives pair AB(CD) */
/*                                    will be performed in reverse order */
/*                                    BA(DC) for efficiency reasons */
/*                   TR1234        =  is .true., if a bra <-> ket */
/*                                    transposition has been applied */
/*                   LEXPx         =  pointers to locate appropriate */
/*                                    section of the exponent array */
/*                                    corresponding to csh x */
/*                   LCCx          =  pointers to locate appropriate */
/*                                    section of the contraction coeff */
/*                                    array corresponding to csh x */
/*                   LCCSEGx       =  pointers to locate appropriate */
/*                                    section of the lowest and highest */
/*                                    primitive index array defining */
/*                                    segmented contraction boundaries */
/*                                    for csh x */
/*                   NXYZHRR       =  maximum dimension of cartesian and */
/*                                    spherical components during the */
/*                                    entire HRR contraction procedure */
/*                   NCOLHRR       =  maximum # of HRR rotation matrix */
/*                                    columns needed to generate the */
/*                                    final HRR rotation matrix */
/*                   NROTHRR       =  maximum # of HRR rotation matrix */
/*                                    elements needed to generate the */
/*                                    final HRR rotation matrix */
/*                   EMPTY         =  int flag, indicating if an */
/*                                    empty batch of integrals is */
/*                                    expected. */
/* ------------------------------------------------------------------------ */
int erd__set_abcd (int ncgto1, int ncgto2,
                   int ncgto3, int ncgto4,
                   int npgto1, int npgto2,
                   int npgto3, int npgto4,
                   int shell1, int shell2,
                   int shell3, int shell4,
                   double x1, double y1, double z1,
                   double x2, double y2, double z2,
                   double x3, double y3, double z3,
                   double x4, double y4, double z4,
                   double *exp1, double *exp2,
                   double *exp3, double *exp4,
                   double *cc1, double *cc2,
                   double *cc3, double *cc4,
                   int spheric, int *ncgtoa, int *ncgtob,
                   int *ncgtoc, int *ncgtod,
                   int *npgtoa, int *npgtob,
                   int *npgtoc, int *npgtod,
                   int *shella, int *shellb,
                   int *shellc, int *shelld,
                   int *shellp, int *shellq, int *shellt,
                   int *mxshell,
                   double *xa, double *ya, double *za,
                   double *xb, double *yb, double *zb,
                   double *xc, double *yc, double *zc,
                   double *xd, double *yd, double *zd,
                   int *atomic, int *atomab, int *atomcd,
                   int *equalab, int *equalcd,
                   double *abx, double *aby, double *abz,
                   double *cdx, double *cdy, double *cdz,
                   int *nabcoor, int *ncdcoor,
                   double *rnabsq, double *rncdsq,
                   double *spnorm, int *nxyza, int *nxyzb,
                   int *nxyzc, int *nxyzd,
                   int *nxyzet, int *nxyzft,
                   int *nxyzp, int * nxyzq,
                   int *nrya, int *nryb,
                   int *nryc, int *nryd,
                   int *indexa, int *indexb,
                   int *indexc, int *indexd,
                   int * swap12, int *swap34,
                   int *swaprs, int * swaptu,
                   int *tr1234, int *lexpa, int * lexpb,
                   int *lexpc, int *lexpd,
                   int * lcca, int *lccb, int *lccc, int *lccd,
                   int *lccsega, int * lccsegb,
                   int *lccsegc, int *lccsegd,
                   int *nxyzhrr, int * ncolhrr,
                   int *nrothrr, int *empty)
{
    int add[3] = { 0, 0, 1 };

    int cc1_offset, cc2_offset, cc3_offset,
        cc4_offset;

    int i, j, m, ngh, ncc1, ncc2, ncc3, nry1, nry2, nry3, nry4,
        ngho, ncol, nrot, nrow;
    int case1, case2;
    int nxyz1, nxyz2, nxyz3, nxyz4;
    int atom12, atom23, atom34;
    int nxyze, nxyzf, nxyzg, nxyzh, nxyzi;
    int equal12, equal34;
    int shellg, shellh, nxyzgo, nxyzho, nhrr2nd, nhrr1st;

    cc1_offset = 1 + npgto1 * 1;
    cc1 -= cc1_offset;
    --exp1;
    cc2_offset = 1 + npgto2 * 1;
    cc2 -= cc2_offset;
    --exp2;
    cc3_offset = 1 + npgto3 * 1;
    cc3 -= cc3_offset;
    --exp3;
    npgto4 = npgto4 - 1 + 1;
    cc4_offset = 1 + npgto4 * 1;
    cc4 -= cc4_offset;
    --exp4;

/*             ...generate all 1,2,3,4 data. Decide as early as */
/*                possible, if a zero batch of integrals is expected. */
    *empty = 0;
    atom12 = ((x1 == x2) && (y1 == y2) && (z1 == z2));
    atom23 = ((x2 == x3) && (y2 == y3) && (z2 == z3));
    atom34 = ((x3 == x4) && (y3 == y4) && (z3 == z4));
    *atomic = atom12 && atom34 && atom23;
    *shellp = shell1 + shell2;
    *shellq = shell3 + shell4;
    *shellt = *shellp + *shellq;

    *mxshell = MAX (shell1, shell2);
    *mxshell = MAX (*mxshell, shell3);
    *mxshell = MAX (*mxshell, shell4);
    case1 = (*shellt % 2 == 1);
    case2 = (spheric && *mxshell) + (*mxshell > *shellt);
    if (*atomic && (case1 || case2))
    {
        *empty = 1;
        return 0;
    }


/*             ...determine csh equality between center pairs 1,2 */
/*                and 3,4 in increasing order of complexity: */
/*                 centers -> shells -> exponents -> ctr coefficients */
    equal12 = atom12;
    if (equal12)
    {
        equal12 = (shell1 == shell2) &&
                  (npgto1 == npgto2) &&
                  (ncgto1 == ncgto2);
        if (equal12)
        {
            for (i = 1; i <= npgto1; ++i)
            {
                equal12 = equal12 && (exp1[i] == exp2[i]);
            }
            if (equal12)
            {
                for (j = 1; j <= ncgto1; ++j)
                {
                    if (equal12)
                    {
                        for (i = 1; i <= npgto1; ++i)
                        {
                            equal12 = equal12 &&
                                (cc1[i + j * npgto1] == cc2[i + j * npgto2]);
                        }
                    }
                }
            }
        }
    }
    equal34 = atom34;
    if (equal34)
    {
        equal34 = (shell3 == shell4) &&
                  (npgto3 == npgto4) &&
                  (ncgto3 == ncgto4);
        if (equal34)
        {
            for (i = 1; i <= npgto3; ++i)
            {
                equal34 = equal34 && (exp3[i] == exp4[i]);
            }
            if (equal34)
            {
                for (j = 1; j <= ncgto3; ++j)
                {
                    if (equal34)
                    {
                        for (i = 1; i <= npgto3; ++i)
                        {
                            equal34 = equal34 &&
                                (cc3[i + j * npgto3] == cc4[i + j * npgto4]);
                        }
                    }
                }
            }
        }
    }

/*             ...set the cartesian and spherical dimensions. In case */
/*                no spherical transformations are wanted, set the */
/*                corresponding dimensions equal to the cartesian ones. */
    nxyz1 = (shell1 + 1) * (shell1 + 2) / 2;
    nxyz2 = (shell2 + 1) * (shell2 + 2) / 2;
    nxyz3 = (shell3 + 1) * (shell3 + 2) / 2;
    nxyz4 = (shell4 + 1) * (shell4 + 2) / 2;
    nry1 = shell1 + shell1 + 1;
    nry2 = shell2 + shell2 + 1;
    nry3 = shell3 + shell3 + 1;
    nry4 = shell4 + shell4 + 1;
    if (!spheric)
    {
        nry1 = nxyz1;
        nry2 = nxyz2;
        nry3 = nxyz3;
        nry4 = nxyz4;
    }

/*             ...decide on the 1 <-> 2 and/or 3 <-> 4 swapping. */
    *swap12 = (shell1 < shell2);
    *swap34 = (shell3 < shell4);
    
/*             ...calculate NXYZHRR for the two possible HRR */
/*                and (if any) cartesian -> spherical transformation */
/*                sequences: */
/*                    i) initial dimension:  NXYZE * NXYZF */
/*                   ii) perform HRR on 34:  NXYZE * NXYZ3 * NXYZ4 */
/*                  iii) cart -> sph on 34:  NXYZE * NRY3 * NRY4 */
/*                   iv) perform HRR on 12:  NXYZ1 * NXYZ2 * NRY3 * NRY4 */
/*                    v) cart -> sph on 12:  NRY1 * NRY2 * NRY3 * NRY4 */
/*                    i) initial dimension:  NXYZE * NXYZF */
/*                   ii) perform HRR on 12:  NXYZ1 * NXYZ2 * NXYZF */
/*                  iii) cart -> sph on 12:  NRY1 * NRY2 * NXYZF */
/*                   iv) perform HRR on 34:  NRY1 * NRY2 * NXYZ3 * NXYZ4 */
/*                    v) cart -> sph on 34:  NRY1 * NRY2 * NRY3 * NRY4 */
/*                The only dimension increasing steps are ii) and iv) */
/*                in both cases. Hence we first find the maximum */
/*                between ii) and iv) for both sequences and then we */
/*                take the overall minimum of these two maxima. */
/*                Since the order of sequence of the HRR on the A,B,C,D */
/*                labels is CD followed by AB, the overall minimum */
/*                will decide if to perform a bra <-> ket transposition */
/*                on the 12/34 labels. */
    *shella = MAX (shell1, shell2);
    *shellb = MIN (shell1, shell2);
    *shellc = MAX (shell3, shell4);
    *shelld = MIN (shell3, shell4);
    nxyze = (*shellp + 1) * (*shellp + 2) * (*shellp + 3) / 6 -
        *shella * (*shella + 1) * (*shella + 2) / 6;
    nxyzf = (*shellq + 1) * (*shellq + 2) * (*shellq + 3) / 6 -
        *shellc * (*shellc + 1) * (*shellc + 2) / 6;
    if (*shellb == 0 && *shelld == 0)
    {
        *nxyzhrr = nxyze * nxyzf;
        *tr1234 = 0;
    }
    else
    {
        nhrr1st = MAX (nxyze * nxyz3 * nxyz4,
                       nxyz1 * nxyz2 * nry3 * nry4);
        nhrr2nd = MAX (nxyzf * nxyz1 * nxyz2,
                       nxyz3 * nxyz4 * nry1 * nry2);
        *nxyzhrr = MIN (nhrr1st, nhrr2nd);
        *tr1234 = nhrr1st > nhrr2nd;
    }

/*             ...according to the previously gathered info, set the */
/*                new A,B,C,D shells, # of primitives + contraction */
/*                coeffs as well as pointers to the alpha exponents */
/*                and contraction coefficients. Also set the info for */
/*                evaluation of the [e0|f0] batches and for the HRR */
/*                steps later on. */
    ncc1 = npgto1 * ncgto1;
    ncc2 = npgto2 * ncgto2;
    ncc3 = npgto3 * ncgto3;
    if (!(*tr1234))
    {
        *nxyzet = nxyze;
        *nxyzft = nxyzf;
        *atomab = atom12;
        *atomcd = atom34;
        *equalab = equal12;
        *equalcd = equal34;
        if (!(*swap12))
        {
            *xa = x1;
            *ya = y1;
            *za = z1;
            *xb = x2;
            *yb = y2;
            *zb = z2;
            *shella = shell1;
            *shellb = shell2;
            *npgtoa = npgto1;
            *npgtob = npgto2;
            *ncgtoa = ncgto1;
            *ncgtob = ncgto2;
            *nxyza = nxyz1;
            *nxyzb = nxyz2;
            *nrya = nry1;
            *nryb = nry2;
            *indexa = 1;
            *indexb = 2;
            *lexpa = 1;
            *lexpb = *lexpa + npgto1;
            *lcca = 1;
            *lccb = *lcca + ncc1;
            *lccsega = 1;
            *lccsegb = *lccsega + ncgto1;
        }
        else
        {
            *xa = x2;
            *ya = y2;
            *za = z2;
            *xb = x1;
            *yb = y1;
            *zb = z1;
            *shella = shell2;
            *shellb = shell1;
            *npgtoa = npgto2;
            *npgtob = npgto1;
            *ncgtoa = ncgto2;
            *ncgtob = ncgto1;
            *nxyza = nxyz2;
            *nxyzb = nxyz1;
            *nrya = nry2;
            *nryb = nry1;
            *indexa = 2;
            *indexb = 1;
            *lexpb = 1;
            *lexpa = *lexpb + npgto1;
            *lccb = 1;
            *lcca = *lccb + ncc1;
            *lccsegb = 1;
            *lccsega = *lccsegb + ncgto1;
        }
        if (!(*swap34))
        {
            *xc = x3;
            *yc = y3;
            *zc = z3;
            *xd = x4;
            *yd = y4;
            *zd = z4;
            *shellc = shell3;
            *shelld = shell4;
            *npgtoc = npgto3;
            *npgtod = npgto4;
            *ncgtoc = ncgto3;
            *ncgtod = ncgto4;
            *nxyzc = nxyz3;
            *nxyzd = nxyz4;
            *nryc = nry3;
            *nryd = nry4;
            *indexc = 3;
            *indexd = 4;
            *lexpc = npgto1 + 1 + npgto2;
            *lexpd = *lexpc + npgto3;
            *lccc = ncc1 + 1 + ncc2;
            *lccd = *lccc + ncc3;
            *lccsegc = ncgto1 + 1 + ncgto2;
            *lccsegd = *lccsegc + ncgto3;
        }
        else
        {
            *xc = x4;
            *yc = y4;
            *zc = z4;
            *xd = x3;
            *yd = y3;
            *zd = z3;
            *shellc = shell4;
            *shelld = shell3;
            *npgtoc = npgto4;
            *npgtod = npgto3;
            *ncgtoc = ncgto4;
            *ncgtod = ncgto3;
            *nxyzc = nxyz4;
            *nxyzd = nxyz3;
            *nryc = nry4;
            *nryd = nry3;
            *indexc = 4;
            *indexd = 3;
            *lexpd = npgto1 + 1 + npgto2;
            *lexpc = *lexpd + npgto3;
            *lccd = ncc1 + 1 + ncc2;
            *lccc = *lccd + ncc3;
            *lccsegd = ncgto1 + 1 + ncgto2;
            *lccsegc = *lccsegd + ncgto3;
        }
    }
    else
    {
        *nxyzet = nxyzf;
        *nxyzft = nxyze;
        *atomab = atom34;
        *atomcd = atom12;
        *equalab = equal34;
        *equalcd = equal12;
        if (!(*swap12))
        {
            *xc = x1;
            *yc = y1;
            *zc = z1;
            *xd = x2;
            *yd = y2;
            *zd = z2;
            *shellc = shell1;
            *shelld = shell2;
            *npgtoc = npgto1;
            *npgtod = npgto2;
            *ncgtoc = ncgto1;
            *ncgtod = ncgto2;
            *nxyzc = nxyz1;
            *nxyzd = nxyz2;
            *nryc = nry1;
            *nryd = nry2;
            *indexc = 1;
            *indexd = 2;
            *lexpc = 1;
            *lexpd = *lexpc + npgto1;
            *lccc = 1;
            *lccd = *lccc + ncc1;
            *lccsegc = 1;
            *lccsegd = *lccsegc + ncgto1;
        }
        else
        {
            *xc = x2;
            *yc = y2;
            *zc = z2;
            *xd = x1;
            *yd = y1;
            *zd = z1;
            *shellc = shell2;
            *shelld = shell1;
            *npgtoc = npgto2;
            *npgtod = npgto1;
            *ncgtoc = ncgto2;
            *ncgtod = ncgto1;
            *nxyzc = nxyz2;
            *nxyzd = nxyz1;
            *nryc = nry2;
            *nryd = nry1;
            *indexc = 2;
            *indexd = 1;
            *lexpd = 1;
            *lexpc = *lexpd + npgto1;
            *lccd = 1;
            *lccc = *lccd + ncc1;
            *lccsegd = 1;
            *lccsegc = *lccsegd + ncgto1;
        }
        if (!(*swap34))
        {
            *xa = x3;
            *ya = y3;
            *za = z3;
            *xb = x4;
            *yb = y4;
            *zb = z4;
            *shella = shell3;
            *shellb = shell4;
            *npgtoa = npgto3;
            *npgtob = npgto4;
            *ncgtoa = ncgto3;
            *ncgtob = ncgto4;
            *nxyza = nxyz3;
            *nxyzb = nxyz4;
            *nrya = nry3;
            *nryb = nry4;
            *indexa = 3;
            *indexb = 4;
            *lexpa = npgto1 + 1 + npgto2;
            *lexpb = *lexpa + npgto3;
            *lcca = ncc1 + 1 + ncc2;
            *lccb = *lcca + ncc3;
            *lccsega = ncgto1 + 1 + ncgto2;
            *lccsegb = *lccsega + ncgto3;
        }
        else
        {
            *xa = x4;
            *ya = y4;
            *za = z4;
            *xb = x3;
            *yb = y3;
            *zb = z3;
            *shella = shell4;
            *shellb = shell3;
            *npgtoa = npgto4;
            *npgtob = npgto3;
            *ncgtoa = ncgto4;
            *ncgtob = ncgto3;
            *nxyza = nxyz4;
            *nxyzb = nxyz3;
            *nrya = nry4;
            *nryb = nry3;
            *indexa = 4;
            *indexb = 3;
            *lexpb = npgto1 + 1 + npgto2;
            *lexpa = *lexpb + npgto3;
            *lccb = ncc1 + 1 + ncc2;
            *lcca = *lccb + ncc3;
            *lccsegb = ncgto1 + 1 + ncgto2;
            *lccsega = *lccsegb + ncgto3;
        }
    }

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
    *swaprs = (*npgtoa > *npgtob);
    *swaptu = (*npgtoc > *npgtod);
    *shellp = *shella + *shellb;
    *shellq = *shellc + *shelld;
    *nxyzp = (*shellp + 1) * (*shellp + 2) / 2;
    *nxyzq = (*shellq + 1) * (*shellq + 2) / 2;
    *spnorm = 1.;
    if (*shella == 1)
    {
        *spnorm += *spnorm;
    }
    if (*shellb == 1)
    {
        *spnorm += *spnorm;
    }
    if (*shellc == 1)
    {
        *spnorm += *spnorm;
    }
    if (*shelld == 1)
    {
        *spnorm += *spnorm;
    }


/*             ...calculate the coordinate differences between centers */
/*                A and B and between centers C and D and calculate the */
/*                square of the magnitude of the distances. Also */
/*                determine the number of non-zero coordinate differences */
/*                for each pair A,B and C,D. */
    if (!(*atomab))
    {
        *abx = *xa - *xb;
        *aby = *ya - *yb;
        *abz = *za - *zb;
        *rnabsq = *abx * *abx + *aby * *aby + *abz * *abz;
        *nabcoor = 3;
        if (fabs (*abx) == 0.0)
        {
            --(*nabcoor);
        }
        if (fabs (*aby) == 0.0)
        {
            --(*nabcoor);
        }
        if (fabs (*abz) == 0.0)
        {
            --(*nabcoor);
        }
    }
    else
    {
        *abx = 0.0;
        *aby = 0.0;
        *abz = 0.0;
        *rnabsq = 0.0;
        *nabcoor = 0;
    }
    if (!(*atomcd))
    {
        *cdx = *xc - *xd;
        *cdy = *yc - *yd;
        *cdz = *zc - *zd;
        *rncdsq = *cdx * *cdx + *cdy * *cdy + *cdz * *cdz;
        *ncdcoor = 3;
        if (fabs (*cdx) == 0.)
        {
            --(*ncdcoor);
        }
        if (fabs (*cdy) == 0.)
        {
            --(*ncdcoor);
        }
        if (fabs (*cdz) == 0.)
        {
            --(*ncdcoor);
        }
    }
    else
    {
        *cdx = 0.;
        *cdy = 0.;
        *cdz = 0.;
        *rncdsq = 0.;
        *ncdcoor = 0;
    }

/*             ...if HRR contractions are to be performed, calculate */
/*                NCOLHRR (maximum # of HRR rotation matrix columns */
/*                needed to generate the final HRR rotation matrices) and */
/*                NROTHRR (maximum # of HRR rotation matrix elements). */
/*                First find maximum values for the HRR on the AB-part. */
    *ncolhrr = 0;
    *nrothrr = 0;
    if (*shellb != 0)
    {
        ngh = *nxyzet;
        nxyzg = *nxyzet;
        nxyzh = 1;
        nxyzgo = *nxyzet;
        nxyzho = 1;
        nxyzi = *nxyzp;
        shellg = *shellp;
        nrow = 1;
        ncol = ngh;
        nrot = ngh;
        for (shellh = 1; shellh <= *shellb; ++shellh)
        {
            nxyzgo -= nxyzi;
            nxyzho = nxyzho + shellh + 1;
            ngho = nxyzgo * nxyzho;
            if (*nabcoor == 3)
            {
                m = shellh / 3 + 1;
                nrow += m * (m + add[shellh % 3]);
            }
            else if (*nabcoor == 2)
            {
                nrow = nrow + shellh / 2 + 1;
            }
            else if (*nabcoor == 1)
            {
                ++nrow;
            }
            ncol = MAX (ngho, ncol);
            nrot = MAX (nrow * ngho, nrot);
            ngh = ngho;
            nxyzg = nxyzgo;
            nxyzh = nxyzho;
            nxyzi = nxyzi - shellg - 1;
            --shellg;
        }
        *ncolhrr = ncol;
        *nrothrr = nrot;
    }


/*             ...next find maximum values for the HRR on the CD-part */
/*                and set overall maximum values. */
    if (*shelld != 0)
    {
        ngh = *nxyzft;
        nxyzg = *nxyzft;
        nxyzh = 1;
        nxyzgo = *nxyzft;
        nxyzho = 1;
        nxyzi = *nxyzq;
        shellg = *shellq;
        nrow = 1;
        ncol = ngh;
        nrot = ngh;
        for (shellh = 1; shellh <= *shelld; ++shellh)
        {
            nxyzgo -= nxyzi;
            nxyzho = nxyzho + shellh + 1;
            ngho = nxyzgo * nxyzho;
            if (*ncdcoor == 3)
            {
                m = shellh / 3 + 1;
                nrow += m * (m + add[shellh % 3]);
            }
            else if (*ncdcoor == 2)
            {
                nrow = nrow + shellh / 2 + 1;
            }
            else if (*ncdcoor == 1)
            {
                ++nrow;
            }
            ncol = MAX (ngho, ncol);
            nrot = MAX (nrow * ngho, nrot);
            ngh = ngho;
            nxyzg = nxyzgo;
            nxyzh = nxyzho;
            nxyzi = nxyzi - shellg - 1;
            --shellg;
        }
        *ncolhrr = MAX (ncol, *ncolhrr);
        *nrothrr = MAX (nrot, *nrothrr);
    }

    return 0;
}


int erd__set_abcd_ (int * ncgto1, int * ncgto2,
                 int * ncgto3, int * ncgto4, int * npgto1,
                 int * npgto2, int * npgto3, int * npgto4,
                 int * shell1, int * shell2, int * shell3,
                 int * shell4, double * x1, double * y1,
                 double * z1, double * x2, double * y2,
                 double * z2, double * x3, double * y3,
                 double * z3, double * x4, double * y4,
                 double * z4, double * exp1, double * exp2,
                 double * exp3, double * exp4, double * cc1,
                 double * cc2, double * cc3, double * cc4,
                 int * spheric, int * ncgtoa, int * ncgtob,
                 int * ncgtoc, int * ncgtod, int * npgtoa,
                 int * npgtob, int * npgtoc, int * npgtod,
                 int * shella, int * shellb, int * shellc,
                 int * shelld, int * shellp, int * shellq,
                 int * shellt, int * mxshell, double * xa,
                 double * ya, double * za, double * xb,
                 double * yb, double * zb, double * xc,
                 double * yc, double * zc, double * xd,
                 double * yd, double * zd, int * atomic,
                 int * atomab, int * atomcd, int * equalab,
                 int * equalcd, double * abx, double * aby,
                 double * abz, double * cdx, double * cdy,
                 double * cdz, int * nabcoor, int * ncdcoor,
                 double * rnabsq, double * rncdsq,
                 double * spnorm, int * nxyza, int * nxyzb,
                 int * nxyzc, int * nxyzd, int * nxyzet,
                 int * nxyzft, int * nxyzp, int * nxyzq,
                 int * nrya, int * nryb, int * nryc,
                 int * nryd, int * indexa, int * indexb,
                 int * indexc, int * indexd, int * swap12,
                 int * swap34, int * swaprs, int * swaptu,
                 int * tr1234, int * lexpa, int * lexpb,
                 int * lexpc, int * lexpd, int * lcca,
                 int * lccb, int * lccc, int * lccd,
                 int * lccsega, int * lccsegb, int * lccsegc,
                 int * lccsegd, int * nxyzhrr, int * ncolhrr,
                 int * nrothrr, int * empty)
{
    erd__set_abcd (*ncgto1, *ncgto2, *ncgto3, *ncgto4,
                   *npgto1, *npgto2, *npgto3, *npgto4,
                   *shell1, *shell2, *shell3, *shell4,
                   *x1, *y1, *z1,
                   *x2, *y2, *z2,
                   *x3, *y3, *z3,
                   *x4, *y4, *z4,
                   exp1, exp2,  exp3, exp4,
                   cc1, cc2, cc3, cc4,
                   *spheric, ncgtoa, ncgtob,
                   ncgtoc, ncgtod,
                   npgtoa, npgtob,
                   npgtoc, npgtod,
                   shella, shellb,
                   shellc, shelld,
                   shellp, shellq, shellt,
                   mxshell,
                   xa, ya, za,
                   xb, yb, zb,
                   xc, yc, zc,
                   xd, yd, zd,
                   atomic, atomab, atomcd,
                   equalab, equalcd,
                   abx, aby, abz,
                   cdx, cdy, cdz,
                   nabcoor, ncdcoor,
                   rnabsq, rncdsq,
                   spnorm, nxyza, nxyzb,
                   nxyzc, nxyzd,
                   nxyzet, nxyzft,
                   nxyzp, nxyzq,
                   nrya, nryb,
                   nryc, nryd,
                   indexa, indexb,
                   indexc, indexd,
                   swap12, swap34,
                   swaprs, swaptu,
                   tr1234, lexpa, lexpb,
                   lexpc, lexpd,
                   lcca, lccb, lccc, lccd,
                   lccsega, lccsegb,
                   lccsegc, lccsegd,
                   nxyzhrr, ncolhrr,
                   nrothrr, empty);
    
    return 0;
}