#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


#include "erd.h"


int erd__e0f0_pcgto_block_ (int * nbatch, int * nint2d,
                         int * atomic, int * atomab, int * atomcd,
                         int * mij, int * mkl, int * mijkl,
                         int * nij, int * nijbeg, int * nijend,
                         int * nkl, int * nklbeg, int * nklend,
                         int * ngqp, int * nmom, int * ngqscr,
                         int * mgqijkl, int * npgtoa,
                         int * npgtob, int * npgtoc, int * npgtod,
                         int * nxyzet, int * nxyzft, int * nxyzp,
                         int * nxyzq, int * shella, int * shellp,
                         int * shellc, int * shellq, double * xa,
                         double * ya, double * za, double * xb,
                         double * yb, double * zb, double * xc,
                         double * yc, double * zc, double * xd,
                         double * yd, double * zd, double * abx,
                         double * aby, double * abz, double * cdx,
                         double * cdy, double * cdz,
                         double * alphaa, double * alphab,
                         double * alphac, double * alphad,
                         double * ftable, int * mgrid,
                         int * ngrid, double * tmax,
                         double * tstep, double * tvstep,
                         int * prima, int * primb, int * primc,
                         int * primd, double * norma,
                         double * normb, double * normc,
                         double * normd, double * rhoab,
                         double * rhocd, double * p, double * px,
                         double * py, double * pz, double * pax,
                         double * pay, double * paz,
                         double * pinvhf, double * scalep,
                         double * q, double * qx, double * qy,
                         double * qz, double * qcx, double * qcy,
                         double * qcz, double * qinvhf,
                         double * scaleq, double * rts,
                         double * wts, double * gqscr,
                         double * tval, double * pqpinv,
                         double * scalepq, double * b00,
                         double * b01, double * b10,
                         double * c00x, double * c00y,
                         double * c00z, double * d00x,
                         double * d00y, double * d00z,
                         double * int2dx, double * int2dy,
                         double * int2dz, double * batch)
{
    /* System generated locals */
    int ftable_dim1, ftable_offset, i_1, i_2;

    /* Builtin functions */
    double sqrt (double);

    static int i_, j, k, l, m, n;
    static int ij, kl, g000, g010, g020, g030, g040, g050, g060;
    static double pqx, pqy, pqz;

    static double expa, expb, expc, expd, pval, qval;

    static double pinv, qinv, pxval, pyval, pzval;
    static int case2d;

    static int caseat;

    static double pscale, invers, pqmult, pqplus, rnpqsq;


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD_E0F0_PCGTO_BLOCK */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD_RYS_ROOTS_WEIGHTS */
/*                ERD_2D_COEFFICIENTS */
/*                ERD_2D_PQ_INTEGRALS */
/*                ERD_INT2D_TO_E000 */
/*                ERD_INT2D_TO_E0F0 */
/*                ERD_2D_ATOM_COEFFICIENTS */
/*                ERD_2D_ATOM_PQ_INTEGRALS */
/*                ERD_ATOM_INT2D_TO_E000 */
/*                ERD_ATOM_INT2D_TO_E0F0 */
/*  DESCRIPTION : This operation calculates a batch of unnormed electron */
/*                repulsion integrals between primitive cartesian */
/*                gaussians for the shell quadruplet range: */

/*                    [E0|F0]     , E = A to P, F = C to Q */
/*                           ijkl */

/*                and the block of ij and kl exponent pairs. The total */
/*                number of eris generated here is thus given by the */
/*                total number of cartesian monomials NXYZET*NXYZFT */
/*                times the total number of exponent pairs MIJKL in the */
/*                present block. */

/*                On exit, the batch elements will be stored as: */

/*                             batch (kl,ij,nxyzt) */


/*                  Input: */

/*                    NBATCH       =  size of the primitive cartesian */
/*                                    integral batch */
/*                    NINT2D       =  space needed for each of the 2D */
/*                                    X,Y,Z integral arrays */
/*                    ATOMIC       =  indicates, if purely atomic */
/*                                    integrals will be evaluated */
/*                    ATOMAB(CD)   =  indicates, if centers A and B */
/*                                    (C and D) coincide */
/*                    MIJ(KL)      =  current # of ij (kl) primitive */
/*                                    index pairs corresponding to */
/*                                    the contracted shell pairs A,B */
/*                                    (C,D) */
/*                    MIJKL        =  current # of ijkl primitive */
/*                                    index quadruplets (= MIJ*MKL) */
/*                    NIJ          =  total # of ij primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair A,B */
/*                    NIJBEG(END)  =  first(last) ij primitive index */
/*                                    defining the ij block */
/*                    NKL          =  total # of kl primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair C,D */
/*                    NKLBEG(END)  =  first(last) kl primitive index */
/*                                    defining the kl block */
/*                    NGQP         =  # of gaussian quadrature points */
/*                                    (roots) */
/*                    NMOM         =  # of necessary moment integrals */
/*                                    to calculate the quadrature roots */
/*                    NGQSCR       =  size of gaussian quadrature */
/*                                    scratch space needed to calculate */
/*                                    all the quadrature roots */
/*                    MGQIJKL      =  # of roots times # of ijkl */
/*                                    quadruplets (= NGQP*MIJKL) */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for contraction shells x = A,B,C,D */
/*                    NXYZE(F)T    =  sum of # of cartesian monomials */
/*                                    for all shells in the range */
/*                                    E = A,...,P=A+B and in the range */
/*                                    F = C,...,Q=C+D */
/*                    NXYZP(Q)     =  # of cartesian monomials for */
/*                                    the P=A+B and Q=C+D shells */
/*                    SHELLx       =  the shell type for contraction */
/*                                    shells x = A,P=A+B,C,Q=C+D */
/*                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers */
/*                                    x = A,B,C,D */
/*                    ABm(CDm)     =  the m=x,y,z-coordinate differences */
/*                                    between centers A and B (C and D) */
/*                    ALPHAx       =  the primitive exponents for */
/*                                    contraction shells x = A,B,C,D */
/*                    FTABLE       =  Fm (T) table for interpolation */
/*                                    in low T region */
/*                    MGRID        =  maximum m in Fm (T) table */
/*                    NGRID        =  # of T's for which Fm (T) table */
/*                                    was set up */
/*                    TMAX         =  maximum T in Fm (T) table */
/*                    TSTEP        =  difference between two consecutive */
/*                                    T's in Fm (T) table */
/*                    TVSTEP       =  Inverse of TSTEP */
/*                    PRIMx        =  i,j,k,l labels of primitives for */
/*                                    the respective contraction shells */
/*                                    x = A,B,C,D */
/*                    NORMx        =  the normalization factors due to */
/*                                    the primitive exponents for the */
/*                                    contraction shells x = A,B,C,D */
/*                    RHOAB(CD)    =  the complete set of NIJ (NKL) */
/*                                    exponential prefactors between */
/*                                    contraction shells A and B */
/*                                    (C and D) */
/*                    P            =  will hold current MIJ exponent */
/*                                    sums for contraction shells A */
/*                                    and B */
/*                    Px           =  will hold current MIJ coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers P=A+B */
/*                    PAx          =  will hold current MIJ coordinate */
/*                                    x=X,Y,Z differences P-A between */
/*                                    centers P and A */
/*                    PINVHF       =  will hold current MIJ values of */
/*                                    1/(2*P), where P are the exponent */
/*                                    sums for contraction shells A */
/*                                    and B */
/*                    SCALEP       =  will hold current MIJ values of */
/*                                    scaling factors related to point P */
/*                    Q            =  will hold current MKL exponent */
/*                                    sums for contraction shells C */
/*                                    and D */
/*                    Qx           =  will hold current MKL coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers Q=C+D */
/*                    QCx          =  will hold current MKL coordinate */
/*                                    x=X,Y,Z differences Q-C between */
/*                                    centers Q and C */
/*                    QINVHF       =  will hold current MKL values of */
/*                                    1/(2*Q), where Q are the exponent */
/*                                    sums for contraction shells C */
/*                                    and D */
/*                    SCALEQ       =  will hold current MKL values of */
/*                                    scaling factors related to point Q */
/*                    RTS          =  will hold all current MGQIJKL */
/*                                    quadrature roots */
/*                    WTS          =  will hold all current MGQIJKL */
/*                                    quadrature weights */
/*                    GQSCR        =  will be used as scratch space */
/*                                    for determining the quadrature */
/*                                    roots and weights */
/*                    TVAL         =  will hold current MIJKL values */
/*                                    of T-exponents defining the Rys */
/*                                    weight functions */
/*                    PQPINV       =  will hold current MIJKL values */
/*                                    of 1/(P+Q), i.e. the inverses */
/*                                    of all total exponent sums */
/*                    SCALEPQ      =  will hold current distinct MIJKL */
/*                                    (expanded to MGQIJKL) values of */
/*                                    the overal scaling factors for */
/*                                    the integrals */
/*                    Bxx          =  will hold the current MGQIJKL */
/*                                    coordinate independent VRR */
/*                                    B-coefficients (xx=00,01,10) */
/*                    C00x         =  will hold the current MGQIJKL */
/*                                    VRR C-coefficients (individual */
/*                                    cartesian components x=X,Y,Z) for */
/*                                    shell expansion on center P */
/*                    D00x         =  will hold the current MGQIJKL */
/*                                    VRR D-coefficients (individual */
/*                                    cartesian components x=X,Y,Z) for */
/*                                    shell expansion on center Q */
/*                    INT2Dx       =  will hold all current 2D integrals */
/*                                    for each cartesian component */
/*                                    (x = X,Y,Z) */

/*                  Output: */

/*                    BATCH        =  current batch of primitive */
/*                                    cartesian [E0|F0] integrals */



/*  AUTHOR      : Norbert Flocke */
/* ------------------------------------------------------------------------ */


/*             ...include files and declare variables. */




/* ------------------------------------------------------------------------ */


/*             ...predetermine 2D integral case. This is done in */
/*                order to distinguish the P- and Q-shell combinations */
/*                for efficient evaluation of the 2D integrals and */
/*                their VRR coefficients. The cases distinguished */
/*                are summarized in the following table, indicating */
/*                the value of CASE2D: */

/*                                Q-shell */
/*                              s    p   >p */
/*                            --------------- */
/*                           | */
/*                         s |  1    4    7 */
/*                           | */
/*                P-shell  p |  2    5    8 */
/*                           | */
/*                        >p |  3    6    9 */
/*                           | */



    /* Parameter adjustments */
    --batch;
    --int2dz;
    --int2dy;
    --int2dx;
    --scalep;
    --pinvhf;
    --paz;
    --pay;
    --pax;
    --pz;
    --py;
    --px;
    --p;
    --primb;
    --prima;
    --scaleq;
    --qinvhf;
    --qcz;
    --qcy;
    --qcx;
    --qz;
    --qy;
    --qx;
    --q;
    --primd;
    --primc;
    --pqpinv;
    --tval;
    --rhoab;
    --rhocd;
    --gqscr;
    --d00z;
    --d00y;
    --d00x;
    --c00z;
    --c00y;
    --c00x;
    --b10;
    --b01;
    --b00;
    --scalepq;
    --wts;
    --rts;
    --norma;
    --alphaa;
    --normb;
    --alphab;
    --normc;
    --alphac;
    --normd;
    --alphad;
    ftable_dim1 = *mgrid - 0 + 1;
    ftable_offset = 0 + ftable_dim1 * 0;
    ftable -= ftable_offset;

    /* Function Body */
    case2d = MIN(2, *shellq) * 3 + MIN(2, *shellp) + 1;


/*             ...predetermine in 'K2' loops the quantities associated */
/*                with the A,B-part and C,D-part. Set the atom equality */
/*                case CASEAT here to exploit simplifications due to */
/*                center equalities later on: */

/*                     CASEAT = 1  -->    atomic (AA|AA) integrals */
/*                            = 2  -->  2-center (AA|CC) integrals */
/*                            = 3  -->  3-center (AB|CC) integrals */
/*                            = 4  -->  3-center (AA|CD) integrals */
/*                            = 5  -->  4-center (AB|CD) integrals */


    caseat = 5;
    if (*atomab)
    {
        m = 0;
        i_1 = *nijend;
        for (ij = *nijbeg; ij <= i_1; ++ij)
        {
            ++m;
            i_ = prima[m];
            j = primb[m];
            pval = alphaa[i_] + alphab[j];
            p[m] = pval;
            px[m] = *xa;
            py[m] = *ya;
            pz[m] = *za;
            pinvhf[m] = .5 / pval;
            scalep[m] = norma[i_] * normb[j];
/* L100: */
        }
        --caseat;
    }
    else
    {
        m = 0;
        i_1 = *nijend;
        for (ij = *nijbeg; ij <= i_1; ++ij)
        {
            ++m;
            i_ = prima[m];
            j = primb[m];
            expa = alphaa[i_];
            expb = alphab[j];
            pval = expa + expb;
            pinv = 1. / pval;
            p[m] = pval;
            pval = -expb * pinv;
            pax[m] = pval * *abx;
            pay[m] = pval * *aby;
            paz[m] = pval * *abz;
            px[m] = pax[m] + *xa;
            py[m] = pay[m] + *ya;
            pz[m] = paz[m] + *za;
            pinvhf[m] = pinv * .5;
            scalep[m] = norma[i_] * normb[j] * rhoab[ij];
/* L110: */
        }
    }
    if (*atomcd)
    {
        m = 0;
        i_1 = *nklend;
        for (kl = *nklbeg; kl <= i_1; ++kl)
        {
            ++m;
            k = primc[m];
            l = primd[m];
            qval = alphac[k] + alphad[l];
            q[m] = qval;
            qx[m] = *xc;
            qy[m] = *yc;
            qz[m] = *zc;
            qinvhf[m] = .5 / qval;
            scaleq[m] = normc[k] * normd[l];
/* L200: */
        }
        caseat += -2;
    }
    else
    {
        m = 0;
        i_1 = *nklend;
        for (kl = *nklbeg; kl <= i_1; ++kl)
        {
            ++m;
            k = primc[m];
            l = primd[m];
            expc = alphac[k];
            expd = alphad[l];
            qval = expc + expd;
            qinv = 1. / qval;
            q[m] = qval;
            qval = -expd * qinv;
            qcx[m] = qval * *cdx;
            qcy[m] = qval * *cdy;
            qcz[m] = qval * *cdz;
            qx[m] = qcx[m] + *xc;
            qy[m] = qcy[m] + *yc;
            qz[m] = qcz[m] + *zc;
            qinvhf[m] = qinv * .5;
            scaleq[m] = normc[k] * normd[l] * rhocd[kl];
/* L220: */
        }
    }
    if (*atomic)
    {
        --caseat;
    }


/*             ...the 'K4' loop over all ij- and kl-exponent pairs */
/*                in present ij and kl block to calculate all T's */
/*                and scaling factors for the cases: */

/*                     CASEAT = 1  -->    atomic (AA|AA) integrals */
/*                            = 2  -->  2-center (AA|CC) integrals */
/*                            = 3  -->  3-center (AB|CC) integrals */
/*                            = 4  -->  3-center (AA|CD) integrals */
/*                            = 5  -->  4-center (AB|CD) integrals */

/*                4-center (AB|CD) integrals are checked first */
/*                (most common occurence in large systems). */


    if (caseat == 5)
    {
        m = 1;
        i_1 = *mij;
        for (ij = 1; ij <= i_1; ++ij)
        {
            pval = p[ij];
            pxval = px[ij];
            pyval = py[ij];
            pzval = pz[ij];
            pscale = scalep[ij];
            i_2 = *mkl;
            for (kl = 1; kl <= i_2; ++kl)
            {
                qval = q[kl];
                pqmult = pval * qval;
                pqplus = pval + qval;
                invers = 1. / pqplus;
                pqx = pxval - qx[kl];
                pqy = pyval - qy[kl];
                pqz = pzval - qz[kl];
                tval[m] = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult *
                    invers;
                pqpinv[m] = invers;
                scalepq[m] = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
                ++m;
/* L5500: */
            }
/* L5000: */
        }
    }
    else if (caseat == 4)
    {
        m = 1;
        i_1 = *mij;
        for (ij = 1; ij <= i_1; ++ij)
        {
            pval = p[ij];
            pscale = scalep[ij];
            i_2 = *mkl;
            for (kl = 1; kl <= i_2; ++kl)
            {
                qval = q[kl];
                pqmult = pval * qval;
                pqplus = pval + qval;
                invers = 1. / pqplus;
                pqx = *xa - qx[kl];
                pqy = *ya - qy[kl];
                pqz = *za - qz[kl];
                tval[m] = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult *
                    invers;
                pqpinv[m] = invers;
                scalepq[m] = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
                ++m;
/* L4400: */
            }
/* L4000: */
        }
    }
    else if (caseat == 3)
    {
        m = 1;
        i_1 = *mij;
        for (ij = 1; ij <= i_1; ++ij)
        {
            pval = p[ij];
            pxval = px[ij];
            pyval = py[ij];
            pzval = pz[ij];
            pqx = pxval - *xc;
            pqy = pyval - *yc;
            pqz = pzval - *zc;
            rnpqsq = pqx * pqx + pqy * pqy + pqz * pqz;
            pscale = scalep[ij];
            i_2 = *mkl;
            for (kl = 1; kl <= i_2; ++kl)
            {
                qval = q[kl];
                pqmult = pval * qval;
                pqplus = pval + qval;
                invers = 1. / pqplus;
                tval[m] = rnpqsq * pqmult * invers;
                pqpinv[m] = invers;
                scalepq[m] = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
                ++m;
/* L3300: */
            }
/* L3000: */
        }
    }
    else if (caseat == 2)
    {
        pqx = *xa - *xc;
        pqy = *ya - *yc;
        pqz = *za - *zc;
        rnpqsq = pqx * pqx + pqy * pqy + pqz * pqz;
        m = 1;
        i_1 = *mij;
        for (ij = 1; ij <= i_1; ++ij)
        {
            pval = p[ij];
            pscale = scalep[ij];
            i_2 = *mkl;
            for (kl = 1; kl <= i_2; ++kl)
            {
                qval = q[kl];
                pqmult = pval * qval;
                pqplus = pval + qval;
                invers = 1. / pqplus;
                tval[m] = rnpqsq * pqmult * invers;
                pqpinv[m] = invers;
                scalepq[m] = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
                ++m;
/* L2200: */
            }
/* L2000: */
        }
    }
    else
    {
        m = 1;
        i_1 = *mij;
        for (ij = 1; ij <= i_1; ++ij)
        {
            pval = p[ij];
            pscale = scalep[ij];
            i_2 = *mkl;
            for (kl = 1; kl <= i_2; ++kl)
            {
                qval = q[kl];
                pqplus = pval + qval;
                tval[m] = 0.;
                pqpinv[m] = 1. / pqplus;
                scalepq[m] = pscale * scaleq[kl] / (pval * qval * sqrt (pval +
                                                                        qval));
                ++m;
/* L1100: */
            }
/* L1000: */
        }
    }


/*             ...if necessary, expand the scaling array size from */
/*                MIJKL to MGQIJKL starting from the last elements. */


    if (*ngqp > 1)
    {
        n = *mgqijkl + 1;
        for (m = *mijkl; m >= 1; --m)
        {
            i_1 = *ngqp;
            for (i_ = 1; i_ <= i_1; ++i_)
            {
                scalepq[n - i_] = scalepq[m];
/* L20: */
            }
            n -= *ngqp;
/* L10: */
        }
    }


/*             ...determine memory allocation offsets for the scratch */
/*                arrays used to calculate the quadrature roots + */
/*                weights: */

/*                   G000 = offset for A coefficients (Jacobi/Laguerre) */
/*                   G010 = offset for B coefficients (Jacobi/Laguerre) */
/*                   G020 = offset for moments (Jacobi/Laguerre) */
/*                   G030 = offset for diagonals of symmetric termat */
/*                   G040 = offset for offdiagonals of symmetric termat */
/*                   G050 = offset for first row intermediates during */
/*                          evaluation of symmetric termat */
/*                   G060 = offset for second row intermediates during */
/*                          evaluation of symmetric termat */


    g000 = 1;
    g010 = g000 + *nmom;
    g020 = g010 + *nmom - 1;
    g030 = g020 + *nmom;
    g040 = g030 + *ngqp;
    g050 = g040 + *ngqp;
    g060 = g050 + *nmom;


/*             ...calculate all roots and weights. Array B00 is passed */
/*                as a scratch array. */


    erd__rys_roots_weights_ (mijkl, mgqijkl, ngqp, nmom, &tval[1], &b00[1],
                              &ftable[ftable_offset], mgrid, ngrid, tmax,
                              tstep, tvstep, &gqscr[g000], &gqscr[g010],
                              &gqscr[g020], &gqscr[g030], &gqscr[g040],
                              &gqscr[g050], &gqscr[g060], &rts[1], &wts[1]);


/*             ...perform the following steps: */

/*                1) generate all VRR coefficients. */

/*                2) construct all 2D PQ x,y,z integrals using all the */
/*                   weights and all the generated VRR coefficients for */
/*                   all exponent quadruples. */

/*                3) assemble the complete [E0|F0] batch for all ij and */
/*                   kl pairs using the 2D integrals. Arrays B00 and B01 */
/*                   are passed as scratch arrays. */

/*                The last step 3) is the most compute intensive and */
/*                separate routines are provided depending on presence */
/*                of s-shells. The gainings are in the innermost loops */
/*                of these routines, which are considerably simplified */
/*                for the special s-shell cases. Note, that the case */
/*                in which both P- and Q-shells are s-shells cannot */
/*                arise, as this case is dealt with in separate routines. */


    if (!(*atomic))
    {
        erd__2d_coefficients_ (mij, mkl, mijkl, ngqp, mgqijkl, atomab,
                                atomcd, &p[1], &q[1], &px[1], &py[1], &pz[1],
                                &qx[1], &qy[1], &qz[1], &pax[1], &pay[1],
                                &paz[1], &qcx[1], &qcy[1], &qcz[1],
                                &pinvhf[1], &qinvhf[1], &pqpinv[1], &rts[1],
                                &case2d, &b00[1], &b01[1], &b10[1], &c00x[1],
                                &c00y[1], &c00z[1], &d00x[1], &d00y[1],
                                &d00z[1]);
        erd__2d_pq_integrals_ (shellp, shellq, mgqijkl, &wts[1], &b00[1],
                                &b01[1], &b10[1], &c00x[1], &c00y[1],
                                &c00z[1], &d00x[1], &d00y[1], &d00z[1],
                                &case2d, &int2dx[1], &int2dy[1], &int2dz[1]);
        if (*shellq == 0)
        {
            erd__int2d_to_e000_ (shella, shellp, ngqp, mijkl, mgqijkl,
                                  nxyzet, nxyzp, &int2dx[1], &int2dy[1],
                                  &int2dz[1], &b00[1], &b01[1], &scalepq[1],
                                  &batch[1]);
        }
        else if (*shellp == 0)
        {
            erd__int2d_to_e000_ (shellc, shellq, ngqp, mijkl, mgqijkl,
                                  nxyzft, nxyzq, &int2dx[1], &int2dy[1],
                                  &int2dz[1], &b00[1], &b01[1], &scalepq[1],
                                  &batch[1]);
        }
        else
        {
            erd__int2d_to_e0f0_ (shella, shellp, shellc, shellq, ngqp, mijkl,
                                  mgqijkl, nxyzet, nxyzft, nxyzp, nxyzq,
                                  &int2dx[1], &int2dy[1], &int2dz[1], &b00[1],
                                  &b01[1], &scalepq[1], &batch[1]);
        }
    }
    else
    {
        erd__2d_atom_coefficients_ (mij, mkl, mijkl, ngqp, mgqijkl, &p[1],
                                     &pinvhf[1], &qinvhf[1], &pqpinv[1],
                                     &rts[1], &case2d, &b00[1], &b01[1],
                                     &b10[1]);
        erd__2d_atom_pq_integrals_ (shellp, shellq, mgqijkl, &wts[1],
                                     &b00[1], &b01[1], &b10[1], &case2d,
                                     &int2dx[1], &int2dy[1], &int2dz[1]);
        if (*shellq == 0)
        {
            erd__atom_int2d_to_e000_ (shella, shellp, ngqp, mijkl, mgqijkl,
                                       nxyzet, nxyzp, &int2dx[1], &int2dy[1],
                                       &int2dz[1], &b00[1], &b01[1],
                                       &scalepq[1], &batch[1]);
        }
        else if (*shellp == 0)
        {
            erd__atom_int2d_to_e000_ (shellc, shellq, ngqp, mijkl, mgqijkl,
                                       nxyzft, nxyzq, &int2dx[1], &int2dy[1],
                                       &int2dz[1], &b00[1], &b01[1],
                                       &scalepq[1], &batch[1]);
        }
        else
        {
            erd__atom_int2d_to_e0f0_ (shella, shellp, shellc, shellq, ngqp,
                                       mijkl, mgqijkl, nxyzet, nxyzft, nxyzp,
                                       nxyzq, &int2dx[1], &int2dy[1],
                                       &int2dz[1], &b00[1], &b01[1],
                                       &scalepq[1], &batch[1]);
        }
    }


/*             ...ready! */


    return 0;
}                               /* erd_e0f0_pcgto_block_ */
