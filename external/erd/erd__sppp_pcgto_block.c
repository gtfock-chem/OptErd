#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "boys.h"
//#define ERD_TABLE_FREE_BOYS_FUNCTIONS

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__SPPP_PCGTO_BLOCK */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation is designed to provide ultrafast block */
/*                evaluation of a batch of normalized electron repulsion */
/*                integrals between s-shell and p-shell primitive */
/*                spherical gaussian type orbitals. */
/*                A batch is defined here as containing all possible */
/*                integrals, that is its dimension is determined by */
/*                the total number of primitive functions (here = 27) */
/*                times the total number of ij and kl exponent pair */
/*                combinations. */
/*                The integrals are ordered in the batch the following */
/*                way (first index varying fastest): */
/*                    batch (nxyz1,nxyz2,nxyz3,nxyz4,kl,ij) */
/*                where ij and kl indicates alpha exponent pairs */
/*                defining the present block. */
/*                The present routine evaluates batches of the type: */
/*                            sppp , pspp , ppsp , ppps */
/*                The cartesian primitive integrals are evaluated using */
/*                the auxillary functions technique, described by */
/*                V.R. Saunders, "An Introduction to Molecular Integral */
/*                Evaluation" in "Computational Techniques in Quantum */
/*                Chemistry and Molecular Physics" edited by */
/*                GHF Diercksen, BT Sutcliff and A Veillard, D. Reidel */
/*                Publ. Comp. Dordrecht (1975), p. 347. */
/*                The cartesian primitive integrals are each evaluated */
/*                explicitely using common multipliers if possible and */
/*                avoiding array addresses. */
/*                  Input: */
/*                    NBATCH       =  size of the primitive cartesian */
/*                                    sppp/pspp/ppsp/ppps integral batch */
/*                    ATOMIC       =  indicates, if purely atomic */
/*                                    integrals will be evaluated */
/*                    ATOM12(34)   =  indicates, if centers 1 and 2 */
/*                                    (3 and 4) coincide */
/*                    MIJ(KL)      =  current # of ij (kl) primitive */
/*                                    index pairs corresponding to */
/*                                    the contracted shell pairs 1,2 */
/*                                    (3,4) */
/*                    NIJ          =  total # of ij primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair 1,2 */
/*                    NIJBEG(END)  =  first(last) ij primitive index */
/*                                    defining the ij block */
/*                    NKL          =  total # of kl primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair 3,4 */
/*                    NKLBEG(END)  =  first(last) kl primitive index */
/*                                    defining the kl block */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for contraction shells x = 1,2,3,4 */
/*                    SHELLx       =  the shell type for contraction */
/*                                    shells x = 1,3,P=1+2 */
/*                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers */
/*                                    x = 1,2,3,4 */
/*                    Xxx,Yxx,Zxx  =  the x,y,z-coordinate differences */
/*                                    between centers xx = 12 and 34 */
/*                    ALPHAx       =  the primitive exponents for */
/*                                    contraction shells x = 1,2,3,4 */
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
/*                                    x = 1,2,3,4 */
/*                    NORMx        =  the normalization factors due to */
/*                                    the primitive exponents for the */
/*                                    contraction shells x = 1,2,3,4 */
/*                    RHO12(34)    =  the complete set of NIJ (NKL) */
/*                                    exponential prefactors between */
/*                                    contraction shells 1 and 2 */
/*                                    (3 and 4) */
/*                    P            =  will hold current MIJ exponent */
/*                                    sums for contraction shells 1 */
/*                                    and 2 */
/*                    Px           =  will hold current MIJ coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers P=1+2 */
/*                    SCALEP       =  will hold current MIJ values of */
/*                                    scaling factors related to point P */
/*                    Q            =  will hold current MKL exponent */
/*                                    sums for contraction shells 3 */
/*                                    and 4 */
/*                    Qx           =  will hold current MKL coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers Q=3+4 */
/*                    SCALEQ       =  will hold current MKL values of */
/*                                    scaling factors related to point Q */
/*                  Output: */
/*                    BATCH        =  current batch of primitive */
/*                                    cartesian sppp/pspp/ppsp/ppps */
/*                                    integrals */
/* ------------------------------------------------------------------------ */
int erd__sppp_pcgto_block (int nbatch, int atomic, int atom12, int atom34,
                           int mij, int mkl,
                           int nij, int nijbeg, int nijend,
                           int nkl, int nklbeg, int nklend,
                           int npgto1, int npgto2,
                           int npgto3, int npgto4,
                           int shell1, int shell3, int shellp,
                           double x1, double y1, double z1,
                           double x2, double y2, double z2,
                           double x3, double y3, double z3,
                           double x4, double y4, double z4,
                           double x12, double y12, double z12,
                           double x34, double y34, double z34,
                           double *alpha1, double *alpha2,
                           double *alpha3, double *alpha4,
                           double *ftable, int mgrid, int ngrid,
                           double tmax, double tstep, double tvstep,
                           int *prim1, int *prim2,
                           int *prim3, int *prim4,
                           double *norm1, double *norm2,
                           double *norm3, double *norm4,
                           double *rho12, double *rho34,
                           double *p, double *px,
                           double *py, double *pz, double *scalep,
                           double *q, double *qx,
                           double *qy, double *qz,
                           double *scaleq, double *batch)
{
    int ftable_dim1, ftable_offset;

    /* Local variables */
    double a;
    double b;
    double c;
    double d;
    double e;
    double f;
    double g;
    double h;
    int i;
    int j;
    int k;
    int l;
    int m;
    double r;
    double t;
    double f0;
    double f1;
    double f2;
    double f3;
    double u0;
    double u1;
    double u2;
    double u3;
    double u4;
    double u5;
    double aa;
    double bb;
    int ij;
    int kl;
    double pqx;
    double pqy;
    double pqz;
    double exp1;
    double exp2;
    double exp3;
    double exp4;
    int tcase;
    double pval;
    double qval;
    double pinv;
    double qinv;
    double tinv;
    double gxxx;
    double gxxy;
    double gxxz;
    double gxyx;
    double gxyy;
    double gxyz;
    double gxzx;
    double gxzy;
    double gxzz;
    double gyxx;
    double gyxy;
    double gyxz;
    double gyyx;
    double gyyy;
    double gyyz;
    double gyzx;
    double gyzy;
    double gyzz;
    double gzxx;
    double gzxy;
    double gzxz;
    double gzyx;
    double gzyy;
    double gzyz;
    double gzzx;
    double gzzy;
    double gzzz;
#ifndef ERD_TABLE_FREE_BOYS_FUNCTIONS
    double t2inv;
#endif
    double xppp1;
    double xppp2;
    double xppp3;
    double xspp1;
    double xpsp1;
    double xpps1;
    double xssp1;
    double xsps1;
    double xpss1;
    double xssp2;
    double xsps2;
    double xpss2;
    double xspp2;
    double xpsp2;
    double scale;
    double xpps2;
    double xspp3;
    double xpsp3;
    double xpps3;
    double xppp4;
    double yssp1;
    double ysps1;
    double ypss1;
    double yssp2;
    double ysps2;
    double ypss2;
    double yspp1;
    double ypsp1;
#ifndef ERD_TABLE_FREE_BOYS_FUNCTIONS
    int tgrid;
#endif
    double ypps1;
    double yspp2;
    double ypsp2;
    double ypps2;
    double yspp3;
    double ypsp3;
    double ypps3;
    double yppp1;
    double yppp2;
    double yppp3;
    double yppp4;
    double zssp1;
    double zsps1;
    double zpss1;
    double zssp2;
    double zsps2;
    double pxval;
    double pyval;
    double pzval;
    double qxval;
    double qyval;
    double qzval;
    double pxsub;
    double pysub;
    double pzsub;
    double qxsub;
    double qysub;
#ifndef ERD_TABLE_FREE_BOYS_FUNCTIONS
    double delta1;
    double delta2;
    double delta3;
    double delta4;
    double delta5;
    double delta6;
#endif
    double qzsub;
    double zpss2;
    double zspp1;
    double zpsp1;
    double zpps1;
    double zspp2;
    double zpsp2;
    double zpps2;
    double zspp3;
    double zpsp3;
    double zpps3;
    double zppp1;
    double zppp2;
    double zppp3;
    double zppp4;
    double pscale;
    double pqpinv;
    double pqmult;
    double pqplus;
    double rnpqsq;

    --batch;
    --scalep;
    --pz;
    --py;
    --px;
    --p;
    --prim2;
    --prim1;
    --scaleq;
    --qz;
    --qy;
    --qx;
    --q;
    --prim4;
    --prim3;
    --rho12;
    --rho34;
    --norm1;
    --alpha1;
    --norm2;
    --alpha2;
    --norm3;
    --alpha3;
    --norm4;
    --alpha4;
    ftable_dim1 = mgrid - 0 + 1;
    ftable_offset = 0 + ftable_dim1 * 0;
    ftable -= ftable_offset;

    if (atomic)
    {
        return 0;
    }
        
    tcase = 5;
    {
        m = 0;
        for (ij = nijbeg; ij <= nijend; ++ij)
        {
            ++m;
            i = prim1[m];
            j = prim2[m];
            exp1 = alpha1[i];
            exp2 = alpha2[j];
            pval = exp1 + exp2;
            p[m] = pval;
            pval = exp1 / pval;
            px[m] = pval * x12 + x2;
            py[m] = pval * y12 + y2;
            pz[m] = pval * z12 + z2;
            scalep[m] = norm1[i] * norm2[j] * rho12[ij];
        }
    }

    {
        m = 0;
        for (kl = nklbeg; kl <= nklend; ++kl)
        {
            ++m;
            k = prim3[m];
            l = prim4[m];
            exp3 = alpha3[k];
            exp4 = alpha4[l];
            qval = exp3 + exp4;
            q[m] = qval;
            qval = exp3 / qval;
            qx[m] = qval * x34 + x4;
            qy[m] = qval * y34 + y4;
            qz[m] = qval * z34 + z4;
            scaleq[m] = norm3[k] * norm4[l] * rho34[kl];
        }
    }

/*             ...jump according to where the s-type function is */
/*                located and the type of 'K4' loop: */
/*                 SHELLP CASE  |        Integral center type */
/*                --------------|------------------------------------ */
/*                    1     1   |  (AA|AA)  atomic     sppp and pspp */
/*                    1     2   |  (AA|CC)  2-center   sppp and pspp */
/*                    1     3   |  (AB|CC)  3-center   sppp and pspp */
/*                    1     4   |  (AA|CD)  3-center   sppp and pspp */
/*                    1     5   |  (AB|CD)  4-center   sppp and pspp */
/*                    2     1   |  (AA|AA)  atomic     ppsp and ppps */
/*                    2     2   |  (AA|CC)  2-center   ppsp and ppps */
/*                    2     3   |  (AB|CC)  3-center   ppsp and ppps */
/*                    2     4   |  (AA|CD)  3-center   ppsp and ppps */
/*                    2     5   |  (AB|CD)  4-center   ppsp and ppps */
    switch ((shellp - 1) * 5 + tcase)
    {
    case 5:
        goto L15;
    case 10:
        goto L25;
    }

  L15:
    if (shell1 == 1)
    {
        pxsub = x1;
        pysub = y1;
        pzsub = z1;
    }
    else
    {
        pxsub = x2;
        pysub = y2;
        pzsub = z2;
    }
    m = 0;
    for (ij = 1; ij <= mij; ++ij)
    {
        pval = p[ij];
        pxval = px[ij];
        pyval = py[ij];
        pzval = pz[ij];
        pscale = scalep[ij];
        xpss1 = pxval - pxsub;
        ypss1 = pyval - pysub;
        zpss1 = pzval - pzsub;
        for (kl = 1; kl <= mkl; ++kl)
        {
            qval = q[kl];
            qinv = 1. / qval;
            qxval = qx[kl];
            qyval = qy[kl];
            qzval = qz[kl];
            pqmult = pval * qval;
            pqplus = pval + qval;
            pqpinv = 1. / pqplus;
            pqx = pxval - qxval;
            pqy = pyval - qyval;
            pqz = pzval - qzval;
            t = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult * pqpinv;
            scale = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
#ifdef ERD_TABLE_FREE_BOYS_FUNCTIONS
            f0 = scale * boys0(t);
            f1 = scale * boys1(t);
            f2 = scale * boys2(t);
            f3 = scale * boys3(t);
#else
            if (t <= tmax)
            {
                tgrid = (int) (t * tvstep + .5);
                delta1 = tgrid * tstep - t;
                delta2 = delta1 * .5;
                delta3 = delta1 * .333333333333333;
                delta4 = delta2 * .5;
                delta5 = delta1 * .2;
                delta6 = delta3 * .5;
                f0 = (((((ftable[tgrid * ftable_dim1 + 6] * delta6 +
                          ftable[tgrid * ftable_dim1 + 5]) * delta5 +
                         ftable[tgrid * ftable_dim1 + 4]) * delta4 +
                        ftable[tgrid * ftable_dim1 + 3]) * delta3 +
                       ftable[tgrid * ftable_dim1 + 2]) * delta2 +
                      ftable[tgrid * ftable_dim1 + 1]) * delta1 +
                    ftable[tgrid * ftable_dim1];
                f1 = (((((ftable[tgrid * ftable_dim1 + 7] * delta6 +
                          ftable[tgrid * ftable_dim1 + 6]) * delta5 +
                         ftable[tgrid * ftable_dim1 + 5]) * delta4 +
                        ftable[tgrid * ftable_dim1 + 4]) * delta3 +
                       ftable[tgrid * ftable_dim1 + 3]) * delta2 +
                      ftable[tgrid * ftable_dim1 + 2]) * delta1 +
                    ftable[tgrid * ftable_dim1 + 1];
                f2 = (((((ftable[tgrid * ftable_dim1 + 8] * delta6 +
                          ftable[tgrid * ftable_dim1 + 7]) * delta5 +
                         ftable[tgrid * ftable_dim1 + 6]) * delta4 +
                        ftable[tgrid * ftable_dim1 + 5]) * delta3 +
                       ftable[tgrid * ftable_dim1 + 4]) * delta2 +
                      ftable[tgrid * ftable_dim1 + 3]) * delta1 +
                    ftable[tgrid * ftable_dim1 + 2];
                f3 = (((((ftable[tgrid * ftable_dim1 + 9] * delta6 +
                          ftable[tgrid * ftable_dim1 + 8]) * delta5 +
                         ftable[tgrid * ftable_dim1 + 7]) * delta4 +
                        ftable[tgrid * ftable_dim1 + 6]) * delta3 +
                       ftable[tgrid * ftable_dim1 + 5]) * delta2 +
                      ftable[tgrid * ftable_dim1 + 4]) * delta1 +
                    ftable[tgrid * ftable_dim1 + 3];
                f0 = scale * f0;
                f1 = scale * f1;
                f2 = scale * f2;
                f3 = scale * f3;
            }
            else
            {
                tinv = 1. / t;
                t2inv = tinv * .5;
                f0 = scale * .5 * sqrt (tinv * 3.141592653589793);
                f1 = t2inv * f0;
                f2 = t2inv * 3. * f1;
                f3 = t2inv * 5. * f2;
            }
#endif
            u0 = pval * pqpinv;
            u1 = -qval * pqpinv;
            u2 = pqpinv * .5;
            u3 = u2 + pqpinv;
            u4 = qinv * .5;
            u5 = u0 * u4;


/*             ...the X-terms. */
            xssp1 = qxval - x4;
            xsps1 = qxval - x3;
            xssp2 = pqx * u0;
            xpss2 = pqx * u1;
            a = xsps1 + xssp1;
            b = xpss1 * xssp2 + u2;
            xspp1 = xsps1 * xssp1 + u4;
            xspp2 = a * xssp2 - u5;
            xspp3 = xssp2 * xssp2;
            xpsp1 = xpss1 * xssp1;
            xpsp2 = xssp1 * xpss2 + b;
            xpsp3 = xssp2 * xpss2;
            xpps1 = xpss1 * xsps1;
            xpps2 = xsps1 * xpss2 + b;
            xppp1 = xpss1 * xspp1;
            xppp2 = xpss1 * xspp2 + xspp1 * xpss2 + a * u2;
            xppp3 = xpss1 * xspp3 + a * xpsp3 + u3 * xssp2;
            xppp4 = xpss2 * xspp3;


/*             ...the Y-terms. */
            yssp1 = qyval - y4;
            ysps1 = qyval - y3;
            yssp2 = pqy * u0;
            ypss2 = pqy * u1;
            a = ysps1 + yssp1;
            b = ypss1 * yssp2 + u2;
            yspp1 = ysps1 * yssp1 + u4;
            yspp2 = a * yssp2 - u5;
            yspp3 = yssp2 * yssp2;
            ypsp1 = ypss1 * yssp1;
            ypsp2 = yssp1 * ypss2 + b;
            ypsp3 = yssp2 * ypss2;
            ypps1 = ypss1 * ysps1;
            ypps2 = ysps1 * ypss2 + b;
            yppp1 = ypss1 * yspp1;
            yppp2 = ypss1 * yspp2 + yspp1 * ypss2 + a * u2;
            yppp3 = ypss1 * yspp3 + a * ypsp3 + u3 * yssp2;
            yppp4 = ypss2 * yspp3;


/*             ...the Z-terms. */
            zssp1 = qzval - z4;
            zsps1 = qzval - z3;
            zssp2 = pqz * u0;
            zpss2 = pqz * u1;
            a = zsps1 + zssp1;
            b = zpss1 * zssp2 + u2;
            zspp1 = zsps1 * zssp1 + u4;
            zspp2 = a * zssp2 - u5;
            zspp3 = zssp2 * zssp2;
            zpsp1 = zpss1 * zssp1;
            zpsp2 = zssp1 * zpss2 + b;
            zpsp3 = zssp2 * zpss2;
            zpps1 = zpss1 * zsps1;
            zpps2 = zsps1 * zpss2 + b;
            zppp1 = zpss1 * zspp1;
            zppp2 = zpss1 * zspp2 + zspp1 * zpss2 + a * u2;
            zppp3 = zpss1 * zspp3 + a * zpsp3 + u3 * zssp2;
            zppp4 = zpss2 * zspp3;


/*             ...assemble the 4-center (AB|CD) type integrals. */
            gxxx = xppp1 * f0 + xppp2 * f1 + xppp3 * f2 + xppp4 * f3;
            gyyy = yppp1 * f0 + yppp2 * f1 + yppp3 * f2 + yppp4 * f3;
            gzzz = zppp1 * f0 + zppp2 * f1 + zppp3 * f2 + zppp4 * f3;
            aa = xpsp3 * f2;
            bb = xpsp3 * f3;
            a = xpps1 * f0 + xpps2 * f1 + aa;
            b = xpps1 * f1 + xpps2 * f2 + bb;
            c = xpsp1 * f0 + xpsp2 * f1 + aa;
            d = xpsp1 * f1 + xpsp2 * f2 + bb;
            e = xspp1 * f0 + xspp2 * f1 + xspp3 * f2;
            f = xspp1 * f1 + xspp2 * f2 + xspp3 * f3;
            gxxy = yssp1 * a + yssp2 * b;
            gxxz = zssp1 * a + zssp2 * b;
            gxyx = ysps1 * c + yssp2 * d;
            gxzx = zsps1 * c + zssp2 * d;
            gyxx = ypss1 * e + ypss2 * f;
            gzxx = zpss1 * e + zpss2 * f;
            aa = ypsp3 * f2;
            bb = ypsp3 * f3;
            a = ypps1 * f0 + ypps2 * f1 + aa;
            b = ypps1 * f1 + ypps2 * f2 + bb;
            c = ypsp1 * f0 + ypsp2 * f1 + aa;
            d = ypsp1 * f1 + ypsp2 * f2 + bb;
            e = yspp1 * f0 + yspp2 * f1 + yspp3 * f2;
            f = yspp1 * f1 + yspp2 * f2 + yspp3 * f3;
            gyyx = xssp1 * a + xssp2 * b;
            gyyz = zssp1 * a + zssp2 * b;
            gyxy = xsps1 * c + xssp2 * d;
            gyzy = zsps1 * c + zssp2 * d;
            gxyy = xpss1 * e + xpss2 * f;
            gzyy = zpss1 * e + zpss2 * f;
            aa = zpsp3 * f2;
            bb = zpsp3 * f3;
            a = zpps1 * f0 + zpps2 * f1 + aa;
            b = zpps1 * f1 + zpps2 * f2 + bb;
            c = zpsp1 * f0 + zpsp2 * f1 + aa;
            d = zpsp1 * f1 + zpsp2 * f2 + bb;
            e = zspp1 * f0 + zspp2 * f1 + zspp3 * f2;
            f = zspp1 * f1 + zspp2 * f2 + zspp3 * f3;
            gzzx = xssp1 * a + xssp2 * b;
            gzzy = yssp1 * a + yssp2 * b;
            gzxz = xsps1 * c + xssp2 * d;
            gzyz = ysps1 * c + yssp2 * d;
            gxzz = xpss1 * e + xpss2 * f;
            gyzz = ypss1 * e + ypss2 * f;
            a = xpss1 * f0 + xpss2 * f1;
            b = xpss1 * f1 + xpss2 * f2;
            c = xpss1 * f2 + xpss2 * f3;
            d = ypss1 * f0 + ypss2 * f1;
            e = ypss1 * f1 + ypss2 * f2;
            f = ypss1 * f2 + ypss2 * f3;
            g = zpss1 * f0 + zpss2 * f1;
            h = zpss1 * f1 + zpss2 * f2;
            r = zpss1 * f2 + zpss2 * f3;
            gxyz = ysps1 * (zssp1 * a + zssp2 * b) + yssp2 * (zssp1 * b +
                                                              zssp2 * c);
            gxzy = zsps1 * (yssp1 * a + yssp2 * b) + zssp2 * (yssp1 * b +
                                                              yssp2 * c);
            gyxz = xsps1 * (zssp1 * d + zssp2 * e) + xssp2 * (zssp1 * e +
                                                                zssp2 * f);
            gyzx = zsps1 * (xssp1 * d + xssp2 * e) + zssp2 * (xssp1 * e +
                                                                xssp2 * f);
            gzxy = xsps1 * (yssp1 * g + yssp2 * h) + xssp2 * (yssp1 * h +
                                                                yssp2 * r);
            gzyx = ysps1 * (xssp1 * g + xssp2 * h) + yssp2 * (xssp1 * h +
                                                                xssp2 * r);
            batch[m + 1] = gxxx;
            batch[m + 2] = gyxx;
            batch[m + 3] = gzxx;
            batch[m + 4] = gxyx;
            batch[m + 5] = gyyx;
            batch[m + 6] = gzyx;
            batch[m + 7] = gxzx;
            batch[m + 8] = gyzx;
            batch[m + 9] = gzzx;
            batch[m + 10] = gxxy;
            batch[m + 11] = gyxy;
            batch[m + 12] = gzxy;
            batch[m + 13] = gxyy;
            batch[m + 14] = gyyy;
            batch[m + 15] = gzyy;
            batch[m + 16] = gxzy;
            batch[m + 17] = gyzy;
            batch[m + 18] = gzzy;
            batch[m + 19] = gxxz;
            batch[m + 20] = gyxz;
            batch[m + 21] = gzxz;
            batch[m + 22] = gxyz;
            batch[m + 23] = gyyz;
            batch[m + 24] = gzyz;
            batch[m + 25] = gxzz;
            batch[m + 26] = gyzz;
            batch[m + 27] = gzzz;
            m += 27;
        }
    }
    return 0;

  L25:
    if (shell3 == 1)
    {
        qxsub = x3;
        qysub = y3;
        qzsub = z3;
    }
    else
    {
        qxsub = x4;
        qysub = y4;
        qzsub = z4;
    }
    m = 0;
    for (ij = 1; ij <= mij; ++ij)
    {
        pval = p[ij];
        pinv = 1. / pval;
        pxval = px[ij];
        pyval = py[ij];
        pzval = pz[ij];
        pscale = scalep[ij];
        u4 = pinv * .5;
        xsps1 = pxval - x2;
        xpss1 = pxval - x1;
        ysps1 = pyval - y2;
        ypss1 = pyval - y1;
        zsps1 = pzval - z2;
        zpss1 = pzval - z1;
        for (kl = 1; kl <= mkl; ++kl)
        {
            qval = q[kl];
            qxval = qx[kl];
            qyval = qy[kl];
            qzval = qz[kl];
            pqmult = pval * qval;
            pqplus = pval + qval;
            pqpinv = 1. / pqplus;
            pqx = pxval - qxval;
            pqy = pyval - qyval;
            pqz = pzval - qzval;
            t = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult * pqpinv;
            scale = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
#ifdef ERD_TABLE_FREE_BOYS_FUNCTIONS
            f0 = scale * boys0(t);
            f1 = scale * boys1(t);
            f2 = scale * boys2(t);
            f3 = scale * boys3(t);
#else
            if (t <= tmax)
            {
                tgrid = (int) (t * tvstep + .5);
                delta1 = tgrid * tstep - t;
                delta2 = delta1 * .5;
                delta3 = delta1 * .333333333333333;
                delta4 = delta2 * .5;
                delta5 = delta1 * .2;
                delta6 = delta3 * .5;
                f0 = (((((ftable[tgrid * ftable_dim1 + 6] * delta6 +
                          ftable[tgrid * ftable_dim1 + 5]) * delta5 +
                         ftable[tgrid * ftable_dim1 + 4]) * delta4 +
                        ftable[tgrid * ftable_dim1 + 3]) * delta3 +
                       ftable[tgrid * ftable_dim1 + 2]) * delta2 +
                      ftable[tgrid * ftable_dim1 + 1]) * delta1 +
                    ftable[tgrid * ftable_dim1];
                f1 = (((((ftable[tgrid * ftable_dim1 + 7] * delta6 +
                          ftable[tgrid * ftable_dim1 + 6]) * delta5 +
                         ftable[tgrid * ftable_dim1 + 5]) * delta4 +
                        ftable[tgrid * ftable_dim1 + 4]) * delta3 +
                       ftable[tgrid * ftable_dim1 + 3]) * delta2 +
                      ftable[tgrid * ftable_dim1 + 2]) * delta1 +
                    ftable[tgrid * ftable_dim1 + 1];
                f2 = (((((ftable[tgrid * ftable_dim1 + 8] * delta6 +
                          ftable[tgrid * ftable_dim1 + 7]) * delta5 +
                         ftable[tgrid * ftable_dim1 + 6]) * delta4 +
                        ftable[tgrid * ftable_dim1 + 5]) * delta3 +
                       ftable[tgrid * ftable_dim1 + 4]) * delta2 +
                      ftable[tgrid * ftable_dim1 + 3]) * delta1 +
                    ftable[tgrid * ftable_dim1 + 2];
                f3 = (((((ftable[tgrid * ftable_dim1 + 9] * delta6 +
                          ftable[tgrid * ftable_dim1 + 8]) * delta5 +
                         ftable[tgrid * ftable_dim1 + 7]) * delta4 +
                        ftable[tgrid * ftable_dim1 + 6]) * delta3 +
                       ftable[tgrid * ftable_dim1 + 5]) * delta2 +
                      ftable[tgrid * ftable_dim1 + 4]) * delta1 +
                    ftable[tgrid * ftable_dim1 + 3];
                f0 = scale * f0;
                f1 = scale * f1;
                f2 = scale * f2;
                f3 = scale * f3;
            }
            else
            {
                tinv = 1. / t;
                t2inv = tinv * .5;
                f0 = scale * .5 * sqrt (tinv * 3.141592653589793);
                f1 = t2inv * f0;
                f2 = t2inv * 3. * f1;
                f3 = t2inv * 5. * f2;
            }
#endif
            u0 = pval * pqpinv;
            u1 = -qval * pqpinv;
            u2 = pqpinv * .5;
            u3 = u2 + pqpinv;
            u5 = u1 * u4;


/*             ...the X-terms. */
            xssp1 = qxval - qxsub;
            xssp2 = pqx * u0;
            xsps2 = pqx * u1;
            a = xpss1 + xsps1;
            b = xssp1 * xsps2 + u2;
            xspp1 = xsps1 * xssp1;
            xspp2 = xsps1 * xssp2 + b;
            xspp3 = xssp2 * xsps2;
            xpsp1 = xpss1 * xssp1;
            xpsp2 = xpss1 * xssp2 + b;
            xpps1 = xpss1 * xsps1 + u4;
            xpps2 = a * xsps2 + u5;
            xpps3 = xsps2 * xsps2;
            xppp1 = xssp1 * xpps1;
            xppp2 = xssp1 * xpps2 + xpps1 * xssp2 + a * u2;
            xppp3 = xssp1 * xpps3 + a * xspp3 + u3 * xsps2;
            xppp4 = xssp2 * xpps3;


/*             ...the Y-terms. */
            yssp1 = qyval - qysub;
            yssp2 = pqy * u0;
            ysps2 = pqy * u1;
            a = ypss1 + ysps1;
            b = yssp1 * ysps2 + u2;
            yspp1 = ysps1 * yssp1;
            yspp2 = ysps1 * yssp2 + b;
            yspp3 = yssp2 * ysps2;
            ypsp1 = ypss1 * yssp1;
            ypsp2 = ypss1 * yssp2 + b;
            ypps1 = ypss1 * ysps1 + u4;
            ypps2 = a * ysps2 + u5;
            ypps3 = ysps2 * ysps2;
            yppp1 = yssp1 * ypps1;
            yppp2 = yssp1 * ypps2 + ypps1 * yssp2 + a * u2;
            yppp3 = yssp1 * ypps3 + a * yspp3 + u3 * ysps2;
            yppp4 = yssp2 * ypps3;


/*             ...the Z-terms. */
            zssp1 = qzval - qzsub;
            zssp2 = pqz * u0;
            zsps2 = pqz * u1;
            a = zpss1 + zsps1;
            b = zssp1 * zsps2 + u2;
            zspp1 = zsps1 * zssp1;
            zspp2 = zsps1 * zssp2 + b;
            zspp3 = zssp2 * zsps2;
            zpsp1 = zpss1 * zssp1;
            zpsp2 = zpss1 * zssp2 + b;
            zpps1 = zpss1 * zsps1 + u4;
            zpps2 = a * zsps2 + u5;
            zpps3 = zsps2 * zsps2;
            zppp1 = zssp1 * zpps1;
            zppp2 = zssp1 * zpps2 + zpps1 * zssp2 + a * u2;
            zppp3 = zssp1 * zpps3 + a * zspp3 + u3 * zsps2;
            zppp4 = zssp2 * zpps3;


/*             ...assemble the 4-center (AB|CD) type integrals. */
            gxxx = xppp1 * f0 + xppp2 * f1 + xppp3 * f2 + xppp4 * f3;
            gyyy = yppp1 * f0 + yppp2 * f1 + yppp3 * f2 + yppp4 * f3;
            gzzz = zppp1 * f0 + zppp2 * f1 + zppp3 * f2 + zppp4 * f3;
            aa = xspp3 * f2;
            bb = xspp3 * f3;
            a = xpps1 * f0 + xpps2 * f1 + xpps3 * f2;
            b = xpps1 * f1 + xpps2 * f2 + xpps3 * f3;
            c = xpsp1 * f0 + xpsp2 * f1 + aa;
            d = xpsp1 * f1 + xpsp2 * f2 + bb;
            e = xspp1 * f0 + xspp2 * f1 + aa;
            f = xspp1 * f1 + xspp2 * f2 + bb;
            gxxy = yssp1 * a + yssp2 * b;
            gxxz = zssp1 * a + zssp2 * b;
            gxyx = ysps1 * c + ysps2 * d;
            gxzx = zsps1 * c + zsps2 * d;
            gyxx = ypss1 * e + ysps2 * f;
            gzxx = zpss1 * e + zsps2 * f;
            aa = yspp3 * f2;
            bb = yspp3 * f3;
            a = ypps1 * f0 + ypps2 * f1 + ypps3 * f2;
            b = ypps1 * f1 + ypps2 * f2 + ypps3 * f3;
            c = ypsp1 * f0 + ypsp2 * f1 + aa;
            d = ypsp1 * f1 + ypsp2 * f2 + bb;
            e = yspp1 * f0 + yspp2 * f1 + aa;
            f = yspp1 * f1 + yspp2 * f2 + bb;
            gyyx = xssp1 * a + xssp2 * b;
            gyyz = zssp1 * a + zssp2 * b;
            gyxy = xsps1 * c + xsps2 * d;
            gyzy = zsps1 * c + zsps2 * d;
            gxyy = xpss1 * e + xsps2 * f;
            gzyy = zpss1 * e + zsps2 * f;
            aa = zspp3 * f2;
            bb = zspp3 * f3;
            a = zpps1 * f0 + zpps2 * f1 + zpps3 * f2;
            b = zpps1 * f1 + zpps2 * f2 + zpps3 * f3;
            c = zpsp1 * f0 + zpsp2 * f1 + aa;
            d = zpsp1 * f1 + zpsp2 * f2 + bb;
            e = zspp1 * f0 + zspp2 * f1 + aa;
            f = zspp1 * f1 + zspp2 * f2 + bb;
            gzzx = xssp1 * a + xssp2 * b;
            gzzy = yssp1 * a + yssp2 * b;
            gzxz = xsps1 * c + xsps2 * d;
            gzyz = ysps1 * c + ysps2 * d;
            gxzz = xpss1 * e + xsps2 * f;
            gyzz = ypss1 * e + ysps2 * f;
            a = xpss1 * f0 + xsps2 * f1;
            b = xpss1 * f1 + xsps2 * f2;
            c = xpss1 * f2 + xsps2 * f3;
            d = ypss1 * f0 + ysps2 * f1;
            e = ypss1 * f1 + ysps2 * f2;
            f = ypss1 * f2 + ysps2 * f3;
            g = zpss1 * f0 + zsps2 * f1;
            h = zpss1 * f1 + zsps2 * f2;
            r = zpss1 * f2 + zsps2 * f3;
            gxyz = ysps1 * (zssp1 * a + zssp2 * b) + ysps2 * (zssp1 * b +
                                                              zssp2 * c);
            gxzy = zsps1 * (yssp1 * a + yssp2 * b) + zsps2 * (yssp1 * b +
                                                              yssp2 * c);
            gyxz = xsps1 * (zssp1 * d + zssp2 * e) + xsps2 * (zssp1 * e +
                                                                zssp2 * f);
            gyzx = zsps1 * (xssp1 * d + xssp2 * e) + zsps2 * (xssp1 * e +
                                                                xssp2 * f);
            gzxy = xsps1 * (yssp1 * g + yssp2 * h) + xsps2 * (yssp1 * h +
                                                                yssp2 * r);
            gzyx = ysps1 * (xssp1 * g + xssp2 * h) + ysps2 * (xssp1 * h +
                                                                xssp2 * r);
            batch[m + 1] = gxxx;
            batch[m + 2] = gyxx;
            batch[m + 3] = gzxx;
            batch[m + 4] = gxyx;
            batch[m + 5] = gyyx;
            batch[m + 6] = gzyx;
            batch[m + 7] = gxzx;
            batch[m + 8] = gyzx;
            batch[m + 9] = gzzx;
            batch[m + 10] = gxxy;
            batch[m + 11] = gyxy;
            batch[m + 12] = gzxy;
            batch[m + 13] = gxyy;
            batch[m + 14] = gyyy;
            batch[m + 15] = gzyy;
            batch[m + 16] = gxzy;
            batch[m + 17] = gyzy;
            batch[m + 18] = gzzy;
            batch[m + 19] = gxxz;
            batch[m + 20] = gyxz;
            batch[m + 21] = gzxz;
            batch[m + 22] = gxyz;
            batch[m + 23] = gyyz;
            batch[m + 24] = gzyz;
            batch[m + 25] = gxzz;
            batch[m + 26] = gyzz;
            batch[m + 27] = gzzz;
            m += 27;
        }
    }

    return 0;
}