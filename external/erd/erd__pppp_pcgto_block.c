#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__PPPP_PCGTO_BLOCK */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation is designed to provide ultrafast block */
/*                evaluation of a batch of normalized electron repulsion */
/*                integrals between p-shell primitive spherical gaussian */
/*                type orbitals. */
/*                A batch is defined here as containing all possible */
/*                integrals;double that is its dimension is determined by */
/*                the total number of primitive functions (here = 81) */
/*                times the total number of ij and kl exponent pair */
/*                combinations. */
/*                The integrals are ordered in the batch the following */
/*                way (first index varying fastest): */
/*                    batch (nxyz1,nxyz2,nxyz3,nxyz4,kl,ij) */
/*                where ij and kl indicates alpha exponent pairs */
/*                defining the present block. */
/*                The present routine evaluates batches of the type: */
/*                                         pppp */
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
/*                                    pppp integral batch */
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
/*                                    cartesian pppp integrals */
/* ------------------------------------------------------------------------ */
int erd__pppp_pcgto_block (int nij, int nkl,
                           double x1, double y1, double z1,
                           double x2, double y2, double z2,
                           double x3, double y3, double z3,
                           double x4, double y4, double z4,
                           double *alpha1, double *alpha2,
                           double *alpha3, double *alpha4,
                           double *ftable, int mgrid,
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

    double zpppp4;
    double zpppp5;
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
    double f4;
    double u0;
    double u1;
    double u2;
    double u3;
    double u4;
    double u5;
    double u6;
    double u7;
    double u8;
    double aa;
    double bb;
    double cc;
    double dd;
    double ee;
    double ff;
    double gg;
    double hh;
    int ij;
    int kl;
    double rr;
    double pqx;
    double pqy;
    double pqz;
    double exp1;
    double exp2;
    double exp3;
    double exp4;
    double pval;
    double qval;
    double tinv;
    double t2inv;
    double scale;
    int tgrid;
    double pxval;
    double pyval;
    double pzval;
    double qxval;
    double qyval;
    double qzval;
    double delta1;
    double delta2;
    double delta3;
    double delta4;
    double delta5;
    double delta6;
    double gxxxx;
    double gxxxy;
    double gxxxz;
    double gxxyx;
    double gxxyy;
    double gxxyz;
    double gxxzx;
    double gxxzy;
    double gxxzz;
    double gxyxx;
    double gxyxy;
    double gxyxz;
    double gxyyx;
    double gxyyy;
    double gxyyz;
    double gxyzx;
    double gxyzy;
    double gxyzz;
    double gxzxx;
    double gxzxy;
    double gxzxz;
    double gxzyx;
    double gxzyy;
    double gxzyz;
    double gxzzx;
    double gxzzy;
    double gxzzz;
    double gyxxx;
    double gyxxy;
    double gyxxz;
    double gyxyx;
    double gyxyy;
    double gyxyz;
    double gyxzx;
    double gyxzy;
    double gyxzz;
    double gyyxx;
    double gyyxy;
    double gyyxz;
    double gyyyx;
    double gyyyy;
    double gyyyz;
    double gyyzx;
    double gyyzy;
    double gyyzz;
    double gyzxx;
    double gyzxy;
    double gyzxz;
    double gyzyx;
    double gyzyy;
    double gyzyz;
    double gyzzx;
    double gyzzy;
    double gyzzz;
    double gzxxx;
    double gzxxy;
    double gzxxz;
    double gzxyx;
    double gzxyy;
    double gzxyz;
    double gzxzx;
    double gzxzy;
    double gzxzz;
    double gzyxx;
    double gzyxy;
    double gzyxz;
    double gzyyx;
    double gzyyy;
    double gzyyz;
    double gzyzx;
    double gzyzy;
    double gzyzz;
    double gzzxx;
    double gzzxy;
    double gzzxz;
    double gzzyx;
    double gzzyy;
    double gzzyz;
    double gzzzx;
    double gzzzy;
    double gzzzz;
    double pqplus;
    double pqmult;
    double pqpinv;
    double pscale;
    double xsssp1;
    double xssps1;
    double xspss1;
    double xpsss1;
    double xsspp1;
    double xspsp1;
    double xpssp1;
    double xspps1;
    double xpsps1;
    double xppss1;
    double xsppp1;
    double xpspp1;
    double xppsp1;
    double xppps1;
    double xpppp1;
    double xsssp2;
    double xspss2;
    double xsspp2;
    double xspsp2;
    double xpssp2;
    double xspps2;
    double xpsps2;
    double xppss2;
    double xsppp2;
    double xpspp2;
    double xppsp2;
    double xppps2;
    double xpppp2;
    double xsspp3;
    double xspsp3;
    double xppss3;
    double xsppp3;
    double xpspp3;
    double xppsp3;
    double xppps3;
    double xpppp3;
    double xsppp4;
    double xppsp4;
    double xpppp4;
    double xpppp5;
    double ysssp1;
    double yssps1;
    double yspss1;
    double ypsss1;
    double ysspp1;
    double yspsp1;
    double ypssp1;
    double yspps1;
    double ypsps1;
    double yppss1;
    double ysppp1;
    double ypspp1;
    double yppsp1;
    double yppps1;
    double ypppp1;
    double ysssp2;
    double yspss2;
    double ysspp2;
    double yspsp2;
    double ypssp2;
    double yspps2;
    double ypsps2;
    double yppss2;
    double ysppp2;
    double ypspp2;
    double yppsp2;
    double yppps2;
    double ypppp2;
    double ysspp3;
    double yspsp3;
    double yppss3;
    double ysppp3;
    double ypspp3;
    double yppsp3;
    double yppps3;
    double ypppp3;   
    double ysppp4;
    double yppsp4;
    double ypppp4;
    double ypppp5;
    double zsssp1;
    double zssps1;
    double zspss1;
    double zpsss1;
    double zsspp1;
    double zspsp1;
    double zpssp1;
    double zspps1;
    double zpsps1;
    double zppss1;
    double zsppp1;
    double zpspp1;
    double zppsp1;
    double zppps1;
    double zpppp1;
    double zsssp2;
    double zspss2;
    double zsspp2;
    double zspsp2;
    double zpssp2;
    double zspps2;
    double zpsps2;
    double zppss2;
    double zsppp2;
    double zpspp2;
    double zppsp2;
    double zppps2;
    double zpppp2;
    double zsspp3;
    double zspsp3;
    double zppss3;
    double zsppp3;
    double zpspp3;
    double zppsp3;
    double zppps3;
    double zpppp3;
    double zsppp4;
    double zppsp4;
    double x12;
    double y12;
    double z12;
    double x34;
    double y34;
    double z34;
    
    ftable_dim1 = mgrid - 0 + 1;
    ftable_offset = 0 + ftable_dim1 * 0;
    ftable -= ftable_offset;
    x12 = x1 - x2;
    y12 = y1 - y2;
    z12 = z1 - z2;
    x34 = x3 - x4;
    y34 = y3 - y4;
    z34 = z3 - z4;
    
    for (ij = 0; ij < nij; ++ij)
    {
        i = prim1[ij];
        j = prim2[ij];
        exp1 = alpha1[i - 1];
        exp2 = alpha2[j - 1];
        pval = exp1 + exp2;
        p[ij] = pval;
        pval = exp1 / pval;
        px[ij] = pval * x12 + x2;
        py[ij] = pval * y12 + y2;
        pz[ij] = pval * z12 + z2;
        scalep[ij] = norm1[i - 1] * norm2[j - 1] * rho12[ij];
    }

    for (kl = 0; kl < nkl; ++kl)
    {
        k = prim3[kl];
        l = prim4[kl];
        exp3 = alpha3[k - 1];
        exp4 = alpha4[l - 1];
        qval = exp3 + exp4;
        q[kl] = qval;
        qval = exp3 / qval;
        qx[kl] = qval * x34 + x4;
        qy[kl] = qval * y34 + y4;
        qz[kl] = qval * z34 + z4;
        scaleq[kl] = norm3[k - 1] * norm4[l - 1] * rho34[kl];
    }

    m = 0;
    for (ij = 0; ij < nij; ++ij)
    {
        pval = p[ij];
        pxval = px[ij];
        pyval = py[ij];
        pzval = pz[ij];
        pscale = scalep[ij];
        u3 = .5 / pval;
        xspss1 = pxval - x2;
        xpsss1 = pxval - x1;
        yspss1 = pyval - y2;
        ypsss1 = pyval - y1;
        zspss1 = pzval - z2;
        zpsss1 = pzval - z1;
        for (kl = 0; kl < nkl; ++kl)
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
            f4 = scale * boys4(t);
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
                f4 = (((((ftable[tgrid * ftable_dim1 + 10] * delta6 +
                          ftable[tgrid * ftable_dim1 + 9]) * delta5 +
                         ftable[tgrid * ftable_dim1 + 8]) * delta4 +
                        ftable[tgrid * ftable_dim1 + 7]) * delta3 +
                       ftable[tgrid * ftable_dim1 + 6]) * delta2 +
                      ftable[tgrid * ftable_dim1 + 5]) * delta1 +
                    ftable[tgrid * ftable_dim1 + 4];
                f0 = scale * f0;
                f1 = scale * f1;
                f2 = scale * f2;
                f3 = scale * f3;
                f4 = scale * f4;
            }
            else
            {
                tinv = 1. / t;
                t2inv = tinv * .5;
                f0 = scale * .5 * sqrt (tinv * 3.141592653589793);
                f1 = t2inv * f0;
                f2 = t2inv * 3. * f1;
                f3 = t2inv * 5. * f2;
                f4 = t2inv * 7. * f3;
            }
#endif
            u0 = pval * pqpinv;
            u1 = -qval * pqpinv;
            u2 = pqpinv * .5;
            u4 = .5 / qval;
            u5 = u2 + pqpinv;
            u6 = u5 + pqpinv;
            u7 = -u0 * u4;
            u8 = u1 * u3;

/*             ...the X-terms (with exception of XSPSS1 and XPSSS1, */
/*                which can be evaluated in the P-loop). */
            xsssp1 = qxval - x4;
            xssps1 = qxval - x3;
            xsssp2 = pqx * u0;
            xspss2 = pqx * u1;
            a = xsssp1 + xssps1;
            b = xspss1 + xpsss1;
            c = a * b;
            xspsp1 = xspss1 * xsssp1;
            xpssp1 = xpsss1 * xsssp1;
            xspps1 = xspss1 * xssps1;
            xpsps1 = xpsss1 * xssps1;
            xsspp1 = xsssp1 * xssps1 + u4;
            xppss1 = xspss1 * xpsss1 + u3;
            d = xsssp1 * xspss2 + u2;
            e = xssps1 * xspss2 + u2;
            xspsp2 = xspss1 * xsssp2 + d;
            xpssp2 = xpsss1 * xsssp2 + d;
            xspps2 = xspss1 * xsssp2 + e;
            xpsps2 = xpsss1 * xsssp2 + e;
            xsspp2 = a * xsssp2 + u7;
            xppss2 = b * xspss2 + u8;
            xspsp3 = xsssp2 * xspss2;
            xsspp3 = xsssp2 * xsssp2;
            xppss3 = xspss2 * xspss2;
            xsppp1 = xsspp1 * xspss1;
            xpspp1 = xsspp1 * xpsss1;
            xppsp1 = xppss1 * xsssp1;
            xppps1 = xppss1 * xssps1;
            d = xspss2 * xsspp1 + a * u2;
            e = xsssp2 * xppss1 + b * u2;
            xsppp2 = xspss1 * xsspp2 + d;
            xpspp2 = xpsss1 * xsspp2 + d;
            xppsp2 = xsssp1 * xppss2 + e;
            xppps2 = xssps1 * xppss2 + e;
            d = a * xspsp3 + xsssp2 * u5;
            e = b * xspsp3 + xspss2 * u5;
            xsppp3 = xspss1 * xsspp3 + d;
            xpspp3 = xpsss1 * xsspp3 + d;
            xppsp3 = xsssp1 * xppss3 + e;
            xppps3 = xssps1 * xppss3 + e;
            xsppp4 = xsspp3 * xspss2;
            xppsp4 = xppss3 * xsssp2;
            d = b * xsssp2 + a * xspss2 + u2;
            xpppp1 = xsspp1 * xppss1;
            xpppp2 = xsspp2 * xppss1 + xppss2 * xsspp1 + c * u2;
            xpppp3 = xppss1 * xsspp3 + xsspp1 * xppss3 + xspsp3 * c + u5 *
                d;
            xpppp4 = xspsp3 * (d + u6);
            xpppp5 = xsspp3 * xppss3;

/*             ...the Y-terms (with exception of YSPSS1 and YPSSS1, */
/*                which can be evaluated in the P-loop). */
            ysssp1 = qyval - y4;
            yssps1 = qyval - y3;
            ysssp2 = pqy * u0;
            yspss2 = pqy * u1;
            a = ysssp1 + yssps1;
            b = yspss1 + ypsss1;
            c = a * b;
            yspsp1 = yspss1 * ysssp1;
            ypssp1 = ypsss1 * ysssp1;
            yspps1 = yspss1 * yssps1;
            ypsps1 = ypsss1 * yssps1;
            ysspp1 = ysssp1 * yssps1 + u4;
            yppss1 = yspss1 * ypsss1 + u3;
            d = ysssp1 * yspss2 + u2;
            e = yssps1 * yspss2 + u2;
            yspsp2 = yspss1 * ysssp2 + d;
            ypssp2 = ypsss1 * ysssp2 + d;
            yspps2 = yspss1 * ysssp2 + e;
            ypsps2 = ypsss1 * ysssp2 + e;
            ysspp2 = a * ysssp2 + u7;
            yppss2 = b * yspss2 + u8;
            yspsp3 = ysssp2 * yspss2;
            ysspp3 = ysssp2 * ysssp2;
            yppss3 = yspss2 * yspss2;
            ysppp1 = ysspp1 * yspss1;
            ypspp1 = ysspp1 * ypsss1;
            yppsp1 = yppss1 * ysssp1;
            yppps1 = yppss1 * yssps1;
            d = yspss2 * ysspp1 + a * u2;
            e = ysssp2 * yppss1 + b * u2;
            ysppp2 = yspss1 * ysspp2 + d;
            ypspp2 = ypsss1 * ysspp2 + d;
            yppsp2 = ysssp1 * yppss2 + e;
            yppps2 = yssps1 * yppss2 + e;
            d = a * yspsp3 + ysssp2 * u5;
            e = b * yspsp3 + yspss2 * u5;
            ysppp3 = yspss1 * ysspp3 + d;
            ypspp3 = ypsss1 * ysspp3 + d;
            yppsp3 = ysssp1 * yppss3 + e;
            yppps3 = yssps1 * yppss3 + e;
            ysppp4 = ysspp3 * yspss2;
            yppsp4 = yppss3 * ysssp2;
            d = b * ysssp2 + a * yspss2 + u2;
            ypppp1 = ysspp1 * yppss1;
            ypppp2 = ysspp2 * yppss1 + yppss2 * ysspp1 + c * u2;
            ypppp3 = yppss1 * ysspp3 + ysspp1 * yppss3 + yspsp3 * c + u5 *
                d;
            ypppp4 = yspsp3 * (d + u6);
            ypppp5 = ysspp3 * yppss3;

/*             ...the Z-terms (with exception of ZSPSS1 and ZPSSS1, */
/*                which can be evaluated in the P-loop). */
            zsssp1 = qzval - z4;
            zssps1 = qzval - z3;
            zsssp2 = pqz * u0;
            zspss2 = pqz * u1;
            a = zsssp1 + zssps1;
            b = zspss1 + zpsss1;
            c = a * b;
            zspsp1 = zspss1 * zsssp1;
            zpssp1 = zpsss1 * zsssp1;
            zspps1 = zspss1 * zssps1;
            zpsps1 = zpsss1 * zssps1;
            zsspp1 = zsssp1 * zssps1 + u4;
            zppss1 = zspss1 * zpsss1 + u3;
            d = zsssp1 * zspss2 + u2;
            e = zssps1 * zspss2 + u2;
            zspsp2 = zspss1 * zsssp2 + d;
            zpssp2 = zpsss1 * zsssp2 + d;
            zspps2 = zspss1 * zsssp2 + e;
            zpsps2 = zpsss1 * zsssp2 + e;
            zsspp2 = a * zsssp2 + u7;
            zppss2 = b * zspss2 + u8;
            zspsp3 = zsssp2 * zspss2;
            zsspp3 = zsssp2 * zsssp2;
            zppss3 = zspss2 * zspss2;
            zsppp1 = zsspp1 * zspss1;
            zpspp1 = zsspp1 * zpsss1;
            zppsp1 = zppss1 * zsssp1;
            zppps1 = zppss1 * zssps1;
            d = zspss2 * zsspp1 + a * u2;
            e = zsssp2 * zppss1 + b * u2;
            zsppp2 = zspss1 * zsspp2 + d;
            zpspp2 = zpsss1 * zsspp2 + d;
            zppsp2 = zsssp1 * zppss2 + e;
            zppps2 = zssps1 * zppss2 + e;
            d = a * zspsp3 + zsssp2 * u5;
            e = b * zspsp3 + zspss2 * u5;
            zsppp3 = zspss1 * zsspp3 + d;
            zpspp3 = zpsss1 * zsspp3 + d;
            zppsp3 = zsssp1 * zppss3 + e;
            zppps3 = zssps1 * zppss3 + e;
            zsppp4 = zsspp3 * zspss2;
            zppsp4 = zppss3 * zsssp2;
            d = b * zsssp2 + a * zspss2 + u2;
            zpppp1 = zsspp1 * zppss1;
            zpppp2 = zsspp2 * zppss1 + zppss2 * zsspp1 + c * u2;
            zpppp3 = zppss1 * zsspp3 + zsspp1 * zppss3 + zspsp3 * c + u5 *
                d;
            zpppp4 = zspsp3 * (d + u6);
            zpppp5 = zsspp3 * zppss3;

/*             ...assemble the 4-center (AB|CD) type integrals. */
            gxxxx = xpppp1 * f0 + xpppp2 * f1 + xpppp3 * f2 + xpppp4 * f3 +
                xpppp5 * f4;
            gyyyy = ypppp1 * f0 + ypppp2 * f1 + ypppp3 * f2 + ypppp4 * f3 +
                ypppp5 * f4;
            gzzzz = zpppp1 * f0 + zpppp2 * f1 + zpppp3 * f2 + zpppp4 * f3 +
                zpppp5 * f4;
            a = xppps1 * f0 + xppps2 * f1 + xppps3 * f2 + xppsp4 * f3;
            b = xppps1 * f1 + xppps2 * f2 + xppps3 * f3 + xppsp4 * f4;
            c = xppsp1 * f0 + xppsp2 * f1 + xppsp3 * f2 + xppsp4 * f3;
            d = xppsp1 * f1 + xppsp2 * f2 + xppsp3 * f3 + xppsp4 * f4;
            e = xpspp1 * f0 + xpspp2 * f1 + xpspp3 * f2 + xsppp4 * f3;
            f = xpspp1 * f1 + xpspp2 * f2 + xpspp3 * f3 + xsppp4 * f4;
            g = xsppp1 * f0 + xsppp2 * f1 + xsppp3 * f2 + xsppp4 * f3;
            h = xsppp1 * f1 + xsppp2 * f2 + xsppp3 * f3 + xsppp4 * f4;
            gxxxy = a * ysssp1 + b * ysssp2;
            gxxyx = c * yssps1 + d * ysssp2;
            gxyxx = e * yspss1 + f * yspss2;
            gyxxx = g * ypsss1 + h * yspss2;
            gxxxz = a * zsssp1 + b * zsssp2;
            gxxzx = c * zssps1 + d * zsssp2;
            gxzxx = e * zspss1 + f * zspss2;
            gzxxx = g * zpsss1 + h * zspss2;
            a = yppps1 * f0 + yppps2 * f1 + yppps3 * f2 + yppsp4 * f3;
            b = yppps1 * f1 + yppps2 * f2 + yppps3 * f3 + yppsp4 * f4;
            c = yppsp1 * f0 + yppsp2 * f1 + yppsp3 * f2 + yppsp4 * f3;
            d = yppsp1 * f1 + yppsp2 * f2 + yppsp3 * f3 + yppsp4 * f4;
            e = ypspp1 * f0 + ypspp2 * f1 + ypspp3 * f2 + ysppp4 * f3;
            f = ypspp1 * f1 + ypspp2 * f2 + ypspp3 * f3 + ysppp4 * f4;
            g = ysppp1 * f0 + ysppp2 * f1 + ysppp3 * f2 + ysppp4 * f3;
            h = ysppp1 * f1 + ysppp2 * f2 + ysppp3 * f3 + ysppp4 * f4;
            gyyyx = a * xsssp1 + b * xsssp2;
            gyyxy = c * xssps1 + d * xsssp2;
            gyxyy = e * xspss1 + f * xspss2;
            gxyyy = g * xpsss1 + h * xspss2;
            gyyyz = a * zsssp1 + b * zsssp2;
            gyyzy = c * zssps1 + d * zsssp2;
            gyzyy = e * zspss1 + f * zspss2;
            gzyyy = g * zpsss1 + h * zspss2;
            a = zppps1 * f0 + zppps2 * f1 + zppps3 * f2 + zppsp4 * f3;
            b = zppps1 * f1 + zppps2 * f2 + zppps3 * f3 + zppsp4 * f4;
            c = zppsp1 * f0 + zppsp2 * f1 + zppsp3 * f2 + zppsp4 * f3;
            d = zppsp1 * f1 + zppsp2 * f2 + zppsp3 * f3 + zppsp4 * f4;
            e = zpspp1 * f0 + zpspp2 * f1 + zpspp3 * f2 + zsppp4 * f3;
            f = zpspp1 * f1 + zpspp2 * f2 + zpspp3 * f3 + zsppp4 * f4;
            g = zsppp1 * f0 + zsppp2 * f1 + zsppp3 * f2 + zsppp4 * f3;
            h = zsppp1 * f1 + zsppp2 * f2 + zsppp3 * f3 + zsppp4 * f4;
            gzzzx = a * xsssp1 + b * xsssp2;
            gzzxz = c * xssps1 + d * xsssp2;
            gzxzz = e * xspss1 + f * xspss2;
            gxzzz = g * xpsss1 + h * xspss2;
            gzzzy = a * ysssp1 + b * ysssp2;
            gzzyz = c * yssps1 + d * ysssp2;
            gzyzz = e * yspss1 + f * yspss2;
            gyzzz = g * ypsss1 + h * yspss2;
            a = xppss1 * f0 + xppss2 * f1 + xppss3 * f2;
            b = xppss1 * f1 + xppss2 * f2 + xppss3 * f3;
            c = xppss1 * f2 + xppss2 * f3 + xppss3 * f4;
            d = yppss1 * f0 + yppss2 * f1 + yppss3 * f2;
            e = yppss1 * f1 + yppss2 * f2 + yppss3 * f3;
            f = yppss1 * f2 + yppss2 * f3 + yppss3 * f4;
            g = zppss1 * f0 + zppss2 * f1 + zppss3 * f2;
            h = zppss1 * f1 + zppss2 * f2 + zppss3 * f3;
            r = zppss1 * f2 + zppss2 * f3 + zppss3 * f4;
            gxxyy = a * ysspp1 + b * ysspp2 + c * ysspp3;
            gxxzz = a * zsspp1 + b * zsspp2 + c * zsspp3;
            gyyxx = d * xsspp1 + e * xsspp2 + f * xsspp3;
            gyyzz = d * zsspp1 + e * zsspp2 + f * zsspp3;
            gzzxx = g * xsspp1 + h * xsspp2 + r * xsspp3;
            gzzyy = g * ysspp1 + h * ysspp2 + r * ysspp3;
            gxxyz = (a * yssps1 + b * ysssp2) * zsssp1 + (b * yssps1 + c *
                                                          ysssp2) * zsssp2;
            gxxzy = (a * zssps1 + b * zsssp2) * ysssp1 + (b * zssps1 + c *
                                                          zsssp2) * ysssp2;
            gyyxz = (d * xssps1 + e * xsssp2) * zsssp1 + (e * xssps1 + f *
                                                            xsssp2) * zsssp2;
            gyyzx = (d * zssps1 + e * zsssp2) * xsssp1 + (e * zssps1 + f *
                                                            zsssp2) * xsssp2;
            gzzxy = (g * xssps1 + h * xsssp2) * ysssp1 + (h * xssps1 +
                                                            r * xsssp2) *
                ysssp2;
            gzzyx =
                (g * yssps1 + h * ysssp2) * xsssp1 + (h * yssps1 +
                                                        r * ysssp2) *
                xsssp2;
            aa = xspsp3 * f2;
            bb = xspsp3 * f3;
            cc = xspsp3 * f4;
            dd = yspsp3 * f2;
            ee = yspsp3 * f3;
            ff = yspsp3 * f4;
            gg = zspsp3 * f2;
            hh = zspsp3 * f3;
            rr = zspsp3 * f4;
            a = xpsps1 * f0 + xpsps2 * f1 + aa;
            b = xpsps1 * f1 + xpsps2 * f2 + bb;
            c = xpsps1 * f2 + xpsps2 * f3 + cc;
            d = ypsps1 * f0 + ypsps2 * f1 + dd;
            e = ypsps1 * f1 + ypsps2 * f2 + ee;
            f = ypsps1 * f2 + ypsps2 * f3 + ff;
            g = zpsps1 * f0 + zpsps2 * f1 + gg;
            h = zpsps1 * f1 + zpsps2 * f2 + hh;
            r = zpsps1 * f2 + zpsps2 * f3 + rr;
            gxyxy = a * yspsp1 + b * yspsp2 + c * yspsp3;
            gxzxz = a * zspsp1 + b * zspsp2 + c * zspsp3;
            gyxyx = d * xspsp1 + e * xspsp2 + f * xspsp3;
            gyzyz = d * zspsp1 + e * zspsp2 + f * zspsp3;
            gzxzx = g * xspsp1 + h * xspsp2 + r * xspsp3;
            gzyzy = g * yspsp1 + h * yspsp2 + r * yspsp3;
            gxyxz = (a * yspss1 + b * yspss2) * zsssp1 + (b * yspss1 + c *
                                                          yspss2) * zsssp2;
            gxzxy = (a * zspss1 + b * zspss2) * ysssp1 + (b * zspss1 + c *
                                                          zspss2) * ysssp2;
            gyxyz = (d * xspss1 + e * xspss2) * zsssp1 + (e * xspss1 + f *
                                                            xspss2) * zsssp2;
            gyzyx = (d * zspss1 + e * zspss2) * xsssp1 + (e * zspss1 + f *
                                                            zspss2) * xsssp2;
            gzxzy = (g * xspss1 + h * xspss2) * ysssp1 + (h * xspss1 +
                                                            r * xspss2) *
                ysssp2;
            gzyzx =
                (g * yspss1 + h * yspss2) * xsssp1 + (h * yspss1 +
                                                        r * yspss2) *
                xsssp2;
            a = xpssp1 * f0 + xpssp2 * f1 + aa;
            b = xpssp1 * f1 + xpssp2 * f2 + bb;
            c = xpssp1 * f2 + xpssp2 * f3 + cc;
            d = ypssp1 * f0 + ypssp2 * f1 + dd;
            e = ypssp1 * f1 + ypssp2 * f2 + ee;
            f = ypssp1 * f2 + ypssp2 * f3 + ff;
            g = zpssp1 * f0 + zpssp2 * f1 + gg;
            h = zpssp1 * f1 + zpssp2 * f2 + hh;
            r = zpssp1 * f2 + zpssp2 * f3 + rr;
            gxyyx = a * yspps1 + b * yspps2 + c * yspsp3;
            gxzzx = a * zspps1 + b * zspps2 + c * zspsp3;
            gyxxy = d * xspps1 + e * xspps2 + f * xspsp3;
            gyzzy = d * zspps1 + e * zspps2 + f * zspsp3;
            gzxxz = g * xspps1 + h * xspps2 + r * xspsp3;
            gzyyz = g * yspps1 + h * yspps2 + r * yspsp3;
            gxyzx = (a * yspss1 + b * yspss2) * zssps1 + (b * yspss1 + c *
                                                          yspss2) * zsssp2;
            gxzyx = (a * zspss1 + b * zspss2) * yssps1 + (b * zspss1 + c *
                                                          zspss2) * ysssp2;
            gyxzy = (d * xspss1 + e * xspss2) * zssps1 + (e * xspss1 + f *
                                                            xspss2) * zsssp2;
            gyzxy = (d * zspss1 + e * zspss2) * xssps1 + (e * zspss1 + f *
                                                            zspss2) * xsssp2;
            gzxyz = (g * xspss1 + h * xspss2) * yssps1 + (h * xspss1 +
                                                            r * xspss2) *
                ysssp2;
            gzyxz =
                (g * yspss1 + h * yspss2) * xssps1 + (h * yspss1 +
                                                        r * yspss2) *
                xsssp2;
            a = xspps1 * f0 + xspps2 * f1 + aa;
            b = xspps1 * f1 + xspps2 * f2 + bb;
            c = xspps1 * f2 + xspps2 * f3 + cc;
            d = yspps1 * f0 + yspps2 * f1 + dd;
            e = yspps1 * f1 + yspps2 * f2 + ee;
            f = yspps1 * f2 + yspps2 * f3 + ff;
            g = zspps1 * f0 + zspps2 * f1 + gg;
            h = zspps1 * f1 + zspps2 * f2 + hh;
            r = zspps1 * f2 + zspps2 * f3 + rr;
            gyxxz = (a * ypsss1 + b * yspss2) * zsssp1 + (b * ypsss1 + c *
                                                          yspss2) * zsssp2;
            gzxxy = (a * zpsss1 + b * zspss2) * ysssp1 + (b * zpsss1 + c *
                                                          zspss2) * ysssp2;
            gxyyz = (d * xpsss1 + e * xspss2) * zsssp1 + (e * xpsss1 + f *
                                                            xspss2) * zsssp2;
            gzyyx = (d * zpsss1 + e * zspss2) * xsssp1 + (e * zpsss1 + f *
                                                            zspss2) * xsssp2;
            gxzzy = (g * xpsss1 + h * xspss2) * ysssp1 + (h * xpsss1 +
                                                            r * xspss2) *
                ysssp2;
            gyzzx =
                (g * ypsss1 + h * yspss2) * xsssp1 + (h * ypsss1 +
                                                        r * yspss2) *
                xsssp2;
            a = xspsp1 * f0 + xspsp2 * f1 + aa;
            b = xspsp1 * f1 + xspsp2 * f2 + bb;
            c = xspsp1 * f2 + xspsp2 * f3 + cc;
            d = yspsp1 * f0 + yspsp2 * f1 + dd;
            e = yspsp1 * f1 + yspsp2 * f2 + ee;
            f = yspsp1 * f2 + yspsp2 * f3 + ff;
            g = zspsp1 * f0 + zspsp2 * f1 + gg;
            h = zspsp1 * f1 + zspsp2 * f2 + hh;
            r = zspsp1 * f2 + zspsp2 * f3 + rr;
            gyxzx = (a * ypsss1 + b * yspss2) * zssps1 + (b * ypsss1 + c *
                                                          yspss2) * zsssp2;
            gzxyx = (a * zpsss1 + b * zspss2) * yssps1 + (b * zpsss1 + c *
                                                          zspss2) * ysssp2;
            gxyzy = (d * xpsss1 + e * xspss2) * zssps1 + (e * xpsss1 + f *
                                                            xspss2) * zsssp2;
            gzyxy = (d * zpsss1 + e * zspss2) * xssps1 + (e * zpsss1 + f *
                                                            zspss2) * xsssp2;
            gxzyz = (g * xpsss1 + h * xspss2) * yssps1 + (h * xpsss1 +
                                                            r * xspss2) *
                ysssp2;
            gyzxz =
                (g * ypsss1 + h * yspss2) * xssps1 + (h * ypsss1 +
                                                        r * yspss2) *
                xsssp2;
            a = xsspp1 * f0 + xsspp2 * f1 + xsspp3 * f2;
            b = xsspp1 * f1 + xsspp2 * f2 + xsspp3 * f3;
            c = xsspp1 * f2 + xsspp2 * f3 + xsspp3 * f4;
            d = ysspp1 * f0 + ysspp2 * f1 + ysspp3 * f2;
            e = ysspp1 * f1 + ysspp2 * f2 + ysspp3 * f3;
            f = ysspp1 * f2 + ysspp2 * f3 + ysspp3 * f4;
            g = zsspp1 * f0 + zsspp2 * f1 + zsspp3 * f2;
            h = zsspp1 * f1 + zsspp2 * f2 + zsspp3 * f3;
            r = zsspp1 * f2 + zsspp2 * f3 + zsspp3 * f4;
            gyzxx = (a * ypsss1 + b * yspss2) * zspss1 + (b * ypsss1 + c *
                                                          yspss2) * zspss2;
            gzyxx = (a * zpsss1 + b * zspss2) * yspss1 + (b * zpsss1 + c *
                                                          zspss2) * yspss2;
            gxzyy = (d * xpsss1 + e * xspss2) * zspss1 + (e * xpsss1 + f *
                                                            xspss2) * zspss2;
            gzxyy = (d * zpsss1 + e * zspss2) * xspss1 + (e * zpsss1 + f *
                                                            zspss2) * xspss2;
            gxyzz = (g * xpsss1 + h * xspss2) * yspss1 + (h * xpsss1 +
                                                            r * xspss2) *yspss2;
            gyxzz =
                (g * ypsss1 + h * yspss2) * xspss1 + (h * ypsss1 +
                                                        r * yspss2) * xspss2;
            batch[m + 0] = gxxxx;
            batch[m + 1] = gyxxx;
            batch[m + 2] = gzxxx;
            batch[m + 3] = gxyxx;
            batch[m + 4] = gyyxx;
            batch[m + 5] = gzyxx;
            batch[m + 6] = gxzxx;
            batch[m + 7] = gyzxx;
            batch[m + 8] = gzzxx;
            batch[m + 9] = gxxyx;
            batch[m + 10] = gyxyx;
            batch[m + 11] = gzxyx;
            batch[m + 12] = gxyyx;
            batch[m + 13] = gyyyx;
            batch[m + 14] = gzyyx;
            batch[m + 15] = gxzyx;
            batch[m + 16] = gyzyx;
            batch[m + 17] = gzzyx;
            batch[m + 18] = gxxzx;
            batch[m + 19] = gyxzx;
            batch[m + 20] = gzxzx;
            batch[m + 21] = gxyzx;
            batch[m + 22] = gyyzx;
            batch[m + 23] = gzyzx;
            batch[m + 24] = gxzzx;
            batch[m + 25] = gyzzx;
            batch[m + 26] = gzzzx;
            batch[m + 27] = gxxxy;
            batch[m + 28] = gyxxy;
            batch[m + 29] = gzxxy;
            batch[m + 30] = gxyxy;
            batch[m + 31] = gyyxy;
            batch[m + 32] = gzyxy;
            batch[m + 33] = gxzxy;
            batch[m + 34] = gyzxy;
            batch[m + 35] = gzzxy;
            batch[m + 36] = gxxyy;
            batch[m + 37] = gyxyy;
            batch[m + 38] = gzxyy;
            batch[m + 39] = gxyyy;
            batch[m + 40] = gyyyy;
            batch[m + 41] = gzyyy;
            batch[m + 42] = gxzyy;
            batch[m + 43] = gyzyy;
            batch[m + 44] = gzzyy;
            batch[m + 45] = gxxzy;
            batch[m + 46] = gyxzy;
            batch[m + 47] = gzxzy;
            batch[m + 48] = gxyzy;
            batch[m + 49] = gyyzy;
            batch[m + 50] = gzyzy;
            batch[m + 51] = gxzzy;
            batch[m + 52] = gyzzy;
            batch[m + 53] = gzzzy;
            batch[m + 54] = gxxxz;
            batch[m + 55] = gyxxz;
            batch[m + 56] = gzxxz;
            batch[m + 57] = gxyxz;
            batch[m + 58] = gyyxz;
            batch[m + 59] = gzyxz;
            batch[m + 60] = gxzxz;
            batch[m + 61] = gyzxz;
            batch[m + 62] = gzzxz;
            batch[m + 63] = gxxyz;
            batch[m + 64] = gyxyz;
            batch[m + 65] = gzxyz;
            batch[m + 66] = gxyyz;
            batch[m + 67] = gyyyz;
            batch[m + 68] = gzyyz;
            batch[m + 69] = gxzyz;
            batch[m + 70] = gyzyz;
            batch[m + 71] = gzzyz;
            batch[m + 72] = gxxzz;
            batch[m + 73] = gyxzz;
            batch[m + 74] = gzxzz;
            batch[m + 75] = gxyzz;
            batch[m + 76] = gyyzz;
            batch[m + 77] = gzyzz;
            batch[m + 78] = gxzzz;
            batch[m + 79] = gyzzz;
            batch[m + 80] = gzzzz;
            m += 81;
        }
    }

    return 0;
}