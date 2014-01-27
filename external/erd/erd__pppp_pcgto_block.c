#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <yepPredefines.h>
#include "boys.h"
//#define ERD_TABLE_FREE_BOYS_FUNCTIONS


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
void erd__pppp_pcgto_block(int nij, int nkl,
                           double x1, double y1, double z1,
                           double x2, double y2, double z2,
                           double x3, double y3, double z3,
                           double x4, double y4, double z4,
                           double *YEP_RESTRICT alpha1, double *YEP_RESTRICT alpha2,
                           double *YEP_RESTRICT alpha3, double *YEP_RESTRICT alpha4,
                           double *YEP_RESTRICT cc1, double *YEP_RESTRICT cc2,
                           double *YEP_RESTRICT cc3, double *YEP_RESTRICT cc4,
                           double *YEP_RESTRICT ftable, int mgrid,
                           double tmax, double tstep, double tvstep,
                           int *YEP_RESTRICT prim1, int *YEP_RESTRICT prim2,
                           int *YEP_RESTRICT prim3, int *YEP_RESTRICT prim4,
                           double *YEP_RESTRICT norm1, double *YEP_RESTRICT norm2,
                           double *YEP_RESTRICT norm3, double *YEP_RESTRICT norm4,
                           double *YEP_RESTRICT rho12, double *YEP_RESTRICT rho34,
                           double *YEP_RESTRICT p, double *YEP_RESTRICT px,
                           double *YEP_RESTRICT py, double *YEP_RESTRICT pz, double *YEP_RESTRICT scalep,
                           double *YEP_RESTRICT q, double *YEP_RESTRICT qx,
                           double *YEP_RESTRICT qy, double *YEP_RESTRICT qz,
                           double *YEP_RESTRICT scaleq, double *YEP_RESTRICT cbatch)
{
    const int ftable_dim1 = mgrid + 1;
    const double x12 = x1 - x2;
    const double y12 = y1 - y2;
    const double z12 = z1 - z2;
    const double x34 = x3 - x4;
    const double y34 = y3 - y4;
    const double z34 = z3 - z4;
    
    for (int ij = 0; ij < nij; ij += 1) {
        const int i = prim1[ij];
        const int j = prim2[ij];
        const double exp1 = alpha1[i];
        const double exp2 = alpha2[j];
        double pval = exp1 + exp2;
        p[ij] = pval;
        pval = exp1 / pval;
        px[ij] = pval * x12 + x2;
        py[ij] = pval * y12 + y2;
        pz[ij] = pval * z12 + z2;
        scalep[ij] = cc1[i] * cc2[j] * norm1[i] * norm2[j] * rho12[ij];
    }

    for (int kl = 0; kl < nkl; kl += 1) {
        const int k = prim3[kl];
        const int l = prim4[kl];
        const double exp3 = alpha3[k];
        const double exp4 = alpha4[l];
        double qval = exp3 + exp4;
        q[kl] = qval;
        qval = exp3 / qval;
        qx[kl] = qval * x34 + x4;
        qy[kl] = qval * y34 + y4;
        qz[kl] = qval * z34 + z4;
        scaleq[kl] = cc3[k] * cc4[l] * norm3[k] * norm4[l] * rho34[kl];
    }

    for (int ij = 0; ij < nij; ij += 1) {
        const double pval = p[ij];
        const double pxval = px[ij];
        const double pyval = py[ij];
        const double pzval = pz[ij];
        const double pscale = scalep[ij];
        const double u3 = .5 / pval;
        const double xspss1 = pxval - x2;
        const double xpsss1 = pxval - x1;
        const double yspss1 = pyval - y2;
        const double ypsss1 = pyval - y1;
        const double zspss1 = pzval - z2;
        const double zpsss1 = pzval - z1;
        for (int kl = 0; kl < nkl; kl += 1) {
            const double qval = q[kl];
            const double qxval = qx[kl];
            const double qyval = qy[kl];
            const double qzval = qz[kl];
            const double pqmult = pval * qval;
            const double pqplus = pval + qval;
            const double pqpinv = 1. / pqplus;
            const double pqx = pxval - qxval;
            const double pqy = pyval - qyval;
            const double pqz = pzval - qzval;
            const double t = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult * pqpinv;
            const double scale = pscale * scaleq[kl] / (pqmult * __builtin_sqrt(pqplus));
#ifdef ERD_TABLE_FREE_BOYS_FUNCTIONS
            const double f0 = scale * boys0(t);
            const double f1 = scale * boys1(t);
            const double f2 = scale * boys2(t);
            const double f3 = scale * boys3(t);
            const double f4 = scale * boys4(t);
#else
            double f0, f1, f2, f3, f4;
            if (t <= tmax) {
                const int tgrid = __builtin_lround(t * tvstep);
                const double delta1 = tgrid * tstep - t;
                const double delta2 = delta1 * 0x1.0000000000000p-1;
                const double delta3 = delta1 * 0x1.5555555555555p-2;
                const double delta4 = delta1 * 0x1.0000000000000p-2;
                const double delta5 = delta1 * 0x1.999999999999Ap-3;
                const double delta6 = delta1 * 0x1.5555555555555p-3;
                f0 = (((((boys_table[tgrid][6] * delta6 +
                    boys_table[tgrid][5]) * delta5 +
                        boys_table[tgrid][4]) * delta4 +
                            boys_table[tgrid][3]) * delta3 +
                                boys_table[tgrid][2]) * delta2 +
                                    boys_table[tgrid][1]) * delta1 +
                                        boys_table[tgrid][0];
                f1 = (((((boys_table[tgrid][7] * delta6 +
                    boys_table[tgrid][6]) * delta5 +
                        boys_table[tgrid][5]) * delta4 +
                            boys_table[tgrid][4]) * delta3 +
                                boys_table[tgrid][3]) * delta2 +
                                    boys_table[tgrid][2]) * delta1 +
                                        boys_table[tgrid][1];
                f2 = (((((boys_table[tgrid][8] * delta6 +
                    boys_table[tgrid][7]) * delta5 +
                        boys_table[tgrid][6]) * delta4 +
                            boys_table[tgrid][5]) * delta3 +
                                boys_table[tgrid][4]) * delta2 +
                                    boys_table[tgrid][3]) * delta1 +
                                        boys_table[tgrid][2];
                f3 = (((((boys_table[tgrid][9] * delta6 +
                    boys_table[tgrid][8]) * delta5 +
                        boys_table[tgrid][7]) * delta4 +
                            boys_table[tgrid][6]) * delta3 +
                                boys_table[tgrid][5]) * delta2 +
                                    boys_table[tgrid][4]) * delta1 +
                                        boys_table[tgrid][3];
                f4 = (((((boys_table[tgrid][10] * delta6 +
                    boys_table[tgrid][9]) * delta5 +
                        boys_table[tgrid][8]) * delta4 +
                            boys_table[tgrid][7]) * delta3 +
                                boys_table[tgrid][6]) * delta2 +
                                    boys_table[tgrid][5]) * delta1 +
                                        boys_table[tgrid][4];
                f0 = scale * f0;
                f1 = scale * f1;
                f2 = scale * f2;
                f3 = scale * f3;
                f4 = scale * f4;
            } else {
                const double tinv = 1. / t;
                const double t2inv = tinv * .5;
                f0 = scale * .5 * __builtin_sqrt(tinv * 3.141592653589793);
                f1 = t2inv * f0;
                f2 = t2inv * 3. * f1;
                f3 = t2inv * 5. * f2;
                f4 = t2inv * 7. * f3;
            }
#endif
            const double u0 = pval * pqpinv;
            const double u1 = -qval * pqpinv;
            const double u2 = pqpinv * .5;
            const double u4 = .5 / qval;
            const double u5 = u2 + pqpinv;
            const double u6 = u5 + pqpinv;
            const double u7 = -u0 * u4;
            const double u8 = u1 * u3;

/*             ...the X-terms (with exception of XSPSS1 and XPSSS1, */
/*                which can be evaluated in the P-loop). */
            const double xsssp1 = qxval - x4;
            const double xssps1 = qxval - x3;
            const double xsssp2 = pqx * u0;
            const double xspss2 = pqx * u1;
            double a = xsssp1 + xssps1;
            double b = xspss1 + xpsss1;
            double c = a * b;
            const double xspsp1 = xspss1 * xsssp1;
            const double xpssp1 = xpsss1 * xsssp1;
            const double xspps1 = xspss1 * xssps1;
            const double xpsps1 = xpsss1 * xssps1;
            const double xsspp1 = xsssp1 * xssps1 + u4;
            const double xppss1 = xspss1 * xpsss1 + u3;
            double d = xsssp1 * xspss2 + u2;
            double e = xssps1 * xspss2 + u2;
            const double xspsp2 = xspss1 * xsssp2 + d;
            const double xpssp2 = xpsss1 * xsssp2 + d;
            const double xspps2 = xspss1 * xsssp2 + e;
            const double xpsps2 = xpsss1 * xsssp2 + e;
            const double xsspp2 = a * xsssp2 + u7;
            const double xppss2 = b * xspss2 + u8;
            const double xspsp3 = xsssp2 * xspss2;
            const double xsspp3 = xsssp2 * xsssp2;
            const double xppss3 = xspss2 * xspss2;
            const double xsppp1 = xsspp1 * xspss1;
            const double xpspp1 = xsspp1 * xpsss1;
            const double xppsp1 = xppss1 * xsssp1;
            const double xppps1 = xppss1 * xssps1;
            d = xspss2 * xsspp1 + a * u2;
            e = xsssp2 * xppss1 + b * u2;
            const double xsppp2 = xspss1 * xsspp2 + d;
            const double xpspp2 = xpsss1 * xsspp2 + d;
            const double xppsp2 = xsssp1 * xppss2 + e;
            const double xppps2 = xssps1 * xppss2 + e;
            d = a * xspsp3 + xsssp2 * u5;
            e = b * xspsp3 + xspss2 * u5;
            const double xsppp3 = xspss1 * xsspp3 + d;
            const double xpspp3 = xpsss1 * xsspp3 + d;
            const double xppsp3 = xsssp1 * xppss3 + e;
            const double xppps3 = xssps1 * xppss3 + e;
            const double xsppp4 = xsspp3 * xspss2;
            const double xppsp4 = xppss3 * xsssp2;
            d = b * xsssp2 + a * xspss2 + u2;
            const double xpppp1 = xsspp1 * xppss1;
            const double xpppp2 = xsspp2 * xppss1 + xppss2 * xsspp1 + c * u2;
            const double xpppp3 = xppss1 * xsspp3 + xsspp1 * xppss3 + xspsp3 * c + u5 * d;
            const double xpppp4 = xspsp3 * (d + u6);
            const double xpppp5 = xsspp3 * xppss3;

/*             ...the Y-terms (with exception of YSPSS1 and YPSSS1, */
/*                which can be evaluated in the P-loop). */
            const double ysssp1 = qyval - y4;
            const double yssps1 = qyval - y3;
            const double ysssp2 = pqy * u0;
            const double yspss2 = pqy * u1;
            a = ysssp1 + yssps1;
            b = yspss1 + ypsss1;
            c = a * b;
            const double yspsp1 = yspss1 * ysssp1;
            const double ypssp1 = ypsss1 * ysssp1;
            const double yspps1 = yspss1 * yssps1;
            const double ypsps1 = ypsss1 * yssps1;
            const double ysspp1 = ysssp1 * yssps1 + u4;
            const double yppss1 = yspss1 * ypsss1 + u3;
            d = ysssp1 * yspss2 + u2;
            e = yssps1 * yspss2 + u2;
            const double yspsp2 = yspss1 * ysssp2 + d;
            const double ypssp2 = ypsss1 * ysssp2 + d;
            const double yspps2 = yspss1 * ysssp2 + e;
            const double ypsps2 = ypsss1 * ysssp2 + e;
            const double ysspp2 = a * ysssp2 + u7;
            const double yppss2 = b * yspss2 + u8;
            const double yspsp3 = ysssp2 * yspss2;
            const double ysspp3 = ysssp2 * ysssp2;
            const double yppss3 = yspss2 * yspss2;
            const double ysppp1 = ysspp1 * yspss1;
            const double ypspp1 = ysspp1 * ypsss1;
            const double yppsp1 = yppss1 * ysssp1;
            const double yppps1 = yppss1 * yssps1;
            d = yspss2 * ysspp1 + a * u2;
            e = ysssp2 * yppss1 + b * u2;
            const double ysppp2 = yspss1 * ysspp2 + d;
            const double ypspp2 = ypsss1 * ysspp2 + d;
            const double yppsp2 = ysssp1 * yppss2 + e;
            const double yppps2 = yssps1 * yppss2 + e;
            d = a * yspsp3 + ysssp2 * u5;
            e = b * yspsp3 + yspss2 * u5;
            const double ysppp3 = yspss1 * ysspp3 + d;
            const double ypspp3 = ypsss1 * ysspp3 + d;
            const double yppsp3 = ysssp1 * yppss3 + e;
            const double yppps3 = yssps1 * yppss3 + e;
            const double ysppp4 = ysspp3 * yspss2;
            const double yppsp4 = yppss3 * ysssp2;
            d = b * ysssp2 + a * yspss2 + u2;
            const double ypppp1 = ysspp1 * yppss1;
            const double ypppp2 = ysspp2 * yppss1 + yppss2 * ysspp1 + c * u2;
            const double ypppp3 = yppss1 * ysspp3 + ysspp1 * yppss3 + yspsp3 * c + u5 * d;
            const double ypppp4 = yspsp3 * (d + u6);
            const double ypppp5 = ysspp3 * yppss3;

/*             ...the Z-terms (with exception of ZSPSS1 and ZPSSS1, */
/*                which can be evaluated in the P-loop). */
            const double zsssp1 = qzval - z4;
            const double zssps1 = qzval - z3;
            const double zsssp2 = pqz * u0;
            const double zspss2 = pqz * u1;
            a = zsssp1 + zssps1;
            b = zspss1 + zpsss1;
            c = a * b;
            const double zspsp1 = zspss1 * zsssp1;
            const double zpssp1 = zpsss1 * zsssp1;
            const double zspps1 = zspss1 * zssps1;
            const double zpsps1 = zpsss1 * zssps1;
            const double zsspp1 = zsssp1 * zssps1 + u4;
            const double zppss1 = zspss1 * zpsss1 + u3;
            d = zsssp1 * zspss2 + u2;
            e = zssps1 * zspss2 + u2;
            const double zspsp2 = zspss1 * zsssp2 + d;
            const double zpssp2 = zpsss1 * zsssp2 + d;
            const double zspps2 = zspss1 * zsssp2 + e;
            const double zpsps2 = zpsss1 * zsssp2 + e;
            const double zsspp2 = a * zsssp2 + u7;
            const double zppss2 = b * zspss2 + u8;
            const double zspsp3 = zsssp2 * zspss2;
            const double zsspp3 = zsssp2 * zsssp2;
            const double zppss3 = zspss2 * zspss2;
            const double zsppp1 = zsspp1 * zspss1;
            const double zpspp1 = zsspp1 * zpsss1;
            const double zppsp1 = zppss1 * zsssp1;
            const double zppps1 = zppss1 * zssps1;
            d = zspss2 * zsspp1 + a * u2;
            e = zsssp2 * zppss1 + b * u2;
            const double zsppp2 = zspss1 * zsspp2 + d;
            const double zpspp2 = zpsss1 * zsspp2 + d;
            const double zppsp2 = zsssp1 * zppss2 + e;
            const double zppps2 = zssps1 * zppss2 + e;
            d = a * zspsp3 + zsssp2 * u5;
            e = b * zspsp3 + zspss2 * u5;
            const double zsppp3 = zspss1 * zsspp3 + d;
            const double zpspp3 = zpsss1 * zsspp3 + d;
            const double zppsp3 = zsssp1 * zppss3 + e;
            const double zppps3 = zssps1 * zppss3 + e;
            const double zsppp4 = zsspp3 * zspss2;
            const double zppsp4 = zppss3 * zsssp2;
            d = b * zsssp2 + a * zspss2 + u2;
            const double zpppp1 = zsspp1 * zppss1;
            const double zpppp2 = zsspp2 * zppss1 + zppss2 * zsspp1 + c * u2;
            const double zpppp3 = zppss1 * zsspp3 + zsspp1 * zppss3 + zspsp3 * c + u5 * d;
            const double zpppp4 = zspsp3 * (d + u6);
            const double zpppp5 = zsspp3 * zppss3;

/*             ...assemble the 4-center (AB|CD) type integrals. */
            const double gxxxx = xpppp1 * f0 + xpppp2 * f1 + xpppp3 * f2 + xpppp4 * f3 + xpppp5 * f4;
            const double gyyyy = ypppp1 * f0 + ypppp2 * f1 + ypppp3 * f2 + ypppp4 * f3 + ypppp5 * f4;
            const double gzzzz = zpppp1 * f0 + zpppp2 * f1 + zpppp3 * f2 + zpppp4 * f3 + zpppp5 * f4;
            a = xppps1 * f0 + xppps2 * f1 + xppps3 * f2 + xppsp4 * f3;
            b = xppps1 * f1 + xppps2 * f2 + xppps3 * f3 + xppsp4 * f4;
            c = xppsp1 * f0 + xppsp2 * f1 + xppsp3 * f2 + xppsp4 * f3;
            d = xppsp1 * f1 + xppsp2 * f2 + xppsp3 * f3 + xppsp4 * f4;
            e = xpspp1 * f0 + xpspp2 * f1 + xpspp3 * f2 + xsppp4 * f3;
            double f = xpspp1 * f1 + xpspp2 * f2 + xpspp3 * f3 + xsppp4 * f4;
            double g = xsppp1 * f0 + xsppp2 * f1 + xsppp3 * f2 + xsppp4 * f3;
            double h = xsppp1 * f1 + xsppp2 * f2 + xsppp3 * f3 + xsppp4 * f4;
            const double gxxxy = a * ysssp1 + b * ysssp2;
            const double gxxyx = c * yssps1 + d * ysssp2;
            const double gxyxx = e * yspss1 + f * yspss2;
            const double gyxxx = g * ypsss1 + h * yspss2;
            const double gxxxz = a * zsssp1 + b * zsssp2;
            const double gxxzx = c * zssps1 + d * zsssp2;
            const double gxzxx = e * zspss1 + f * zspss2;
            const double gzxxx = g * zpsss1 + h * zspss2;
            a = yppps1 * f0 + yppps2 * f1 + yppps3 * f2 + yppsp4 * f3;
            b = yppps1 * f1 + yppps2 * f2 + yppps3 * f3 + yppsp4 * f4;
            c = yppsp1 * f0 + yppsp2 * f1 + yppsp3 * f2 + yppsp4 * f3;
            d = yppsp1 * f1 + yppsp2 * f2 + yppsp3 * f3 + yppsp4 * f4;
            e = ypspp1 * f0 + ypspp2 * f1 + ypspp3 * f2 + ysppp4 * f3;
            f = ypspp1 * f1 + ypspp2 * f2 + ypspp3 * f3 + ysppp4 * f4;
            g = ysppp1 * f0 + ysppp2 * f1 + ysppp3 * f2 + ysppp4 * f3;
            h = ysppp1 * f1 + ysppp2 * f2 + ysppp3 * f3 + ysppp4 * f4;
            const double gyyyx = a * xsssp1 + b * xsssp2;
            const double gyyxy = c * xssps1 + d * xsssp2;
            const double gyxyy = e * xspss1 + f * xspss2;
            const double gxyyy = g * xpsss1 + h * xspss2;
            const double gyyyz = a * zsssp1 + b * zsssp2;
            const double gyyzy = c * zssps1 + d * zsssp2;
            const double gyzyy = e * zspss1 + f * zspss2;
            const double gzyyy = g * zpsss1 + h * zspss2;
            a = zppps1 * f0 + zppps2 * f1 + zppps3 * f2 + zppsp4 * f3;
            b = zppps1 * f1 + zppps2 * f2 + zppps3 * f3 + zppsp4 * f4;
            c = zppsp1 * f0 + zppsp2 * f1 + zppsp3 * f2 + zppsp4 * f3;
            d = zppsp1 * f1 + zppsp2 * f2 + zppsp3 * f3 + zppsp4 * f4;
            e = zpspp1 * f0 + zpspp2 * f1 + zpspp3 * f2 + zsppp4 * f3;
            f = zpspp1 * f1 + zpspp2 * f2 + zpspp3 * f3 + zsppp4 * f4;
            g = zsppp1 * f0 + zsppp2 * f1 + zsppp3 * f2 + zsppp4 * f3;
            h = zsppp1 * f1 + zsppp2 * f2 + zsppp3 * f3 + zsppp4 * f4;
            const double gzzzx = a * xsssp1 + b * xsssp2;
            const double gzzxz = c * xssps1 + d * xsssp2;
            const double gzxzz = e * xspss1 + f * xspss2;
            const double gxzzz = g * xpsss1 + h * xspss2;
            const double gzzzy = a * ysssp1 + b * ysssp2;
            const double gzzyz = c * yssps1 + d * ysssp2;
            const double gzyzz = e * yspss1 + f * yspss2;
            const double gyzzz = g * ypsss1 + h * yspss2;
            a = xppss1 * f0 + xppss2 * f1 + xppss3 * f2;
            b = xppss1 * f1 + xppss2 * f2 + xppss3 * f3;
            c = xppss1 * f2 + xppss2 * f3 + xppss3 * f4;
            d = yppss1 * f0 + yppss2 * f1 + yppss3 * f2;
            e = yppss1 * f1 + yppss2 * f2 + yppss3 * f3;
            f = yppss1 * f2 + yppss2 * f3 + yppss3 * f4;
            g = zppss1 * f0 + zppss2 * f1 + zppss3 * f2;
            h = zppss1 * f1 + zppss2 * f2 + zppss3 * f3;
            double r = zppss1 * f2 + zppss2 * f3 + zppss3 * f4;
            const double gxxyy = a * ysspp1 + b * ysspp2 + c * ysspp3;
            const double gxxzz = a * zsspp1 + b * zsspp2 + c * zsspp3;
            const double gyyxx = d * xsspp1 + e * xsspp2 + f * xsspp3;
            const double gyyzz = d * zsspp1 + e * zsspp2 + f * zsspp3;
            const double gzzxx = g * xsspp1 + h * xsspp2 + r * xsspp3;
            const double gzzyy = g * ysspp1 + h * ysspp2 + r * ysspp3;
            const double gxxyz = (a * yssps1 + b * ysssp2) * zsssp1 + (b * yssps1 + c * ysssp2) * zsssp2;
            const double gxxzy = (a * zssps1 + b * zsssp2) * ysssp1 + (b * zssps1 + c * zsssp2) * ysssp2;
            const double gyyxz = (d * xssps1 + e * xsssp2) * zsssp1 + (e * xssps1 + f * xsssp2) * zsssp2;
            const double gyyzx = (d * zssps1 + e * zsssp2) * xsssp1 + (e * zssps1 + f * zsssp2) * xsssp2;
            const double gzzxy = (g * xssps1 + h * xsssp2) * ysssp1 + (h * xssps1 + r * xsssp2) * ysssp2;
            const double gzzyx = (g * yssps1 + h * ysssp2) * xsssp1 + (h * yssps1 + r * ysssp2) * xsssp2;
            double aa = xspsp3 * f2;
            double bb = xspsp3 * f3;
            double cc = xspsp3 * f4;
            double dd = yspsp3 * f2;
            double ee = yspsp3 * f3;
            double ff = yspsp3 * f4;
            double gg = zspsp3 * f2;
            double hh = zspsp3 * f3;
            double rr = zspsp3 * f4;
            a = xpsps1 * f0 + xpsps2 * f1 + aa;
            b = xpsps1 * f1 + xpsps2 * f2 + bb;
            c = xpsps1 * f2 + xpsps2 * f3 + cc;
            d = ypsps1 * f0 + ypsps2 * f1 + dd;
            e = ypsps1 * f1 + ypsps2 * f2 + ee;
            f = ypsps1 * f2 + ypsps2 * f3 + ff;
            g = zpsps1 * f0 + zpsps2 * f1 + gg;
            h = zpsps1 * f1 + zpsps2 * f2 + hh;
            r = zpsps1 * f2 + zpsps2 * f3 + rr;
            const double gxyxy = a * yspsp1 + b * yspsp2 + c * yspsp3;
            const double gxzxz = a * zspsp1 + b * zspsp2 + c * zspsp3;
            const double gyxyx = d * xspsp1 + e * xspsp2 + f * xspsp3;
            const double gyzyz = d * zspsp1 + e * zspsp2 + f * zspsp3;
            const double gzxzx = g * xspsp1 + h * xspsp2 + r * xspsp3;
            const double gzyzy = g * yspsp1 + h * yspsp2 + r * yspsp3;
            const double gxyxz = (a * yspss1 + b * yspss2) * zsssp1 + (b * yspss1 + c * yspss2) * zsssp2;
            const double gxzxy = (a * zspss1 + b * zspss2) * ysssp1 + (b * zspss1 + c * zspss2) * ysssp2;
            const double gyxyz = (d * xspss1 + e * xspss2) * zsssp1 + (e * xspss1 + f * xspss2) * zsssp2;
            const double gyzyx = (d * zspss1 + e * zspss2) * xsssp1 + (e * zspss1 + f * zspss2) * xsssp2;
            const double gzxzy = (g * xspss1 + h * xspss2) * ysssp1 + (h * xspss1 + r * xspss2) * ysssp2;
            const double gzyzx = (g * yspss1 + h * yspss2) * xsssp1 + (h * yspss1 + r * yspss2) * xsssp2;
            a = xpssp1 * f0 + xpssp2 * f1 + aa;
            b = xpssp1 * f1 + xpssp2 * f2 + bb;
            c = xpssp1 * f2 + xpssp2 * f3 + cc;
            d = ypssp1 * f0 + ypssp2 * f1 + dd;
            e = ypssp1 * f1 + ypssp2 * f2 + ee;
            f = ypssp1 * f2 + ypssp2 * f3 + ff;
            g = zpssp1 * f0 + zpssp2 * f1 + gg;
            h = zpssp1 * f1 + zpssp2 * f2 + hh;
            r = zpssp1 * f2 + zpssp2 * f3 + rr;
            const double gxyyx = a * yspps1 + b * yspps2 + c * yspsp3;
            const double gxzzx = a * zspps1 + b * zspps2 + c * zspsp3;
            const double gyxxy = d * xspps1 + e * xspps2 + f * xspsp3;
            const double gyzzy = d * zspps1 + e * zspps2 + f * zspsp3;
            const double gzxxz = g * xspps1 + h * xspps2 + r * xspsp3;
            const double gzyyz = g * yspps1 + h * yspps2 + r * yspsp3;
            const double gxyzx = (a * yspss1 + b * yspss2) * zssps1 + (b * yspss1 + c * yspss2) * zsssp2;
            const double gxzyx = (a * zspss1 + b * zspss2) * yssps1 + (b * zspss1 + c * zspss2) * ysssp2;
            const double gyxzy = (d * xspss1 + e * xspss2) * zssps1 + (e * xspss1 + f * xspss2) * zsssp2;
            const double gyzxy = (d * zspss1 + e * zspss2) * xssps1 + (e * zspss1 + f * zspss2) * xsssp2;
            const double gzxyz = (g * xspss1 + h * xspss2) * yssps1 + (h * xspss1 + r * xspss2) * ysssp2;
            const double gzyxz = (g * yspss1 + h * yspss2) * xssps1 + (h * yspss1 + r * yspss2) * xsssp2;
            a = xspps1 * f0 + xspps2 * f1 + aa;
            b = xspps1 * f1 + xspps2 * f2 + bb;
            c = xspps1 * f2 + xspps2 * f3 + cc;
            d = yspps1 * f0 + yspps2 * f1 + dd;
            e = yspps1 * f1 + yspps2 * f2 + ee;
            f = yspps1 * f2 + yspps2 * f3 + ff;
            g = zspps1 * f0 + zspps2 * f1 + gg;
            h = zspps1 * f1 + zspps2 * f2 + hh;
            r = zspps1 * f2 + zspps2 * f3 + rr;
            const double gyxxz = (a * ypsss1 + b * yspss2) * zsssp1 + (b * ypsss1 + c * yspss2) * zsssp2;
            const double gzxxy = (a * zpsss1 + b * zspss2) * ysssp1 + (b * zpsss1 + c * zspss2) * ysssp2;
            const double gxyyz = (d * xpsss1 + e * xspss2) * zsssp1 + (e * xpsss1 + f * xspss2) * zsssp2;
            const double gzyyx = (d * zpsss1 + e * zspss2) * xsssp1 + (e * zpsss1 + f * zspss2) * xsssp2;
            const double gxzzy = (g * xpsss1 + h * xspss2) * ysssp1 + (h * xpsss1 + r * xspss2) * ysssp2;
            const double gyzzx = (g * ypsss1 + h * yspss2) * xsssp1 + (h * ypsss1 + r * yspss2) * xsssp2;
            a = xspsp1 * f0 + xspsp2 * f1 + aa;
            b = xspsp1 * f1 + xspsp2 * f2 + bb;
            c = xspsp1 * f2 + xspsp2 * f3 + cc;
            d = yspsp1 * f0 + yspsp2 * f1 + dd;
            e = yspsp1 * f1 + yspsp2 * f2 + ee;
            f = yspsp1 * f2 + yspsp2 * f3 + ff;
            g = zspsp1 * f0 + zspsp2 * f1 + gg;
            h = zspsp1 * f1 + zspsp2 * f2 + hh;
            r = zspsp1 * f2 + zspsp2 * f3 + rr;
            const double gyxzx = (a * ypsss1 + b * yspss2) * zssps1 + (b * ypsss1 + c * yspss2) * zsssp2;
            const double gzxyx = (a * zpsss1 + b * zspss2) * yssps1 + (b * zpsss1 + c * zspss2) * ysssp2;
            const double gxyzy = (d * xpsss1 + e * xspss2) * zssps1 + (e * xpsss1 + f * xspss2) * zsssp2;
            const double gzyxy = (d * zpsss1 + e * zspss2) * xssps1 + (e * zpsss1 + f * zspss2) * xsssp2;
            const double gxzyz = (g * xpsss1 + h * xspss2) * yssps1 + (h * xpsss1 + r * xspss2) * ysssp2;
            const double gyzxz = (g * ypsss1 + h * yspss2) * xssps1 + (h * ypsss1 + r * yspss2) * xsssp2;
            a = xsspp1 * f0 + xsspp2 * f1 + xsspp3 * f2;
            b = xsspp1 * f1 + xsspp2 * f2 + xsspp3 * f3;
            c = xsspp1 * f2 + xsspp2 * f3 + xsspp3 * f4;
            d = ysspp1 * f0 + ysspp2 * f1 + ysspp3 * f2;
            e = ysspp1 * f1 + ysspp2 * f2 + ysspp3 * f3;
            f = ysspp1 * f2 + ysspp2 * f3 + ysspp3 * f4;
            g = zsspp1 * f0 + zsspp2 * f1 + zsspp3 * f2;
            h = zsspp1 * f1 + zsspp2 * f2 + zsspp3 * f3;
            r = zsspp1 * f2 + zsspp2 * f3 + zsspp3 * f4;
            const double gyzxx = (a * ypsss1 + b * yspss2) * zspss1 + (b * ypsss1 + c * yspss2) * zspss2;
            const double gzyxx = (a * zpsss1 + b * zspss2) * yspss1 + (b * zpsss1 + c * zspss2) * yspss2;
            const double gxzyy = (d * xpsss1 + e * xspss2) * zspss1 + (e * xpsss1 + f * xspss2) * zspss2;
            const double gzxyy = (d * zpsss1 + e * zspss2) * xspss1 + (e * zpsss1 + f * zspss2) * xspss2;
            const double gxyzz = (g * xpsss1 + h * xspss2) * yspss1 + (h * xpsss1 + r * xspss2) *yspss2;
            const double gyxzz = (g * ypsss1 + h * yspss2) * xspss1 + (h * ypsss1 + r * yspss2) * xspss2;
            cbatch[0] += gxxxx;
            cbatch[1] += gyxxx;
            cbatch[2] += gzxxx;
            cbatch[3] += gxyxx;
            cbatch[4] += gyyxx;
            cbatch[5] += gzyxx;
            cbatch[6] += gxzxx;
            cbatch[7] += gyzxx;
            cbatch[8] += gzzxx;
            cbatch[9] += gxxyx;
            cbatch[10] += gyxyx;
            cbatch[11] += gzxyx;
            cbatch[12] += gxyyx;
            cbatch[13] += gyyyx;
            cbatch[14] += gzyyx;
            cbatch[15] += gxzyx;
            cbatch[16] += gyzyx;
            cbatch[17] += gzzyx;
            cbatch[18] += gxxzx;
            cbatch[19] += gyxzx;
            cbatch[20] += gzxzx;
            cbatch[21] += gxyzx;
            cbatch[22] += gyyzx;
            cbatch[23] += gzyzx;
            cbatch[24] += gxzzx;
            cbatch[25] += gyzzx;
            cbatch[26] += gzzzx;
            cbatch[27] += gxxxy;
            cbatch[28] += gyxxy;
            cbatch[29] += gzxxy;
            cbatch[30] += gxyxy;
            cbatch[31] += gyyxy;
            cbatch[32] += gzyxy;
            cbatch[33] += gxzxy;
            cbatch[34] += gyzxy;
            cbatch[35] += gzzxy;
            cbatch[36] += gxxyy;
            cbatch[37] += gyxyy;
            cbatch[38] += gzxyy;
            cbatch[39] += gxyyy;
            cbatch[40] += gyyyy;
            cbatch[41] += gzyyy;
            cbatch[42] += gxzyy;
            cbatch[43] += gyzyy;
            cbatch[44] += gzzyy;
            cbatch[45] += gxxzy;
            cbatch[46] += gyxzy;
            cbatch[47] += gzxzy;
            cbatch[48] += gxyzy;
            cbatch[49] += gyyzy;
            cbatch[50] += gzyzy;
            cbatch[51] += gxzzy;
            cbatch[52] += gyzzy;
            cbatch[53] += gzzzy;
            cbatch[54] += gxxxz;
            cbatch[55] += gyxxz;
            cbatch[56] += gzxxz;
            cbatch[57] += gxyxz;
            cbatch[58] += gyyxz;
            cbatch[59] += gzyxz;
            cbatch[60] += gxzxz;
            cbatch[61] += gyzxz;
            cbatch[62] += gzzxz;
            cbatch[63] += gxxyz;
            cbatch[64] += gyxyz;
            cbatch[65] += gzxyz;
            cbatch[66] += gxyyz;
            cbatch[67] += gyyyz;
            cbatch[68] += gzyyz;
            cbatch[69] += gxzyz;
            cbatch[70] += gyzyz;
            cbatch[71] += gzzyz;
            cbatch[72] += gxxzz;
            cbatch[73] += gyxzz;
            cbatch[74] += gzxzz;
            cbatch[75] += gxyzz;
            cbatch[76] += gyyzz;
            cbatch[77] += gzyzz;
            cbatch[78] += gxzzz;
            cbatch[79] += gyzzz;
            cbatch[80] += gzzzz;
        }
    }
}