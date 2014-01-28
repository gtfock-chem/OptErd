#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <yepPredefines.h>
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
void erd__sppp_pcgto_block (int nij, int nkl,
                           int shell1, int shell3, int shellp,
                           double x1, double y1, double z1,
                           double x2, double y2, double z2,
                           double x3, double y3, double z3,
                           double x4, double y4, double z4,
                           double *YEP_RESTRICT alpha1, double *YEP_RESTRICT alpha2,
                           double *YEP_RESTRICT alpha3, double *YEP_RESTRICT alpha4,
                           double *YEP_RESTRICT cc1, double *YEP_RESTRICT cc2,
                           double *YEP_RESTRICT cc3, double *YEP_RESTRICT cc4,
                           int *YEP_RESTRICT prim1, int *YEP_RESTRICT prim2,
                           int *YEP_RESTRICT prim3, int *YEP_RESTRICT prim4,
                           double *YEP_RESTRICT norm1, double *YEP_RESTRICT norm2,
                           double *YEP_RESTRICT norm3, double *YEP_RESTRICT norm4,
                           double *YEP_RESTRICT rho12, double *YEP_RESTRICT rho34,
                           double *YEP_RESTRICT p, double *YEP_RESTRICT px,
                           double *YEP_RESTRICT py, double *YEP_RESTRICT pz, double *YEP_RESTRICT scalep,
                           double *YEP_RESTRICT q, double *YEP_RESTRICT qx,
                           double *YEP_RESTRICT qy, double *YEP_RESTRICT qz, double *YEP_RESTRICT scaleq, double *YEP_RESTRICT cbatch)
{
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

    if (shellp == 1) {
        // 1     5   |  (AB|CD)  4-center   sppp and pspp
        double pxsub, pysub, pzsub;
        if (shell1 == 1) {
            pxsub = x1;
            pysub = y1;
            pzsub = z1;
        } else {
            pxsub = x2;
            pysub = y2;
            pzsub = z2;
        }
        for (int ij = 0; ij < nij; ij += 1) {
            const double pval = p[ij];
            const double pxval = px[ij];
            const double pyval = py[ij];
            const double pzval = pz[ij];
            const double pscale = scalep[ij];
            const double xpss1 = pxval - pxsub;
            const double ypss1 = pyval - pysub;
            const double zpss1 = pzval - pzsub;
            for (int kl = 0; kl < nkl; kl += 1) {
                const double qval = q[kl];
                const double qinv = 1. / qval;
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
                const double f0 = scale * boys0 (t);
                const double f1 = scale * boys1 (t);
                const double f2 = scale * boys2 (t);
                const double f3 = scale * boys3 (t);
#else
                double f0, f1, f2, f3;
                if (t <= tmax) {
                    const int tgrid = __builtin_lround(t * tvstep);
                    const double delta1 = tgrid * tstep - t;
                    const double delta2 = delta1 * .5;
                    const double delta3 = delta1 * .333333333333333;
                    const double delta4 = delta2 * .5;
                    const double delta5 = delta1 * .2;
                    const double delta6 = delta3 * .5;
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
                    f0 = scale * f0;
                    f1 = scale * f1;
                    f2 = scale * f2;
                    f3 = scale * f3;
                } else {
                    const double tinv = 1. / t;
                    const double t2inv = tinv * .5;
                    f0 = scale * .5 * __builtin_sqrt(tinv * 3.141592653589793);
                    f1 = t2inv * f0;
                    f2 = t2inv * 3. * f1;
                    f3 = t2inv * 5. * f2;
                }
#endif
                const double u0 = pval * pqpinv;
                const double u1 = -qval * pqpinv;
                const double u2 = pqpinv * .5;
                const double u3 = u2 + pqpinv;
                const double u4 = qinv * .5;
                const double u5 = u0 * u4;

/*             ...the X-terms. */
                const double xssp1 = qxval - x4;
                const double xsps1 = qxval - x3;
                const double xssp2 = pqx * u0;
                const double xpss2 = pqx * u1;
                double a = xsps1 + xssp1;
                double b = xpss1 * xssp2 + u2;
                const double xspp1 = xsps1 * xssp1 + u4;
                const double xspp2 = a * xssp2 - u5;
                const double xspp3 = xssp2 * xssp2;
                const double xpsp1 = xpss1 * xssp1;
                const double xpsp2 = xssp1 * xpss2 + b;
                const double xpsp3 = xssp2 * xpss2;
                const double xpps1 = xpss1 * xsps1;
                const double xpps2 = xsps1 * xpss2 + b;
                const double xppp1 = xpss1 * xspp1;
                const double xppp2 = xpss1 * xspp2 + xspp1 * xpss2 + a * u2;
                const double xppp3 = xpss1 * xspp3 + a * xpsp3 + u3 * xssp2;
                const double xppp4 = xpss2 * xspp3;

/*             ...the Y-terms. */
                const double yssp1 = qyval - y4;
                const double ysps1 = qyval - y3;
                const double yssp2 = pqy * u0;
                const double ypss2 = pqy * u1;
                a = ysps1 + yssp1;
                b = ypss1 * yssp2 + u2;
                const double yspp1 = ysps1 * yssp1 + u4;
                const double yspp2 = a * yssp2 - u5;
                const double yspp3 = yssp2 * yssp2;
                const double ypsp1 = ypss1 * yssp1;
                const double ypsp2 = yssp1 * ypss2 + b;
                const double ypsp3 = yssp2 * ypss2;
                const double ypps1 = ypss1 * ysps1;
                const double ypps2 = ysps1 * ypss2 + b;
                const double yppp1 = ypss1 * yspp1;
                const double yppp2 = ypss1 * yspp2 + yspp1 * ypss2 + a * u2;
                const double yppp3 = ypss1 * yspp3 + a * ypsp3 + u3 * yssp2;
                const double yppp4 = ypss2 * yspp3;

/*             ...the Z-terms. */
                const double zssp1 = qzval - z4;
                const double zsps1 = qzval - z3;
                const double zssp2 = pqz * u0;
                const double zpss2 = pqz * u1;
                a = zsps1 + zssp1;
                b = zpss1 * zssp2 + u2;
                const double zspp1 = zsps1 * zssp1 + u4;
                const double zspp2 = a * zssp2 - u5;
                const double zspp3 = zssp2 * zssp2;
                const double zpsp1 = zpss1 * zssp1;
                const double zpsp2 = zssp1 * zpss2 + b;
                const double zpsp3 = zssp2 * zpss2;
                const double zpps1 = zpss1 * zsps1;
                const double zpps2 = zsps1 * zpss2 + b;
                const double zppp1 = zpss1 * zspp1;
                const double zppp2 = zpss1 * zspp2 + zspp1 * zpss2 + a * u2;
                const double zppp3 = zpss1 * zspp3 + a * zpsp3 + u3 * zssp2;
                const double zppp4 = zpss2 * zspp3;

/*             ...assemble the 4-center (AB|CD) type integrals. */
                const double gxxx = xppp1 * f0 + xppp2 * f1 + xppp3 * f2 + xppp4 * f3;
                const double gyyy = yppp1 * f0 + yppp2 * f1 + yppp3 * f2 + yppp4 * f3;
                const double gzzz = zppp1 * f0 + zppp2 * f1 + zppp3 * f2 + zppp4 * f3;
                double aa = xpsp3 * f2;
                double bb = xpsp3 * f3;
                a = xpps1 * f0 + xpps2 * f1 + aa;
                b = xpps1 * f1 + xpps2 * f2 + bb;
                double c = xpsp1 * f0 + xpsp2 * f1 + aa;
                double d = xpsp1 * f1 + xpsp2 * f2 + bb;
                double e = xspp1 * f0 + xspp2 * f1 + xspp3 * f2;
                double f = xspp1 * f1 + xspp2 * f2 + xspp3 * f3;
                const double gxxy = yssp1 * a + yssp2 * b;
                const double gxxz = zssp1 * a + zssp2 * b;
                const double gxyx = ysps1 * c + yssp2 * d;
                const double gxzx = zsps1 * c + zssp2 * d;
                const double gyxx = ypss1 * e + ypss2 * f;
                const double gzxx = zpss1 * e + zpss2 * f;
                aa = ypsp3 * f2;
                bb = ypsp3 * f3;
                a = ypps1 * f0 + ypps2 * f1 + aa;
                b = ypps1 * f1 + ypps2 * f2 + bb;
                c = ypsp1 * f0 + ypsp2 * f1 + aa;
                d = ypsp1 * f1 + ypsp2 * f2 + bb;
                e = yspp1 * f0 + yspp2 * f1 + yspp3 * f2;
                f = yspp1 * f1 + yspp2 * f2 + yspp3 * f3;
                const double gyyx = xssp1 * a + xssp2 * b;
                const double gyyz = zssp1 * a + zssp2 * b;
                const double gyxy = xsps1 * c + xssp2 * d;
                const double gyzy = zsps1 * c + zssp2 * d;
                const double gxyy = xpss1 * e + xpss2 * f;
                const double gzyy = zpss1 * e + zpss2 * f;
                aa = zpsp3 * f2;
                bb = zpsp3 * f3;
                a = zpps1 * f0 + zpps2 * f1 + aa;
                b = zpps1 * f1 + zpps2 * f2 + bb;
                c = zpsp1 * f0 + zpsp2 * f1 + aa;
                d = zpsp1 * f1 + zpsp2 * f2 + bb;
                e = zspp1 * f0 + zspp2 * f1 + zspp3 * f2;
                f = zspp1 * f1 + zspp2 * f2 + zspp3 * f3;
                const double gzzx = xssp1 * a + xssp2 * b;
                const double gzzy = yssp1 * a + yssp2 * b;
                const double gzxz = xsps1 * c + xssp2 * d;
                const double gzyz = ysps1 * c + yssp2 * d;
                const double gxzz = xpss1 * e + xpss2 * f;
                const double gyzz = ypss1 * e + ypss2 * f;
                a = xpss1 * f0 + xpss2 * f1;
                b = xpss1 * f1 + xpss2 * f2;
                c = xpss1 * f2 + xpss2 * f3;
                d = ypss1 * f0 + ypss2 * f1;
                e = ypss1 * f1 + ypss2 * f2;
                f = ypss1 * f2 + ypss2 * f3;
                double g = zpss1 * f0 + zpss2 * f1;
                double h = zpss1 * f1 + zpss2 * f2;
                double r = zpss1 * f2 + zpss2 * f3;
                const double gxyz = ysps1 * (zssp1 * a + zssp2 * b) + yssp2 * (zssp1 * b + zssp2 * c);
                const double gxzy = zsps1 * (yssp1 * a + yssp2 * b) + zssp2 * (yssp1 * b + yssp2 * c);
                const double gyxz = xsps1 * (zssp1 * d + zssp2 * e) + xssp2 * (zssp1 * e + zssp2 * f);
                const double gyzx = zsps1 * (xssp1 * d + xssp2 * e) + zssp2 * (xssp1 * e + xssp2 * f);
                const double gzxy = xsps1 * (yssp1 * g + yssp2 * h) + xssp2 * (yssp1 * h + yssp2 * r);
                const double gzyx = ysps1 * (xssp1 * g + xssp2 * h) + yssp2 * (xssp1 * h + xssp2 * r);
                cbatch[0] += gxxx;
                cbatch[1] += gyxx;
                cbatch[2] += gzxx;
                cbatch[3] += gxyx;
                cbatch[4] += gyyx;
                cbatch[5] += gzyx;
                cbatch[6] += gxzx;
                cbatch[7] += gyzx;
                cbatch[8] += gzzx;
                cbatch[9] += gxxy;
                cbatch[10] += gyxy;
                cbatch[11] += gzxy;
                cbatch[12] += gxyy;
                cbatch[13] += gyyy;
                cbatch[14] += gzyy;
                cbatch[15] += gxzy;
                cbatch[16] += gyzy;
                cbatch[17] += gzzy;
                cbatch[18] += gxxz;
                cbatch[19] += gyxz;
                cbatch[20] += gzxz;
                cbatch[21] += gxyz;
                cbatch[22] += gyyz;
                cbatch[23] += gzyz;
                cbatch[24] += gxzz;
                cbatch[25] += gyzz;
                cbatch[26] += gzzz;
            }
        }
    } else {
        // 2     5   |  (AB|CD)  4-center   ppsp and ppps
        double qxsub, qysub, qzsub;
        if (shell3 == 1) {
            qxsub = x3;
            qysub = y3;
            qzsub = z3;
        } else {
            qxsub = x4;
            qysub = y4;
            qzsub = z4;
        }
        for (int ij = 0; ij < nij; ij += 1) {
            const double pval = p[ij];
            const double pinv = 1. / pval;
            const double pxval = px[ij];
            const double pyval = py[ij];
            const double pzval = pz[ij];
            const double pscale = scalep[ij];
            const double u4 = pinv * .5;
            const double xsps1 = pxval - x2;
            const double xpss1 = pxval - x1;
            const double ysps1 = pyval - y2;
            const double ypss1 = pyval - y1;
            const double zsps1 = pzval - z2;
            const double zpss1 = pzval - z1;
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
                const double f0 = scale * boys0 (t);
                const double f1 = scale * boys1 (t);
                const double f2 = scale * boys2 (t);
                const double f3 = scale * boys3 (t);
#else
                double f0, f1, f2, f3;
                if (t <= tmax) {
                    const int tgrid = __builtin_lround(t * tvstep);
                    const double delta1 = tgrid * tstep - t;
                    const double delta2 = delta1 * .5;
                    const double delta3 = delta1 * .333333333333333;
                    const double delta4 = delta2 * .5;
                    const double delta5 = delta1 * .2;
                    const double delta6 = delta3 * .5;
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
                    f0 = scale * f0;
                    f1 = scale * f1;
                    f2 = scale * f2;
                    f3 = scale * f3;
                } else {
                    const double tinv = 1. / t;
                    const double t2inv = tinv * .5;
                    f0 = scale * .5 * __builtin_sqrt(tinv * 3.141592653589793);
                    f1 = t2inv * f0;
                    f2 = t2inv * 3. * f1;
                    f3 = t2inv * 5. * f2;
                }
#endif
                const double u0 = pval * pqpinv;
                const double u1 = -qval * pqpinv;
                const double u2 = pqpinv * .5;
                const double u3 = u2 + pqpinv;
                const double u5 = u1 * u4;

/*             ...the X-terms. */
                const double xssp1 = qxval - qxsub;
                const double xssp2 = pqx * u0;
                const double xsps2 = pqx * u1;
                double a = xpss1 + xsps1;
                double b = xssp1 * xsps2 + u2;
                const double xspp1 = xsps1 * xssp1;
                const double xspp2 = xsps1 * xssp2 + b;
                const double xspp3 = xssp2 * xsps2;
                const double xpsp1 = xpss1 * xssp1;
                const double xpsp2 = xpss1 * xssp2 + b;
                const double xpps1 = xpss1 * xsps1 + u4;
                const double xpps2 = a * xsps2 + u5;
                const double xpps3 = xsps2 * xsps2;
                const double xppp1 = xssp1 * xpps1;
                const double xppp2 = xssp1 * xpps2 + xpps1 * xssp2 + a * u2;
                const double xppp3 = xssp1 * xpps3 + a * xspp3 + u3 * xsps2;
                const double xppp4 = xssp2 * xpps3;

/*             ...the Y-terms. */
                const double yssp1 = qyval - qysub;
                const double yssp2 = pqy * u0;
                const double ysps2 = pqy * u1;
                a = ypss1 + ysps1;
                b = yssp1 * ysps2 + u2;
                const double yspp1 = ysps1 * yssp1;
                const double yspp2 = ysps1 * yssp2 + b;
                const double yspp3 = yssp2 * ysps2;
                const double ypsp1 = ypss1 * yssp1;
                const double ypsp2 = ypss1 * yssp2 + b;
                const double ypps1 = ypss1 * ysps1 + u4;
                const double ypps2 = a * ysps2 + u5;
                const double ypps3 = ysps2 * ysps2;
                const double yppp1 = yssp1 * ypps1;
                const double yppp2 = yssp1 * ypps2 + ypps1 * yssp2 + a * u2;
                const double yppp3 = yssp1 * ypps3 + a * yspp3 + u3 * ysps2;
                const double yppp4 = yssp2 * ypps3;

/*             ...the Z-terms. */
                const double zssp1 = qzval - qzsub;
                const double zssp2 = pqz * u0;
                const double zsps2 = pqz * u1;
                a = zpss1 + zsps1;
                b = zssp1 * zsps2 + u2;
                const double zspp1 = zsps1 * zssp1;
                const double zspp2 = zsps1 * zssp2 + b;
                const double zspp3 = zssp2 * zsps2;
                const double zpsp1 = zpss1 * zssp1;
                const double zpsp2 = zpss1 * zssp2 + b;
                const double zpps1 = zpss1 * zsps1 + u4;
                const double zpps2 = a * zsps2 + u5;
                const double zpps3 = zsps2 * zsps2;
                const double zppp1 = zssp1 * zpps1;
                const double zppp2 = zssp1 * zpps2 + zpps1 * zssp2 + a * u2;
                const double zppp3 = zssp1 * zpps3 + a * zspp3 + u3 * zsps2;
                const double zppp4 = zssp2 * zpps3;

/*             ...assemble the 4-center (AB|CD) type integrals. */
                const double gxxx = xppp1 * f0 + xppp2 * f1 + xppp3 * f2 + xppp4 * f3;
                const double gyyy = yppp1 * f0 + yppp2 * f1 + yppp3 * f2 + yppp4 * f3;
                const double gzzz = zppp1 * f0 + zppp2 * f1 + zppp3 * f2 + zppp4 * f3;
                double aa = xspp3 * f2;
                double bb = xspp3 * f3;
                a = xpps1 * f0 + xpps2 * f1 + xpps3 * f2;
                b = xpps1 * f1 + xpps2 * f2 + xpps3 * f3;
                double c = xpsp1 * f0 + xpsp2 * f1 + aa;
                double d = xpsp1 * f1 + xpsp2 * f2 + bb;
                double e = xspp1 * f0 + xspp2 * f1 + aa;
                double f = xspp1 * f1 + xspp2 * f2 + bb;
                const double gxxy = yssp1 * a + yssp2 * b;
                const double gxxz = zssp1 * a + zssp2 * b;
                const double gxyx = ysps1 * c + ysps2 * d;
                const double gxzx = zsps1 * c + zsps2 * d;
                const double gyxx = ypss1 * e + ysps2 * f;
                const double gzxx = zpss1 * e + zsps2 * f;
                aa = yspp3 * f2;
                bb = yspp3 * f3;
                a = ypps1 * f0 + ypps2 * f1 + ypps3 * f2;
                b = ypps1 * f1 + ypps2 * f2 + ypps3 * f3;
                c = ypsp1 * f0 + ypsp2 * f1 + aa;
                d = ypsp1 * f1 + ypsp2 * f2 + bb;
                e = yspp1 * f0 + yspp2 * f1 + aa;
                f = yspp1 * f1 + yspp2 * f2 + bb;
                const double gyyx = xssp1 * a + xssp2 * b;
                const double gyyz = zssp1 * a + zssp2 * b;
                const double gyxy = xsps1 * c + xsps2 * d;
                const double gyzy = zsps1 * c + zsps2 * d;
                const double gxyy = xpss1 * e + xsps2 * f;
                const double gzyy = zpss1 * e + zsps2 * f;
                aa = zspp3 * f2;
                bb = zspp3 * f3;
                a = zpps1 * f0 + zpps2 * f1 + zpps3 * f2;
                b = zpps1 * f1 + zpps2 * f2 + zpps3 * f3;
                c = zpsp1 * f0 + zpsp2 * f1 + aa;
                d = zpsp1 * f1 + zpsp2 * f2 + bb;
                e = zspp1 * f0 + zspp2 * f1 + aa;
                f = zspp1 * f1 + zspp2 * f2 + bb;
                const double gzzx = xssp1 * a + xssp2 * b;
                const double gzzy = yssp1 * a + yssp2 * b;
                const double gzxz = xsps1 * c + xsps2 * d;
                const double gzyz = ysps1 * c + ysps2 * d;
                const double gxzz = xpss1 * e + xsps2 * f;
                const double gyzz = ypss1 * e + ysps2 * f;
                a = xpss1 * f0 + xsps2 * f1;
                b = xpss1 * f1 + xsps2 * f2;
                c = xpss1 * f2 + xsps2 * f3;
                d = ypss1 * f0 + ysps2 * f1;
                e = ypss1 * f1 + ysps2 * f2;
                f = ypss1 * f2 + ysps2 * f3;
                double g = zpss1 * f0 + zsps2 * f1;
                double h = zpss1 * f1 + zsps2 * f2;
                double r = zpss1 * f2 + zsps2 * f3;
                const double gxyz = ysps1 * (zssp1 * a + zssp2 * b) + ysps2 * (zssp1 * b + zssp2 * c);
                const double gxzy = zsps1 * (yssp1 * a + yssp2 * b) + zsps2 * (yssp1 * b + yssp2 * c);
                const double gyxz = xsps1 * (zssp1 * d + zssp2 * e) + xsps2 * (zssp1 * e + zssp2 * f);
                const double gyzx = zsps1 * (xssp1 * d + xssp2 * e) + zsps2 * (xssp1 * e + xssp2 * f);
                const double gzxy = xsps1 * (yssp1 * g + yssp2 * h) + xsps2 * (yssp1 * h + yssp2 * r);
                const double gzyx = ysps1 * (xssp1 * g + xssp2 * h) + ysps2 * (xssp1 * h + xssp2 * r);
                cbatch[0] += gxxx;
                cbatch[1] += gyxx;
                cbatch[2] += gzxx;
                cbatch[3] += gxyx;
                cbatch[4] += gyyx;
                cbatch[5] += gzyx;
                cbatch[6] += gxzx;
                cbatch[7] += gyzx;
                cbatch[8] += gzzx;
                cbatch[9] += gxxy;
                cbatch[10] += gyxy;
                cbatch[11] += gzxy;
                cbatch[12] += gxyy;
                cbatch[13] += gyyy;
                cbatch[14] += gzyy;
                cbatch[15] += gxzy;
                cbatch[16] += gyzy;
                cbatch[17] += gzzy;
                cbatch[18] += gxxz;
                cbatch[19] += gyxz;
                cbatch[20] += gzxz;
                cbatch[21] += gxyz;
                cbatch[22] += gyyz;
                cbatch[23] += gzyz;
                cbatch[24] += gxzz;
                cbatch[25] += gyzz;
                cbatch[26] += gzzz;
            }
        }
    }
}
