#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <yepPredefines.h>
#include "boys.h"
//#define ERD_TABLE_FREE_BOYS_FUNCTIONS

#pragma offload_attribute(push, target(mic))

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__SSSP_PCGTO_BLOCK */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation is designed to provide ultrafast block */
/*                evaluation of a batch of normalized electron repulsion */
/*                integrals between s-shell and p-shell primitive */
/*                spherical gaussian type orbitals. */
/*                A batch is defined here as containing all possible */
/*                integrals, that is its dimension is determined by */
/*                the total number of primitive functions (here = 3) */
/*                times the total number of ij and kl exponent pair */
/*                combinations. */
/*                The integrals are ordered in the batch the following */
/*                way (first index varying fastest): */
/*                    batch (nxyz1,nxyz2,nxyz3,nxyz4,kl,ij) */
/*                where ij and kl indicates alpha exponent pairs */
/*                defining the present block. */
/*                The present routine evaluates batches of the type: */
/*                        sssp , ssps , spss , psss */
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
/*                                    sssp/ssps/spss/psss integral batch */
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
/*                                    cartesian sssp/ssps/spss/psss */
/*                                    integrals */
/* ------------------------------------------------------------------------ */
void erd__sssp_pcgto_block (int nij, int nkl,
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

    if (shellp == 0) {
        // 0     5   |  (AB|CD)  4-center   sssp and ssps
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
            const double pxval = px[ij];
            const double pyval = py[ij];
            const double pzval = pz[ij];
            const double pscale = scalep[ij];
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
                double f1 = scale * boys1 (t);
#else
                double f0, f1;
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
                    f0 = scale * f0;
                    f1 = scale * f1;
                } else {
                    const double tinv = 1. / t;
                    f0 = scale * .5 * __builtin_sqrt(tinv * 3.141592653589793);
                    f1 = tinv * .5 * f0;
                }
#endif
                f1 = f1 * pval * pqpinv;
                cbatch[0] += (qxval - qxsub) * f0 + pqx * f1;
                cbatch[1] += (qyval - qysub) * f0 + pqy * f1;
                cbatch[2] += (qzval - qzsub) * f0 + pqz * f1;
            }
        }
    } else {
        // 1     5   |  (AB|CD)  4-center   spss and psss
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
            const double xp = pxval - pxsub;
            const double yp = pyval - pysub;
            const double zp = pzval - pzsub;
            const double pscale = scalep[ij];
            for (int kl = 0; kl < nkl; kl += 1) {
                const double qval = q[kl];
                const double pqmult = pval * qval;
                const double pqplus = pval + qval;
                const double pqpinv = 1. / pqplus;
                const double pqx = pxval - qx[kl];
                const double pqy = pyval - qy[kl];
                const double pqz = pzval - qz[kl];
                const double t = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult * pqpinv;
                const double scale = pscale * scaleq[kl] / (pqmult * __builtin_sqrt(pqplus));
#ifdef ERD_TABLE_FREE_BOYS_FUNCTIONS
                const double f0 = scale * boys0 (t);
                double f1 = scale * boys1 (t);
#else
                double f0, f1;
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
                    f0 = scale * f0;
                    f1 = scale * f1;
                } else {
                    const double tinv = 1. / t;
                    f0 = scale * .5 * __builtin_sqrt(tinv * 3.141592653589793);
                    f1 = tinv * .5 * f0;
                }
#endif
                f1 = -f1 * qval * pqpinv;
                cbatch[0] += xp * f0 + pqx * f1;
                cbatch[1] += yp * f0 + pqy * f1;
                cbatch[2] += zp * f0 + pqz * f1;
            }
        }
    }
}

#pragma offload_attribute(pop)
