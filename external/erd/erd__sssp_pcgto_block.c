#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "boys.h"
//#define ERD_TABLE_FREE_BOYS_FUNCTIONS

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
int erd__sssp_pcgto_block (int atomic, int nij, int nkl,
                           int shell1, int shell3, int shellp,
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

    int i;
    int j;
    int k;
    int l;
    int m;
    double t;
    double f0;
    double f1;
    int ij, kl;
    double xp;
    double yp;
    double zp;
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
    double scale;
    int tgrid;
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
    double delta1;
    double delta2;
    double delta3;
    double delta4;
    double delta5;
    double delta6;
    double qzsub;
    double pscale;
    double pqpinv;
    double pqmult;
    double pqplus;
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
    
    if (atomic)
        return 0;
    
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

    // 0     5   |  (AB|CD)  4-center   sssp and ssps
    if (shellp == 0)
    {
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
        for (ij = 0; ij < nij; ++ij)
        {
            pval = p[ij];
            pxval = px[ij];
            pyval = py[ij];
            pzval = pz[ij];
            pscale = scalep[ij];
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
                    f0 = scale * f0;
                    f1 = scale * f1;
                }
                else
                {
                    tinv = 1. / t;
                    f0 = scale * .5 * sqrt (tinv * 3.141592653589793);
                    f1 = tinv * .5 * f0;
                }
#endif
                f1 = f1 * pval * pqpinv;
                batch[m] = (qxval - qxsub) * f0 + pqx * f1;
                batch[m + 1] = (qyval - qysub) * f0 + pqy * f1;
                batch[m + 2] = (qzval - qzsub) * f0 + pqz * f1;
                m += 3;
            }
        }
    }
    // 1     5   |  (AB|CD)  4-center   spss and psss
    else
    {
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
        for (ij = 0; ij < nij; ++ij)
        {
            pval = p[ij];
            pxval = px[ij];
            pyval = py[ij];
            pzval = pz[ij];
            xp = pxval - pxsub;
            yp = pyval - pysub;
            zp = pzval - pzsub;
            pscale = scalep[ij];
            for (kl = 0; kl < nkl; ++kl)
            {
                qval = q[kl];
                pqmult = pval * qval;
                pqplus = pval + qval;
                pqpinv = 1. / pqplus;
                pqx = pxval - qx[kl];
                pqy = pyval - qy[kl];
                pqz = pzval - qz[kl];
                t = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult * pqpinv;
                scale = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
#ifdef ERD_TABLE_FREE_BOYS_FUNCTIONS
                f0 = scale * boys0(t);
                f1 = scale * boys1(t);
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
                    f0 = scale * f0;
                    f1 = scale * f1;
                }
                else
                {
                    tinv = 1. / t;
                    f0 = scale * .5 * sqrt (tinv * 3.141592653589793);
                    f1 = tinv * .5 * f0;
                }
#endif
                f1 = -f1 * qval * pqpinv;
                batch[m] = xp * f0 + pqx * f1;
                batch[m + 1] = yp * f0 + pqy * f1;
                batch[m + 2] = zp * f0 + pqz * f1;
                m += 3;
            }
        }
    }

    return 0;
}