#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__SSSS_PCGTO_BLOCK */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation is designed to provide ultrafast block */
/*                evaluation of a batch of normalized electron repulsion */
/*                integrals between s-shell primitive gaussian type */
/*                orbitals. */
/*                A batch is defined here as containing all possible */
/*                integrals, that is its dimension is determined by */
/*                the total number of primitive functions (here = 1) */
/*                times the total number of ij and kl exponent pair */
/*                combinations. */
/*                The integrals are ordered in the batch the following */
/*                way (first index varying fastest): */
/*                    batch (nxyz1,nxyz2,nxyz3,nxyz4,kl,ij) */
/*                where ij and kl indicates alpha exponent pairs */
/*                defining the present block. */
/*                The present routine evaluates batches of the type: */
/*                                     ssss */
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
/*                                    ssss integral batch */
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
/*                                    cartesian ssss integrals */
/* ------------------------------------------------------------------------ */
int erd__ssss_pcgto_block (int nbatch, int atomic, int atom12, int atom34,
                           int mij, int mkl,
                           int nij, int nijbeg, int nijend,
                           int nkl, int nklbeg, int nklend,
                           int npgto1, int npgto2, int npgto3, int npgto4,
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
                           double *p, double *px, double *py, double *pz,
                           double *scalep, double *q, double *qx,
                           double *qy, double *qz,
                           double *scaleq, double *batch)
{
    int ftable_dim1, ftable_offset;

    /* Local variables */
    int i;
    int j;
    int k;
    int l;
    int m;
    double t;
    double f0;
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
    double scale;
    double delta;
    int tgrid;
    double pxval;
    double pyval;
    double pzval;
    double pscale;
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
    ftable_dim1 = mgrid + 1;
    ftable_offset = 0 + ftable_dim1 * 0;
    ftable -= ftable_offset;

    /* Function Body */
    tcase = 5;
    if (atomic || atom12)
    {
        m = 0;
        for (ij = nijbeg; ij <= nijend; ij++)
        {
            m++;
            i = prim1[m];
            j = prim2[m];
            p[m] = alpha1[i] + alpha2[j];
            scalep[m] = norm1[i] * norm2[j];
        }
        tcase--;
    }
    else
    {
        m = 0;
        for (ij = nijbeg; ij <= nijend; ij++)
        {
            m++;
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
    if (atomic || atom34)
    {
        m = 0;
        for (kl = nklbeg; kl <= nklend; ++kl)
        {
            m++;
            k = prim3[m];
            l = prim4[m];
            q[m] = alpha3[k] + alpha4[l];
            scaleq[m] = norm3[k] * norm4[l];
        }
        tcase += -2;
    }
    else
    {
        m = 0;
        for (kl = nklbeg; kl <= nklend; ++kl)
        {
            m++;
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
    if (atomic)
    {
        tcase--;
    }


/*             ...jump according to type of 'K4' loop: */
/*                       TCASE |  Integral center type */
/*                     --------|----------------------- */
/*                         1   |     (AA|AA)  atomic */
/*                         2   |     (AA|CC)  2-center */
/*                         3   |     (AB|CC)  3-center */
/*                         4   |     (AA|CD)  3-center */
/*                         5   |     (AB|CD)  4-center */
    switch (tcase)
    {
    case 1:
        goto L10;
    case 2:
        goto L20;
    case 3:
        goto L30;
    case 4:
        goto L40;
    case 5:
        goto L50;
    }


  L10:
    m = 0;
    for (ij = 1; ij <= mij; ++ij)
    {
        pval = p[ij];
        pscale = scalep[ij];
        for (kl = 1; kl <= mkl; ++kl)
        {
            qval = q[kl];
            ++m;
            batch[m] =
                pscale * scaleq[kl] / (pval * qval * sqrt (pval + qval));
        }
    }
    return 0;

  L20:
    pqx = x1 - x3;
    pqy = y1 - y3;
    pqz = z1 - z3;
    rnpqsq = pqx * pqx + pqy * pqy + pqz * pqz;
    m = 0;
    for (ij = 1; ij <= mij; ++ij)
    {
        pval = p[ij];
        pscale = scalep[ij];
        for (kl = 1; kl <= mkl; ++kl)
        {
            qval = q[kl];
            pqmult = pval * qval;
            pqplus = pval + qval;
            t = rnpqsq * pqmult / pqplus;
            scale = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
            if (t <= tmax)
            {
                tgrid = (int) (t * tvstep + .5);
                delta = tgrid * tstep - t;
                f0 = (((((ftable[tgrid * ftable_dim1 + 6] * delta * 0.166666666666667 +
                          ftable[tgrid * ftable_dim1 + 5]) * delta * 0.2 +
                         ftable[tgrid * ftable_dim1 + 4]) * delta * 0.25 +
                        ftable[tgrid * ftable_dim1 +
                               3]) * delta * 0.333333333333333 +
                       ftable[tgrid * ftable_dim1 + 2]) * delta * 0.5 +
                      ftable[tgrid * ftable_dim1 + 1]) * delta +
                    ftable[tgrid * ftable_dim1];
            }
            else
            {
                f0 = sqrt (M_PI / t) * 0.5;
            }
            ++m;
            batch[m] = scale * f0;
        }
    }
    return 0;

  L30:
    m = 0;
    for (ij = 1; ij <= mij; ++ij)
    {
        pval = p[ij];
        pqx = px[ij] - x3;
        pqy = py[ij] - y3;
        pqz = pz[ij] - z3;
        rnpqsq = pqx * pqx + pqy * pqy + pqz * pqz;
        pscale = scalep[ij];
        for (kl = 1; kl <= mkl; ++kl)
        {
            qval = q[kl];
            pqmult = pval * qval;
            pqplus = pval + qval;
            t = rnpqsq * pqmult / pqplus;
            scale = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
            if (t <= tmax)
            {
                tgrid = (int) (t * tvstep + .5);
                delta = tgrid * tstep - t;
                f0 = (((((ftable[tgrid * ftable_dim1 + 6] * delta *
                          .166666666666667 + ftable[tgrid * ftable_dim1 +
                                                    5]) * delta * .2 +
                         ftable[tgrid * ftable_dim1 + 4]) * delta * .25 +
                        ftable[tgrid * ftable_dim1 +
                               3]) * delta * .333333333333333 +
                       ftable[tgrid * ftable_dim1 + 2]) * delta * .5 +
                      ftable[tgrid * ftable_dim1 + 1]) * delta +
                    ftable[tgrid * ftable_dim1];
            }
            else
            {
                f0 = sqrt (M_PI / t) * .5;
            }
            ++m;
            batch[m] = scale * f0;
        }
    }
    return 0;

  L40:
    m = 0;
    for (ij = 1; ij <= mij; ++ij)
    {
        pval = p[ij];
        pscale = scalep[ij];
        for (kl = 1; kl <= mkl; ++kl)
        {
            qval = q[kl];
            pqmult = pval * qval;
            pqplus = pval + qval;
            pqx = x1 - qx[kl];
            pqy = y1 - qy[kl];
            pqz = z1 - qz[kl];
            t = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult / pqplus;
            scale = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
            if (t <= tmax)
            {
                tgrid = (int) (t * tvstep + .5);
                delta = tgrid * tstep - t;
                f0 = (((((ftable[tgrid * ftable_dim1 + 6] * delta *
                          .166666666666667 + ftable[tgrid * ftable_dim1 +
                                                    5]) * delta * .2 +
                         ftable[tgrid * ftable_dim1 + 4]) * delta * .25 +
                        ftable[tgrid * ftable_dim1 +
                               3]) * delta * .333333333333333 +
                       ftable[tgrid * ftable_dim1 + 2]) * delta * .5 +
                      ftable[tgrid * ftable_dim1 + 1]) * delta +
                    ftable[tgrid * ftable_dim1];
            }
            else
            {
                f0 = sqrt (M_PI / t) * .5;
            }
            ++m;
            batch[m] = scale * f0;
        }
    }
    return 0;

  L50:
    m = 0;
    for (ij = 1; ij <= mij; ++ij)
    {
        pval = p[ij];
        pxval = px[ij];
        pyval = py[ij];
        pzval = pz[ij];
        pscale = scalep[ij];
        for (kl = 1; kl <= mkl; ++kl)
        {
            qval = q[kl];
            pqmult = pval * qval;
            pqplus = pval + qval;
            pqx = pxval - qx[kl];
            pqy = pyval - qy[kl];
            pqz = pzval - qz[kl];
            t = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult / pqplus;
            scale = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
            if (t <= tmax)
            {
                tgrid = (int) (t * tvstep + .5);
                delta = tgrid * tstep - t;
                f0 = (((((ftable[tgrid * ftable_dim1 + 6] * delta *
                          .166666666666667 + ftable[tgrid * ftable_dim1 +
                                                    5]) * delta * .2 +
                         ftable[tgrid * ftable_dim1 + 4]) * delta * .25 +
                        ftable[tgrid * ftable_dim1 +
                               3]) * delta * .333333333333333 +
                       ftable[tgrid * ftable_dim1 + 2]) * delta * .5 +
                      ftable[tgrid * ftable_dim1 + 1]) * delta +
                    ftable[tgrid * ftable_dim1];
            }
            else
            {
                f0 = sqrt (M_PI / t) * .5;
            }
            ++m;
            batch[m] = scale * f0;
        }
    }

    return 0;
}