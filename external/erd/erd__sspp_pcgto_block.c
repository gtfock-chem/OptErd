#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__SSPP_PCGTO_BLOCK */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation is designed to provide ultrafast block */
/*                evaluation of a batch of normalized electron repulsion */
/*                integrals between s-shell and p-shell primitive */
/*                spherical gaussian type orbitals. */
/*                A batch is defined here as containing all possible */
/*                integrals, that is its dimension is determined by */
/*                the total number of primitive functions (here = 9) */
/*                times the total number of ij and kl exponent pair */
/*                combinations. */
/*                The integrals are ordered in the batch the following */
/*                way (first index varying fastest): */
/*                    batch (nxyz1,nxyz2,nxyz3,nxyz4,kl,ij) */
/*                where ij and kl indicates alpha exponent pairs */
/*                defining the present block. */
/*                The present routine evaluates batches of the type: */
/*                      sspp , spsp , pssp , spps , psps , ppss */
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
/*                                    sspp/spsp/pssp/spps/psps/ppss */
/*                                    integral batch */
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
/*                                    cartesian sspp/spsp/pssp/spps/ */
/*                                    psps/ppss integrals */
/* ------------------------------------------------------------------------ */
int erd__sspp_pcgto_block (int nbatch, int atomic, int atom12, int atom34,
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

    int i;
    int j;
    int k;
    int l;
    int m;
    double t;
    double f0;
    double f1;
    double f2;
    int ij;
    int kl;
    double u0;
    double u1;   
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
    int tcase;
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
    double rnpqsq;
    double xps1;
    double yps1;
    double zps1;
    double xsp1;
    double ysp1;
    double zsp1;
    double xsp2;
    double ysp2;
    double zsp2;
    double xps2;
    double yps2;
    double zps2;
    double a;
    double b;
    double c;
    double d;
    double e;
    double f;
    double g;

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
    

/*             ...jump according to where the s-type functions are */
/*                located and the type of 'K4' loop: */
/*                SHELLP CASE  |        Integral center type */
/*               --------------|------------------------------------ */
/*                   0     1   |  (AA|AA)  atomic     sspp */
/*                   0     2   |  (AA|CC)  2-center   sspp */
/*                   0     3   |  (AB|CC)  3-center   sspp */
/*                   0     4   |  (AA|CD)  3-center   sspp */
/*                   0     5   |  (AB|CD)  4-center   sspp */
/*                   1     1   |  (AA|AA)  atomic     spsp,spps,pssp,psps */
/*                   1     2   |  (AA|CC)  2-center   spsp,spps,pssp,psps */
/*                   1     3   |  (AB|CC)  3-center   spsp,spps,pssp,psps */
/*                   1     4   |  (AA|CD)  3-center   spsp,spps,pssp,psps */
/*                   1     5   |  (AB|CD)  4-center   spsp,spps,pssp,psps */
/*                   2     1   |  (AA|AA)  atomic     ppss */
/*                   2     2   |  (AA|CC)  2-center   ppss */
/*                   2     3   |  (AB|CC)  3-center   ppss */
/*                   2     4   |  (AA|CD)  3-center   ppss */
/*                   2     5   |  (AB|CD)  4-center   ppss */


    switch (shellp * 5 + tcase)
    {
    case 5:
        goto L15;
    case 10:
        goto L25;
    case 15:
        goto L35;
    }

  L15:
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
                f0 = scale * f0;
                f1 = scale * f1;
                f2 = scale * f2;
            }
            else
            {
                tinv = 1. / t;
                t2inv = tinv * .5;
                f0 = scale * .5 * sqrt (tinv * 3.141592653589793);
                f1 = t2inv * f0;
                f2 = t2inv * 3. * f1;
            }
            u0 = pval * pqpinv;
            u1 = (f0 - u0 * f1) / (qval + qval);
            xps1 = qxval - x3;
            yps1 = qyval - y3;
            zps1 = qzval - z3;
            xsp1 = qxval - x4;
            ysp1 = qyval - y4;
            zsp1 = qzval - z4;
            xsp2 = pqx * u0;
            ysp2 = pqy * u0;
            zsp2 = pqz * u0;
            a = xps1 * f0 + xsp2 * f1;
            b = xps1 * f1 + xsp2 * f2;
            c = yps1 * f0 + ysp2 * f1;
            d = yps1 * f1 + ysp2 * f2;
            e = zps1 * f0 + zsp2 * f1;
            f = zps1 * f1 + zsp2 * f2;
            batch[m + 1] = xsp1 * a + xsp2 * b + u1;
            batch[m + 2] = xsp1 * c + xsp2 * d;
            batch[m + 3] = xsp1 * e + xsp2 * f;
            batch[m + 4] = ysp1 * a + ysp2 * b;
            batch[m + 5] = ysp1 * c + ysp2 * d + u1;
            batch[m + 6] = ysp1 * e + ysp2 * f;
            batch[m + 7] = zsp1 * a + zsp2 * b;
            batch[m + 8] = zsp1 * c + zsp2 * d;
            batch[m + 9] = zsp1 * e + zsp2 * f + u1;
            m += 9;
        }
    }
    return 0;

  L25:
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
        pxval = px[ij];
        pyval = py[ij];
        pzval = pz[ij];
        pscale = scalep[ij];
        xps1 = pxval - pxsub;
        yps1 = pyval - pysub;
        zps1 = pzval - pzsub;
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
                f0 = scale * f0;
                f1 = scale * f1;
                f2 = scale * f2;
            }
            else
            {
                tinv = 1. / t;
                t2inv = tinv * .5;
                f0 = scale * .5 * sqrt (tinv * 3.141592653589793);
                f1 = t2inv * f0;
                f2 = t2inv * 3. * f1;
            }
            u0 = pval * pqpinv;
            u1 = -qval * pqpinv;
            xsp1 = qxval - qxsub;
            ysp1 = qyval - qysub;
            zsp1 = qzval - qzsub;
            xsp2 = pqx * u0;
            ysp2 = pqy * u0;
            zsp2 = pqz * u0;
            xps2 = pqx * u1;
            yps2 = pqy * u1;
            zps2 = pqz * u1;
            a = xps1 * f0 + xps2 * f1;
            b = xps1 * f1 + xps2 * f2;
            c = yps1 * f0 + yps2 * f1;
            d = yps1 * f1 + yps2 * f2;
            e = zps1 * f0 + zps2 * f1;
            f = zps1 * f1 + zps2 * f2;
            g = f1 * .5 * pqpinv;
            batch[m + 1] = xsp1 * a + xsp2 * b + g;
            batch[m + 2] = xsp1 * c + xsp2 * d;
            batch[m + 3] = xsp1 * e + xsp2 * f;
            batch[m + 4] = ysp1 * a + ysp2 * b;
            batch[m + 5] = ysp1 * c + ysp2 * d + g;
            batch[m + 6] = ysp1 * e + ysp2 * f;
            batch[m + 7] = zsp1 * a + zsp2 * b;
            batch[m + 8] = zsp1 * c + zsp2 * d;
            batch[m + 9] = zsp1 * e + zsp2 * f + g;
            m += 9;
        }
    }
    return 0;

  L35:
    m = 0;
    for (ij = 1; ij <= mij; ++ij)
    {
        pval = p[ij];
        pxval = px[ij];
        pyval = py[ij];
        pzval = pz[ij];
        pscale = scalep[ij];
        xps1 = pxval - x1;
        yps1 = pyval - y1;
        zps1 = pzval - z1;
        xsp1 = pxval - x2;
        ysp1 = pyval - y2;
        zsp1 = pzval - z2;
        u1 = 1. / (pval + pval);
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
                f0 = scale * f0;
                f1 = scale * f1;
                f2 = scale * f2;
            }
            else
            {
                tinv = 1. / t;
                t2inv = tinv * .5;
                f0 = scale * .5 * sqrt (tinv * 3.141592653589793);
                f1 = t2inv * f0;
                f2 = t2inv * 3. * f1;
            }
            u0 = -qval * pqpinv;
            xsp2 = pqx * u0;
            ysp2 = pqy * u0;
            zsp2 = pqz * u0;
            a = xps1 * f0 + xsp2 * f1;
            b = xps1 * f1 + xsp2 * f2;
            c = yps1 * f0 + ysp2 * f1;
            d = yps1 * f1 + ysp2 * f2;
            e = zps1 * f0 + zsp2 * f1;
            f = zps1 * f1 + zsp2 * f2;
            g = u1 * (f0 + u0 * f1);
            batch[m + 1] = xsp1 * a + xsp2 * b + g;
            batch[m + 2] = xsp1 * c + xsp2 * d;
            batch[m + 3] = xsp1 * e + xsp2 * f;
            batch[m + 4] = ysp1 * a + ysp2 * b;
            batch[m + 5] = ysp1 * c + ysp2 * d + g;
            batch[m + 6] = ysp1 * e + ysp2 * f;
            batch[m + 7] = zsp1 * a + zsp2 * b;
            batch[m + 8] = zsp1 * c + zsp2 * d;
            batch[m + 9] = zsp1 * e + zsp2 * f + g;
            m += 9;
        }
    }

    return 0;
}