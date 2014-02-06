#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__2D_COEFFICIENTS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation evaluates the VRR coefficients for */
/*                the 2D integrals for the present set of NGQP roots */
/*                corresponding to all i,j,k,l exponent quadruplets. */
/*                Not all coefficients are needed in case there are */
/*                s- or p-shells present on either P or Q side. Also */
/*                observe that the cartesian distance components PAX, */
/*                PAY and PAZ are all equal to zero in case the two */
/*                atomic centers A and B coincide, in which case they */
/*                need not to be addressed inside the algorithm. */
/*                Likewise for QCX,QCY and QCZ, if centers C and D */
/*                coincide. */
/*                  Input: */
/*                    MIJ(KL)      =  current # of ij (kl) primitive */
/*                                    index pairs corresponding to */
/*                                    the csh pairs A,B (C,D) */
/*                    MIJKL        =  current # of ijkl primitive */
/*                                    index quadruplets (= MIJ*MKL) */
/*                    NGQP         =  # of gaussian quadrature points */
/*                                    (roots) */
/*                    MGQIJKL      =  # of roots times # of ijkl */
/*                                    quadruplets (= NGQP*MIJKL) */
/*                    ATOMAB(CD)   =  is true, if centers A and B */
/*                                    (C and D) coincide */
/*                    P(Q)         =  current MIJ (MKL) exponent sums */
/*                                    for csh A and B (C and D) */
/*                    Px(Qx)       =  current MIJ (MKL) coordinates */
/*                                    x=X,Y,Z for gaussian product */
/*                                    centers P=A+B (Q=C+D) */
/*                    PAx(QCx)     =  current MIJ (MKL) coordinate */
/*                                    x=X,Y,Z differences P-A (Q-C) */
/*                                    between centers P and A (Q and C) */
/*                    P(Q)INVHF    =  current MIJ (MKL) values of */
/*                                    1/(2*P(Q)), where P and Q are */
/*                                    the corresponding exponent sums */
/*                                    for csh A and B (C and D) */
/*                    PQPINV       =  current MIJKL values of 1/(P+Q) */
/*                    RTS          =  current MGQIJKL values of all */
/*                                    quadrature roots */
/*                    CASE2D       =  int value within the range */
/*                                    from 1 to 9, indicating which */
/*                                    2D coefficient evaluation case */
/*                                    is present to trigger specific */
/*                                    simplified sections of the code */
/*                  Output: */
/*                    Bxx          =  the coordinate independent */
/*                                    B-coefficients (xx=00,01,10) */
/*                    C00x         =  the C-coefficients (individual */
/*                                    cartesian components x=X,Y,Z) for */
/*                                    shell expansion on center P */
/*                    D00x         =  the D-coefficients (individual */
/*                                    cartesian components x=X,Y,Z) for */
/*                                    shell expansion on center Q */
/* ------------------------------------------------------------------------ */
int erd__2d_coefficients (int mij, int mkl, int ngqp,
                          double *p, double *q,
                          double *px, double *py, double *pz,
                          double *qx, double *qy, double *qz,
                          double *pax, double *pay, double *paz,
                          double *qcx, double *qcy, double *qcz,
                          double *pinvhf, double *qinvhf, double *pqpinv,
                          double *rts, int case2d,
                          double *b00, double *b01, double *b10,
                          double *c00x, double *c00y, double *c00z,
                          double *d00x, double *d00y, double *d00z)
{
    int m, n, ij, ng, kl;
    double pij, pqx, pqy, pqz, pxij, pyij, pzij, root, twop, twoq,
        paxij, payij, pazij, qcxkl, qcykl, qczkl, proot, qroot, twopq,
        pscale, qscale;
    goto L9;
    
    switch (case2d)
    {
    case 1:
        goto L1;
    case 2:
        goto L2;
    case 3:
        goto L5;
    case 4:
        goto L3;
    case 5:
        goto L4;
    case 6:
        goto L7;
    case 7:
        goto L6;
    case 8:
        goto L8;
    case 9:
        goto L9;
    }


/*             ...the case P = s-shell and Q = s-shell. */
/*                (no coefficients here) */
  L1:
    return 0;


/*             ...the case P = p-shell and Q = s-shell. */
/*                (no B00,B01,B10,D00) */
  L2:
    {
        m = 0;
        n = 0;
        for (ij = 0; ij < mij; ++ij)
        {
            pxij = px[ij];
            pyij = py[ij];
            pzij = pz[ij];
            paxij = pax[ij];
            payij = pay[ij];
            pazij = paz[ij];
            for (kl = 0; kl < mkl; ++kl)
            {
                pqx = pxij - qx[kl];
                pqy = pyij - qy[kl];
                pqz = pzij - qz[kl];           
                qscale = q[kl] * pqpinv[m];
                ++m;
                for (ng = 0; ng < ngqp; ++ng)
                {   
                    qroot = qscale * rts[n];
                    c00x[n] = paxij - qroot * pqx;
                    c00y[n] = payij - qroot * pqy;
                    c00z[n] = pazij - qroot * pqz;
                    ++n;
                }
            }
        }
    }
    return 0;


/*             ...the case P = s-shell and Q = p-shell. */
/*                (no B00,B01,B10,C00) */
  L3:
    {
        m = 0;
        n = 0;
        for (ij = 0; ij < mij; ++ij)
        {
            pij = p[ij];
            pxij = px[ij];
            pyij = py[ij];
            pzij = pz[ij];
            for (kl = 0; kl < mkl; ++kl)
            {
                qcxkl = qcx[kl];
                qcykl = qcy[kl];
                qczkl = qcz[kl];
                pqx = pxij - qx[kl];
                pqy = pyij - qy[kl];
                pqz = pzij - qz[kl];
                pscale = pij * pqpinv[m];
                ++m;
                for (ng = 0; ng < ngqp; ++ng)
                {
                    proot = pscale * rts[n];
                    d00x[n] = qcxkl + proot * pqx;
                    d00y[n] = qcykl + proot * pqy;
                    d00z[n] = qczkl + proot * pqz;
                    ++n;
                }
            }
        }
    }
    return 0;


/*             ...the case P = p-shell and Q = p-shell. */
/*                (no B01,B10) */
  L4:
    {
        m = 0;
        n = 0;
        for (ij = 0; ij < mij; ++ij)
        {
            pij = p[ij];
            pxij = px[ij];
            pyij = py[ij];
            pzij = pz[ij];
            paxij = pax[ij];
            payij = pay[ij];
            pazij = paz[ij];
            for (kl = 0; kl < mkl; ++kl)
            {
                qcxkl = qcx[kl];
                qcykl = qcy[kl];
                qczkl = qcz[kl];
                pqx = pxij - qx[kl];
                pqy = pyij - qy[kl];
                pqz = pzij - qz[kl];
                twopq = pqpinv[m] * .5;
                pscale = pij * pqpinv[m];
                qscale = 1. - pscale;
                ++m;
                for (ng = 0; ng < ngqp; ++ng)
                {
                    root = rts[n];
                    proot = pscale * root;
                    qroot = qscale * root;
                    b00[n] = root * twopq;
                    c00x[n] = paxij - qroot * pqx;
                    c00y[n] = payij - qroot * pqy;
                    c00z[n] = pazij - qroot * pqz;
                    d00x[n] = qcxkl + proot * pqx;
                    d00y[n] = qcykl + proot * pqy;
                    d00z[n] = qczkl + proot * pqz;
                    ++n;
                }
            }
        }
    }
    return 0;


/*             ...the case P > p-shell and Q = s-shell. */
/*                (no B00,B01,D00) */
  L5:
    {
        m = 0;
        n = 0;
        for (ij = 0; ij < mij; ++ij)
        {
            pij = p[ij];
            pxij = px[ij];
            pyij = py[ij];
            pzij = pz[ij];
            paxij = pax[ij];
            payij = pay[ij];
            pazij = paz[ij];
            twop = pinvhf[ij];
            for (kl = 0; kl < mkl; ++kl)
            {
                pqx = pxij - qx[kl];
                pqy = pyij - qy[kl];
                pqz = pzij - qz[kl];
                qscale = 1. - pij * pqpinv[m];
                ++m;               
                for (ng = 0; ng < ngqp; ++ng)
                {
                    qroot = qscale * rts[n];
                    b10[n] = (1. - qroot) * twop;
                    c00x[n] = paxij - qroot * pqx;
                    c00y[n] = payij - qroot * pqy;
                    c00z[n] = pazij - qroot * pqz;
                    ++n;
                }
            }
        }
    }
    return 0;


/*             ...the case P = s-shell and Q > p-shell. */
/*                (no B00,B10,C00) */
  L6:
    {
        m = 0;
        n = 0;
        for (ij = 0; ij < mij; ++ij)
        {
            pij = p[ij];
            pxij = px[ij];
            pyij = py[ij];
            pzij = pz[ij];
            for (kl = 0; kl < mkl; ++kl)
            {
                qcxkl = qcx[kl];
                qcykl = qcy[kl];
                qczkl = qcz[kl];
                twoq = qinvhf[kl];
                pqx = pxij - qx[kl];
                pqy = pyij - qy[kl];
                pqz = pzij - qz[kl];
                pscale = pij * pqpinv[m];
                ++m;
                for (ng = 0; ng < ngqp; ++ng)
                {
                    proot = pscale * rts[n];
                    b01[n] = (1. - proot) * twoq;
                    d00x[n] = qcxkl + proot * pqx;
                    d00y[n] = qcykl + proot * pqy;
                    d00z[n] = qczkl + proot * pqz;
                    ++n;
                }
            }
        }
    }
    return 0;


/*             ...the case P > p-shell and Q = p-shell. */
/*                (no B01) */
  L7:
    {
        m = 0;
        n = 0;
        for (ij = 0; ij < mij; ++ij)
        {
            pij = p[ij];
            pxij = px[ij];
            pyij = py[ij];
            pzij = pz[ij];
            paxij = pax[ij];
            payij = pay[ij];
            pazij = paz[ij];
            twop = pinvhf[ij];
            for (kl = 0; kl < mkl; ++kl)
            {
                qcxkl = qcx[kl];
                qcykl = qcy[kl];
                qczkl = qcz[kl];
                pqx = pxij - qx[kl];
                pqy = pyij - qy[kl];
                pqz = pzij - qz[kl];
                twopq = pqpinv[m] * .5;
                pscale = pij * pqpinv[m];
                qscale = 1. - pscale;
                ++m;
                for (ng = 0; ng < ngqp; ++ng)
                {
                    root = rts[n];
                    proot = pscale * root;
                    qroot = qscale * root;
                    b00[n] = root * twopq;
                    b10[n] = (1. - qroot) * twop;
                    c00x[n] = paxij - qroot * pqx;
                    c00y[n] = payij - qroot * pqy;
                    c00z[n] = pazij - qroot * pqz;
                    d00x[n] = qcxkl + proot * pqx;
                    d00y[n] = qcykl + proot * pqy;
                    d00z[n] = qczkl + proot * pqz;
                    ++n;
                }
            }
        }
    }
    return 0;


/*             ...the case P = p-shell and Q > p-shell. */
/*                (no B10) */
  L8:
    {
        m = 0;
        n = 0;
        for (ij = 0; ij < mij; ++ij)
        {
            pij = p[ij];
            pxij = px[ij];
            pyij = py[ij];
            pzij = pz[ij];
            paxij = pax[ij];
            payij = pay[ij];
            pazij = paz[ij];
            for (kl = 0; kl < mkl; ++kl)
            {
                qcxkl = qcx[kl];
                qcykl = qcy[kl];
                qczkl = qcz[kl];
                twoq = qinvhf[kl];
                pqx = pxij - qx[kl];
                pqy = pyij - qy[kl];
                pqz = pzij - qz[kl];
                twopq = pqpinv[m] * .5;
                pscale = pij * pqpinv[m];
                qscale = 1. - pscale;
                ++m;
                for (ng = 0; ng < ngqp; ++ng)
                {
                    root = rts[n];
                    proot = pscale * root;
                    qroot = qscale * root;
                    b00[n] = root * twopq;
                    b01[n] = (1. - proot) * twoq;
                    c00x[n] = paxij - qroot * pqx;
                    c00y[n] = payij - qroot * pqy;
                    c00z[n] = pazij - qroot * pqz;
                    d00x[n] = qcxkl + proot * pqx;
                    d00y[n] = qcykl + proot * pqy;
                    d00z[n] = qczkl + proot * pqz;
                    ++n;
                }
            }
        }
    }
    return 0;


/*             ...the case P > p-shell and Q > p-shell. */
  L9:
    {
        m = 0;
        n = 0;
        for (ij = 0; ij < mij; ++ij)
        {
            pij = p[ij];
            pxij = px[ij];
            pyij = py[ij];
            pzij = pz[ij];
            paxij = pax[ij];
            payij = pay[ij];
            pazij = paz[ij];
            twop = pinvhf[ij];
            for (kl = 0; kl < mkl; ++kl)
            {
                qcxkl = qcx[kl];
                qcykl = qcy[kl];
                qczkl = qcz[kl];
                twoq = qinvhf[kl];
                pqx = pxij - qx[kl];
                pqy = pyij - qy[kl];
                pqz = pzij - qz[kl];
                twopq = pqpinv[m] * .5;
                pscale = pij * pqpinv[m];
                qscale = 1. - pscale;
                ++m;
                for (ng = 0; ng < ngqp; ++ng)
                {
                    root = rts[n];
                    proot = pscale * root;
                    qroot = qscale * root;
                    b00[n] = root * twopq;
                    b01[n] = (1. - proot) * twoq;
                    b10[n] = (1. - qroot) * twop;
                    c00x[n] = paxij - qroot * pqx;
                    c00y[n] = payij - qroot * pqy;
                    c00z[n] = pazij - qroot * pqz;
                    d00x[n] = qcxkl + proot * pqx;
                    d00y[n] = qcykl + proot * pqy;
                    d00z[n] = qczkl + proot * pqz;
                    ++n;
                }
            }
        }
    }


    return 0;
}


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
