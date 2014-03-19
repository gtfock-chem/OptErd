#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"
#include "erdutil.h"

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
ERD_OFFLOAD void erd__2d_coefficients(uint32_t mij, uint32_t mkl, uint32_t ngqp,
    const double *restrict p, const double *restrict q,
    const double *restrict px, const double *restrict py, const double *restrict pz,
    const double *restrict qx, const double *restrict qy, const double *restrict qz,
    const double *restrict pax, const double *restrict pay, const double *restrict paz,
    const double *restrict qcx, const double *restrict qcy, const double *restrict qcz,
    const double *restrict pinvhf, const double *restrict qinvhf, const double *restrict pqpinv,
    const double *restrict rts, uint32_t case2d,
    double *restrict b00, double *restrict b01, double *restrict b10,
    double *restrict c00x, double *restrict c00y, double *restrict c00z,
    double *restrict d00x, double *restrict d00y, double *restrict d00z)
{
#if 0
    switch (case2d) {
        case 1:
            /* ...the case P = s-shell and Q = s-shell. (no coefficients here) */
            return;
        case 2:
        {
            /* ...the case P = p-shell and Q = s-shell. (no B00,B01,B10,D00) */
            uint32_t m = 0;
            uint32_t n = 0;
            for (uint32_t ij = 0; ij < mij; ij++) {
                const double pxij = px[ij];
                const double pyij = py[ij];
                const double pzij = pz[ij];
                const double paxij = pax[ij];
                const double payij = pay[ij];
                const double pazij = paz[ij];
                for (uint32_t kl = 0; kl < mkl; kl++) {
                    const double pqx = pxij - qx[kl];
                    const double pqy = pyij - qy[kl];
                    const double pqz = pzij - qz[kl];           
                    const double qscale = q[kl] * pqpinv[m];
                    ++m;
                    for (uint32_t ng = 0; ng < ngqp; ng++) {   
                        const double qroot = qscale * rts[n];
                        c00x[n] = paxij - qroot * pqx;
                        c00y[n] = payij - qroot * pqy;
                        c00z[n] = pazij - qroot * pqz;
                        ++n;
                    }
                }
            }
            return;
        }
        case 3:
        {
            /* ...the case P > p-shell and Q = s-shell. (no B00,B01,D00) */
            uint32_t m = 0;
            uint32_t n = 0;
            for (uint32_t ij = 0; ij < mij; ++ij) {
                const double pij = p[ij];
                const double pxij = px[ij];
                const double pyij = py[ij];
                const double pzij = pz[ij];
                const double paxij = pax[ij];
                const double payij = pay[ij];
                const double pazij = paz[ij];
                const double twop = pinvhf[ij];
                for (uint32_t kl = 0; kl < mkl; ++kl) {
                    const double pqx = pxij - qx[kl];
                    const double pqy = pyij - qy[kl];
                    const double pqz = pzij - qz[kl];
                    const double qscale = 1.0 - pij * pqpinv[m];
                    ++m;               
                    for (uint32_t ng = 0; ng < ngqp; ++ng) {
                        const double qroot = qscale * rts[n];
                        b10[n] = (1.0 - qroot) * twop;
                        c00x[n] = paxij - qroot * pqx;
                        c00y[n] = payij - qroot * pqy;
                        c00z[n] = pazij - qroot * pqz;
                        ++n;
                    }
                }
            }
            return;
        }
        case 4:
        {
            /* ...the case P = s-shell and Q = p-shell. (no B00,B01,B10,C00) */
            uint32_t m = 0;
            uint32_t n = 0;
            for (uint32_t ij = 0; ij < mij; ++ij) {
                const double pij = p[ij];
                const double pxij = px[ij];
                const double pyij = py[ij];
                const double pzij = pz[ij];
                for (uint32_t kl = 0; kl < mkl; ++kl) {
                    const double qcxkl = qcx[kl];
                    const double qcykl = qcy[kl];
                    const double qczkl = qcz[kl];
                    const double pqx = pxij - qx[kl];
                    const double pqy = pyij - qy[kl];
                    const double pqz = pzij - qz[kl];
                    const double pscale = pij * pqpinv[m];
                    ++m;
                    for (uint32_t ng = 0; ng < ngqp; ++ng) {
                        const double proot = pscale * rts[n];
                        d00x[n] = qcxkl + proot * pqx;
                        d00y[n] = qcykl + proot * pqy;
                        d00z[n] = qczkl + proot * pqz;
                        ++n;
                    }
                }
            }
            return;
        }
        case 5:
        {
            /* ...the case P = p-shell and Q = p-shell. (no B01,B10) */
            uint32_t m = 0;
            uint32_t n = 0;
            for (uint32_t ij = 0; ij < mij; ++ij) {
                const double pij = p[ij];
                const double pxij = px[ij];
                const double pyij = py[ij];
                const double pzij = pz[ij];
                const double paxij = pax[ij];
                const double payij = pay[ij];
                const double pazij = paz[ij];
                for (uint32_t kl = 0; kl < mkl; ++kl) {
                    const double qcxkl = qcx[kl];
                    const double qcykl = qcy[kl];
                    const double qczkl = qcz[kl];
                    const double pqx = pxij - qx[kl];
                    const double pqy = pyij - qy[kl];
                    const double pqz = pzij - qz[kl];
                    const double twopq = pqpinv[m] * 0.5;
                    const double pscale = pij * pqpinv[m];
                    const double qscale = 1.0 - pscale;
                    ++m;
                    for (uint32_t ng = 0; ng < ngqp; ++ng) {
                        const double root = rts[n];
                        const double proot = pscale * root;
                        const double qroot = qscale * root;
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
            return;
        }
        case 6:
        {
            /* ...the case P > p-shell and Q = p-shell. (no B01) */
            uint32_t m = 0;
            uint32_t n = 0;
            for (uint32_t ij = 0; ij < mij; ++ij) {
                const double pij = p[ij];
                const double pxij = px[ij];
                const double pyij = py[ij];
                const double pzij = pz[ij];
                const double paxij = pax[ij];
                const double payij = pay[ij];
                const double pazij = paz[ij];
                const double twop = pinvhf[ij];
                for (uint32_t kl = 0; kl < mkl; ++kl) {
                    const double qcxkl = qcx[kl];
                    const double qcykl = qcy[kl];
                    const double qczkl = qcz[kl];
                    const double pqx = pxij - qx[kl];
                    const double pqy = pyij - qy[kl];
                    const double pqz = pzij - qz[kl];
                    const double twopq = pqpinv[m] * 0.5;
                    const double pscale = pij * pqpinv[m];
                    const double qscale = 1.0 - pscale;
                    ++m;
                    for (uint32_t ng = 0; ng < ngqp; ++ng) {
                        const double root = rts[n];
                        const double proot = pscale * root;
                        const double qroot = qscale * root;
                        b00[n] = root * twopq;
                        b10[n] = (1.0 - qroot) * twop;
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
            return;
        }
        case 7:
        {
            /* ...the case P = s-shell and Q > p-shell. (no B00,B10,C00) */
            uint32_t m = 0;
            uint32_t n = 0;
            for (uint32_t ij = 0; ij < mij; ++ij) {
                const double pij = p[ij];
                const double pxij = px[ij];
                const double pyij = py[ij];
                const double pzij = pz[ij];
                for (uint32_t kl = 0; kl < mkl; ++kl) {
                    const double qcxkl = qcx[kl];
                    const double qcykl = qcy[kl];
                    const double qczkl = qcz[kl];
                    const double twoq = qinvhf[kl];
                    const double pqx = pxij - qx[kl];
                    const double pqy = pyij - qy[kl];
                    const double pqz = pzij - qz[kl];
                    const double pscale = pij * pqpinv[m];
                    ++m;
                    for (uint32_t ng = 0; ng < ngqp; ++ng) {
                        const double proot = pscale * rts[n];
                        b01[n] = (1.0 - proot) * twoq;
                        d00x[n] = qcxkl + proot * pqx;
                        d00y[n] = qcykl + proot * pqy;
                        d00z[n] = qczkl + proot * pqz;
                        ++n;
                    }
                }
            }
            return;
        }
        case 8:
        {
            /* ...the case P = p-shell and Q > p-shell. (no B10) */
            uint32_t m = 0;
            uint32_t n = 0;
            for (uint32_t ij = 0; ij < mij; ++ij) {
                const double pij = p[ij];
                const double pxij = px[ij];
                const double pyij = py[ij];
                const double pzij = pz[ij];
                const double paxij = pax[ij];
                const double payij = pay[ij];
                const double pazij = paz[ij];
                for (uint32_t kl = 0; kl < mkl; ++kl) {
                    const double qcxkl = qcx[kl];
                    const double qcykl = qcy[kl];
                    const double qczkl = qcz[kl];
                    const double twoq = qinvhf[kl];
                    const double pqx = pxij - qx[kl];
                    const double pqy = pyij - qy[kl];
                    const double pqz = pzij - qz[kl];
                    const double twopq = pqpinv[m] * 0.5;
                    const double pscale = pij * pqpinv[m];
                    const double qscale = 1.0 - pscale;
                    ++m;
                    for (uint32_t ng = 0; ng < ngqp; ++ng) {
                        const double root = rts[n];
                        const double proot = pscale * root;
                        const double qroot = qscale * root;
                        b00[n] = root * twopq;
                        b01[n] = (1.0 - proot) * twoq;
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
            return;
        }
        case 9:
        {
            /* ...the case P > p-shell and Q > p-shell. */
            uint32_t m = 0;
            uint32_t n = 0;
            for (uint32_t ij = 0; ij < mij; ++ij) {
                const double pij = p[ij];
                const double pxij = px[ij];
                const double pyij = py[ij];
                const double pzij = pz[ij];
                const double paxij = pax[ij];
                const double payij = pay[ij];
                const double pazij = paz[ij];
                const double twop = pinvhf[ij];
                for (uint32_t kl = 0; kl < mkl; ++kl) {
                    const double qcxkl = qcx[kl];
                    const double qcykl = qcy[kl];
                    const double qczkl = qcz[kl];
                    const double twoq = qinvhf[kl];
                    const double pqx = pxij - qx[kl];
                    const double pqy = pyij - qy[kl];
                    const double pqz = pzij - qz[kl];
                    const double twopq = pqpinv[m] * 0.5;
                    const double pscale = pij * pqpinv[m];
                    const double qscale = 1.0 - pscale;
                    ++m;
                    for (uint32_t ng = 0; ng < ngqp; ++ng) {
                        const double root = rts[n];
                        const double proot = pscale * root;
                        const double qroot = qscale * root;
                        b00[n] = root * twopq;
                        b01[n] = (1.0 - proot) * twoq;
                        b10[n] = (1.0 - qroot) * twop;
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
            return;
        }
    }
#else
    /* General case only */
    uint32_t m = 0;
    uint32_t n = 0;
    for (uint32_t ij = 0; ij < mij; ++ij) {
        const double pij = p[ij];
        const double pxij = px[ij];
        const double pyij = py[ij];
        const double pzij = pz[ij];
        const double paxij = pax[ij];
        const double payij = pay[ij];
        const double pazij = paz[ij];
        const double twop = pinvhf[ij];
        for (uint32_t kl = 0; kl < mkl; ++kl) {
            const double qcxkl = qcx[kl];
            const double qcykl = qcy[kl];
            const double qczkl = qcz[kl];
            const double twoq = qinvhf[kl];
            const double pqx = pxij - qx[kl];
            const double pqy = pyij - qy[kl];
            const double pqz = pzij - qz[kl];
            const double twopq = pqpinv[m] * 0.5;
            const double pscale = pij * pqpinv[m];
            const double qscale = 1.0 - pscale;
            ++m;
            for (uint32_t ng = 0; ng < ngqp; ++ng) {
                const double root = rts[n];
                const double proot = pscale * root;
                const double qroot = qscale * root;
                b00[n] = root * twopq;
                b01[n] = (1.0 - proot) * twoq;
                b10[n] = (1.0 - qroot) * twop;
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
#endif
}
