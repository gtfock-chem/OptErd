#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <immintrin.h>
#include "erd.h"

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif


#if defined (__MIC__)
#define ERD_2D_COEFF(n) \
    \
__m512d pqx512 = _mm512_sub_pd(pxij512, qxkl512); \
__m512d pqy512 = _mm512_sub_pd(pyij512, qykl512); \
__m512d pqz512 = _mm512_sub_pd(pzij512, qzkl512); \
\
__m512d twopq512 = _mm512_mul_pd(pqpinv512, Op5_512); \
__m512d pscale512 = _mm512_mul_pd(pij512, pqpinv512); \
__m512d qscale512 = _mm512_sub_pd(one512, pscale512); \
__m512d root512 = _mm512_load_pd(&rts[n]); \
__m512d proot512 = _mm512_mul_pd(pscale512, root512); \
__m512d qroot512 = _mm512_mul_pd(qscale512, root512); \
\
__m512d b00_512 = _mm512_mul_pd(root512, twopq512); \
_mm512_store_pd(&b00[n], b00_512); \
\
__m512d b01_512 = _mm512_sub_pd(one512, proot512); \
b01_512 = _mm512_mul_pd(b01_512, twoq512); \
_mm512_store_pd(&b01[n], b01_512); \
\
__m512d b10_512 = _mm512_sub_pd(one512, qroot512); \
b10_512 = _mm512_mul_pd(b10_512, twop512); \
_mm512_store_pd(&b10[n], b10_512); \
\
__m512d qroot_pqx512 = _mm512_mul_pd(qroot512, pqx512); \
__m512d c00x512 = _mm512_sub_pd(pxij512, qroot_pqx512); \
c00x512 = _mm512_sub_pd(c00x512, xa512); \
_mm512_store_pd(&c00x[n], c00x512); \
\
__m512d qroot_pqy512 = _mm512_mul_pd(qroot512, pqy512); \
__m512d c00y512 = _mm512_sub_pd(pyij512, qroot_pqy512); \
c00y512 = _mm512_sub_pd(c00y512, ya512); \
_mm512_store_pd(&c00y[n], c00y512); \
\
__m512d qroot_pqz512 = _mm512_mul_pd(qroot512, pqz512); \
__m512d c00z512 = _mm512_sub_pd(pzij512, qroot_pqz512); \
c00z512 = _mm512_sub_pd(c00z512, za512); \
_mm512_store_pd(&c00z[n], c00z512); \
\
__m512d proot_pqx512 = _mm512_mul_pd(proot512, pqx512); \
__m512d d00x512 = _mm512_add_pd(qxkl512, proot_pqx512); \
d00x512 = _mm512_sub_pd(d00x512, xc512); \
_mm512_store_pd(&d00x[n], d00x512); \
\
__m512d proot_pqy512 = _mm512_mul_pd(proot512, pqy512); \
__m512d d00y512 = _mm512_add_pd(qykl512, proot_pqy512); \
d00y512 = _mm512_sub_pd(d00y512, yc512); \
_mm512_store_pd(&d00y[n], d00y512); \
\
__m512d proot_pqz512 = _mm512_mul_pd(proot512, pqz512); \
__m512d d00z512 = _mm512_add_pd(qzkl512, proot_pqz512); \
d00z512 = _mm512_sub_pd(d00z512, zc512); \
_mm512_store_pd(&d00z[n], d00z512);

#elif defined (__AVX__)

#define ERD_2D_COEFF(n) \
    __m256d pqx256 = _mm256_sub_pd(pxij256, qxkl256); \
__m256d pqy256 = _mm256_sub_pd(pyij256, qykl256); \
__m256d pqz256 = _mm256_sub_pd(pzij256, qzkl256); \
\
__m256d twopq256 = _mm256_mul_pd(pqpinv256, Op5_256); \
__m256d pscale256 = _mm256_mul_pd(pij256, pqpinv256); \
__m256d qscale256 = _mm256_sub_pd(one256, pscale256); \
__m256d root256 = _mm256_load_pd(&rts[n]); \
__m256d proot256 = _mm256_mul_pd(pscale256, root256); \
__m256d qroot256 = _mm256_mul_pd(qscale256, root256); \
\
__m256d b00_256 = _mm256_mul_pd(root256, twopq256); \
_mm256_store_pd(&b00[n], b00_256); \
\
__m256d b01_256 = _mm256_sub_pd(one256, proot256); \
b01_256 = _mm256_mul_pd(b01_256, twoq256); \
_mm256_store_pd(&b01[n], b01_256); \
\
__m256d b10_256 = _mm256_sub_pd(one256, qroot256); \
b10_256 = _mm256_mul_pd(b10_256, twop256); \
_mm256_store_pd(&b10[n], b10_256); \
\
__m256d qroot_pqx256 = _mm256_mul_pd(qroot256, pqx256); \
__m256d c00x256 = _mm256_sub_pd(pxij256, qroot_pqx256); \
c00x256 = _mm256_sub_pd(c00x256, xa256); \
_mm256_store_pd(&c00x[n], c00x256); \
\
__m256d qroot_pqy256 = _mm256_mul_pd(qroot256, pqy256); \
__m256d c00y256 = _mm256_sub_pd(pyij256, qroot_pqy256); \
c00y256 = _mm256_sub_pd(c00y256, ya256); \
_mm256_store_pd(&c00y[n], c00y256); \
\
__m256d qroot_pqz256 = _mm256_mul_pd(qroot256, pqz256); \
__m256d c00z256 = _mm256_sub_pd(pzij256, qroot_pqz256); \
c00z256 = _mm256_sub_pd(c00z256, za256); \
_mm256_store_pd(&c00z[n], c00z256); \
\
__m256d proot_pqx256 = _mm256_mul_pd(proot256, pqx256); \
__m256d d00x256 = _mm256_add_pd(qxkl256, proot_pqx256); \
d00x256 = _mm256_sub_pd(d00x256, xc256); \
_mm256_store_pd(&d00x[n], d00x256); \
\
__m256d proot_pqy256 = _mm256_mul_pd(proot256, pqy256); \
__m256d d00y256 = _mm256_add_pd(qykl256, proot_pqy256); \
d00y256 = _mm256_sub_pd(d00y256, yc256); \
_mm256_store_pd(&d00y[n], d00y256); \
\
__m256d proot_pqz256 = _mm256_mul_pd(proot256, pqz256); \
__m256d d00z256 = _mm256_add_pd(qzkl256, proot_pqz256); \
d00z256 = _mm256_sub_pd(d00z256, zc256); \
_mm256_store_pd(&d00z[n], d00z256);

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
        double  xa, double  ya, double  za,
        double  xc, double  yc, double  zc,
        double *pinvhf, double *qinvhf, double *pqpinv,
        double *rts, int case2d,
        double *b00, double *b01, double *b10,
        double *c00x, double *c00y, double *c00z,
        double *d00x, double *d00y, double *d00z)
{
    int m, n, ij, ng, kl;
    goto L9;

#if 0
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
#endif

    /*             ...the case P > p-shell and Q > p-shell. */
L9:
    {
        //printf("mij = %d, mkl = %d, ngqp = %d\n", mij, mkl, ngqp);
        //printf("scalar read = %d, vec read = %d\n", mij * 6 + mij*mkl*6, (mij*mkl*ngqp/4 + 1)*12);
#if defined (__MIC__)
        m = 0;
        n = 0;
        int n1 = 0;

#if 0
        double Op5 = 0.5;
        double one = 1.0;
        __m512d Op5_512 = _mm512_extload_pd(&Op5, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        __m512d one512  = _mm512_extload_pd(&one, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        __m512d xa512 = _mm512_extload_pd(&xa, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        __m512d ya512 = _mm512_extload_pd(&ya, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        __m512d za512 = _mm512_extload_pd(&za, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        __m512d xc512 = _mm512_extload_pd(&xc, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        __m512d yc512 = _mm512_extload_pd(&yc, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        __m512d zc512 = _mm512_extload_pd(&zc, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
#else
        __m512d Op5_512 = _mm512_set1_pd(0.5);
        __m512d one512  = _mm512_set1_pd(1.0);
        __m512d xa512 = _mm512_set1_pd(xa);
        __m512d ya512 = _mm512_set1_pd(ya);
        __m512d za512 = _mm512_set1_pd(za);
        __m512d xc512 = _mm512_set1_pd(xc);
        __m512d yc512 = _mm512_set1_pd(yc);
        __m512d zc512 = _mm512_set1_pd(zc);
#endif
        __mmask8 ij_mask = 0xff;
        __mmask8 kl_mask = 0xff;
        __m512d pij512;
        __m512d pxij512, pyij512, pzij512;
        __m512d qxkl512, qykl512, qzkl512;
        __m512d twop512, twoq512;
        __m512d pqpinv512;

        for (ij = 0; ij < mij; ++ij)
        {
            pij512  = _mm512_mask_extload_pd(pij512 , ij_mask, &p[ij] , _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
            pxij512 = _mm512_mask_extload_pd(pxij512, ij_mask, &px[ij], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
            pyij512 = _mm512_mask_extload_pd(pyij512, ij_mask, &py[ij], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
            pzij512 = _mm512_mask_extload_pd(pzij512, ij_mask, &pz[ij], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

            twop512 = _mm512_mask_extload_pd(twop512, ij_mask, &pinvhf[ij], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

            for (kl = 0; kl < mkl; ++kl)
            {
                qxkl512 = _mm512_mask_extload_pd(qxkl512, kl_mask, &qx[kl], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                qykl512 = _mm512_mask_extload_pd(qykl512, kl_mask, &qy[kl], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                qzkl512 = _mm512_mask_extload_pd(qzkl512, kl_mask, &qz[kl], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

                twoq512   = _mm512_mask_extload_pd(twoq512  , kl_mask, &qinvhf[kl], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                pqpinv512 = _mm512_mask_extload_pd(pqpinv512, kl_mask, &pqpinv[m], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

                for (ng = 0; ng < ngqp; ++ng)
                {
                    n1++;
                    ij_mask = ij_mask << 1;
                    kl_mask = kl_mask << 1;
                    if(n1 == SIMDW)
                    {
                        ERD_2D_COEFF(n)
                            n += SIMDW;
                        n1 = 0;
                        ij_mask = 0xff;
                        pij512  = _mm512_extload_pd(&p[ij] , _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                        pxij512 = _mm512_extload_pd(&px[ij], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                        pyij512 = _mm512_extload_pd(&py[ij], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                        pzij512 = _mm512_extload_pd(&pz[ij], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                        twop512 = _mm512_extload_pd(&pinvhf[ij], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

                        kl_mask = 0xff;
                        qxkl512 = _mm512_extload_pd(&qx[kl], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                        qykl512 = _mm512_extload_pd(&qy[kl], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                        qzkl512 = _mm512_extload_pd(&qz[kl], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                        twoq512 = _mm512_extload_pd(&qinvhf[kl], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                        pqpinv512 = _mm512_extload_pd(&pqpinv[m], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                    }
                } /* ng = 0; ng < ngqp; ++ng */
                ++m;
            } /* kl = 0; kl < mkl; ++kl */
        } /* ij = 0; ij < mij; ++ij */

        if(n1 > 0)
        {
            ERD_2D_COEFF(n)
        }

#elif defined (__AVX__)
        m = 0;
        n = 0;
        int n1 = 0;

#if 0
        __m256d Op5_256 = _mm256_set1_pd(0.5);
        __m256d one256  = _mm256_set1_pd(1.0);

        __m256d xa256 = _mm256_broadcast_sd(&xa);
        __m256d ya256 = _mm256_broadcast_sd(&ya);
        __m256d za256 = _mm256_broadcast_sd(&za);
        __m256d xc256 = _mm256_broadcast_sd(&xc);
        __m256d yc256 = _mm256_broadcast_sd(&yc);
        __m256d zc256 = _mm256_broadcast_sd(&zc);
#else
        __m256d Op5_256 = _mm256_set1_pd(0.5);
        __m256d one256  = _mm256_set1_pd(1.0);

        __m256d xa256 = _mm256_set1_pd(xa);
        __m256d ya256 = _mm256_set1_pd(ya);
        __m256d za256 = _mm256_set1_pd(za);
        __m256d xc256 = _mm256_set1_pd(xc);
        __m256d yc256 = _mm256_set1_pd(yc);
        __m256d zc256 = _mm256_set1_pd(zc);
#endif

        __m256d pij256;
        __m256d pxij256;
        __m256d pyij256;
        __m256d pzij256;
        __m256d qxkl256;
        __m256d qykl256;
        __m256d qzkl256;

        __m256d twop256;
        __m256d twoq256;
        __m256d pqpinv256;
        __m256d v0123 = _mm256_set_pd(3.0, 2.0, 1.0, 0.0);

        for (ij = 0; ij < mij; ++ij)
        {
            double n1d = n1;
            __m256d n1_256  = _mm256_broadcast_sd(&n1d);
            __m256d mask_ij = _mm256_cmp_pd(v0123, n1_256, _CMP_GE_OQ);
            __m256d temp256 = _mm256_broadcast_sd(&p[ij]);
            pij256  = _mm256_blendv_pd(pij256, temp256, mask_ij);
            temp256 = _mm256_broadcast_sd(&px[ij]);
            pxij256  = _mm256_blendv_pd(pxij256, temp256, mask_ij);
            temp256 = _mm256_broadcast_sd(&py[ij]);
            pyij256  = _mm256_blendv_pd(pyij256, temp256, mask_ij);
            temp256 = _mm256_broadcast_sd(&pz[ij]);
            pzij256  = _mm256_blendv_pd(pzij256, temp256, mask_ij);
            temp256 = _mm256_broadcast_sd(&pinvhf[ij]);
            twop256 = _mm256_blendv_pd(twop256, temp256, mask_ij);

            for (kl = 0; kl < mkl; ++kl)
            {
                n1d = n1;
                n1_256  = _mm256_broadcast_sd(&n1d);
                __m256d mask_kl = _mm256_cmp_pd(v0123, n1_256, _CMP_GE_OQ);
                temp256 = _mm256_broadcast_sd(&qx[kl]);
                qxkl256 = _mm256_blendv_pd(qxkl256, temp256, mask_kl);
                temp256 = _mm256_broadcast_sd(&qy[kl]);
                qykl256 = _mm256_blendv_pd(qykl256, temp256, mask_kl);
                temp256 = _mm256_broadcast_sd(&qz[kl]);
                qzkl256 = _mm256_blendv_pd(qzkl256, temp256, mask_kl);
                temp256 = _mm256_broadcast_sd(&qinvhf[kl]);
                twoq256 = _mm256_blendv_pd(twoq256, temp256, mask_kl);
                temp256 = _mm256_broadcast_sd(&pqpinv[m]);
                pqpinv256 = _mm256_blendv_pd(pqpinv256, temp256, mask_kl);

                for (ng = 0; ng < ngqp; ++ng)
                {
                    n1++;
                    if(n1 == SIMDW)
                    {
                        ERD_2D_COEFF(n)
                            n += SIMDW;
                        n1 = 0;
                        pij256 = _mm256_broadcast_sd(&p[ij]);
                        pxij256 = _mm256_broadcast_sd(&px[ij]);
                        pyij256 = _mm256_broadcast_sd(&py[ij]);
                        pzij256 = _mm256_broadcast_sd(&pz[ij]);
                        twop256 = _mm256_broadcast_sd(&pinvhf[ij]);
                        qxkl256 = _mm256_broadcast_sd(&qx[kl]);
                        qykl256 = _mm256_broadcast_sd(&qy[kl]);
                        qzkl256 = _mm256_broadcast_sd(&qz[kl]);
                        twoq256 = _mm256_broadcast_sd(&qinvhf[kl]);
                        pqpinv256 = _mm256_broadcast_sd(&pqpinv[m]);
                    }
                } /* ng = 0; ng < ngqp; ++ng */
                ++m;
            } /* kl = 0; kl < mkl; ++kl */
        } /* ij = 0; ij < mij; ++ij */

        if(n1 > 0)
        {
            ERD_2D_COEFF(n)
        }

#else
        double pij, pqx, pqy, pqz, pxij, pyij, pzij, qxkl, qykl, qzkl, root, twop, twoq,
               paxij, payij, pazij, qcxkl, qcykl, qczkl, proot, qroot, twopq,
               pscale, qscale;
        m = 0;
        n = 0;
        for (ij = 0; ij < mij; ++ij)
        {
            pij = p[ij];
            pxij = px[ij];
            pyij = py[ij];
            pzij = pz[ij];
            paxij = pxij - xa;
            payij = pyij - ya;
            pazij = pzij - za;
            twop = pinvhf[ij];
            for (kl = 0; kl < mkl; ++kl)
            {
                qxkl = qx[kl];
                qykl = qy[kl];
                qzkl = qz[kl];
                qcxkl = qxkl - xc;
                qcykl = qykl - yc;
                qczkl = qzkl - zc;
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

#endif
    }

    return 0;
}


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
