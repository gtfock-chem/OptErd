#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "erd.h"
uint64_t erd__rys_roots_weights_ticks[256];
uint64_t erd__2d_coefficients_ticks[256];
uint64_t erd__2d_pq_integrals_ticks[256];
uint64_t erd__int2d_to_e0f0_ticks[256];

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD_E0F0_PCGTO_BLOCK */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD_RYS_ROOTS_WEIGHTS */
/*                ERD_2D_COEFFICIENTS */
/*                ERD_2D_PQ_INTEGRALS */
/*                ERD_INT2D_TO_E000 */
/*                ERD_INT2D_TO_E0F0 */
/*                ERD_2D_ATOM_COEFFICIENTS */
/*                ERD_2D_ATOM_PQ_INTEGRALS */
/*                ERD_ATOM_INT2D_TO_E000 */
/*                ERD_ATOM_INT2D_TO_E0F0 */
/*  DESCRIPTION : This operation calculates a batch of unnormed electron */
/*                repulsion integrals between primitive cartesian */
/*                gaussians for the shell quadruplet range: */
/*                    [E0|F0]     , E = A to P, F = C to Q */
/*                           ijkl */
/*                and the block of ij and kl exponent pairs. The total */
/*                number of eris generated here is thus given by the */
/*                total number of cartesian monomials NXYZET*NXYZFT */
/*                times the total number of exponent pairs MIJKL in the */
/*                present block. */
/*                On exit, the batch elements will be stored as: */
/*                             batch (kl,ij,nxyzt) */
/*                  Input: */
/*                    NBATCH       =  size of the primitive cartesian */
/*                                    integral batch */
/*                    NINT2D       =  space needed for each of the 2D */
/*                                    X,Y,Z integral arrays */
/*                    ATOMIC       =  indicates, if purely atomic */
/*                                    integrals will be evaluated */
/*                    ATOMAB(CD)   =  indicates, if centers A and B */
/*                                    (C and D) coincide */
/*                    MIJ(KL)      =  current # of ij (kl) primitive */
/*                                    index pairs corresponding to */
/*                                    the contracted shell pairs A,B */
/*                                    (C,D) */
/*                    MIJKL        =  current # of ijkl primitive */
/*                                    index quadruplets (= MIJ*MKL) */
/*                    NIJ          =  total # of ij primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair A,B */
/*                    NIJBEG(END)  =  first(last) ij primitive index */
/*                                    defining the ij block */
/*                    NKL          =  total # of kl primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair C,D */
/*                    NKLBEG(END)  =  first(last) kl primitive index */
/*                                    defining the kl block */
/*                    NGQP         =  # of gaussian quadrature points */
/*                                    (roots) */
/*                    NMOM         =  # of necessary moment integrals */
/*                                    to calculate the quadrature roots */
/*                    NGQSCR       =  size of gaussian quadrature */
/*                                    scratch space needed to calculate */
/*                                    all the quadrature roots */
/*                    MGQIJKL      =  # of roots times # of ijkl */
/*                                    quadruplets (= NGQP*MIJKL) */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for contraction shells x = A,B,C,D */
/*                    NXYZE(F)T    =  sum of # of cartesian monomials */
/*                                    for all shells in the range */
/*                                    E = A,...,P=A+B and in the range */
/*                                    F = C,...,Q=C+D */
/*                    NXYZP(Q)     =  # of cartesian monomials for */
/*                                    the P=A+B and Q=C+D shells */
/*                    SHELLx       =  the shell type for contraction */
/*                                    shells x = A,P=A+B,C,Q=C+D */
/*                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers */
/*                                    x = A,B,C,D */
/*                    ABm(CDm)     =  the m=x,y,z-coordinate differences */
/*                                    between centers A and B (C and D) */
/*                    ALPHAx       =  the primitive exponents for */
/*                                    contraction shells x = A,B,C,D */
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
/*                                    x = A,B,C,D */
/*                    NORMx        =  the normalization factors due to */
/*                                    the primitive exponents for the */
/*                                    contraction shells x = A,B,C,D */
/*                    RHOAB(CD)    =  the complete set of NIJ (NKL) */
/*                                    exponential prefactors between */
/*                                    contraction shells A and B */
/*                                    (C and D) */
/*                    P            =  will hold current MIJ exponent */
/*                                    sums for contraction shells A */
/*                                    and B */
/*                    Px           =  will hold current MIJ coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers P=A+B */
/*                    PAx          =  will hold current MIJ coordinate */
/*                                    x=X,Y,Z differences P-A between */
/*                                    centers P and A */
/*                    PINVHF       =  will hold current MIJ values of */
/*                                    1/(2*P), where P are the exponent */
/*                                    sums for contraction shells A */
/*                                    and B */
/*                    SCALEP       =  will hold current MIJ values of */
/*                                    scaling factors related to point P */
/*                    Q            =  will hold current MKL exponent */
/*                                    sums for contraction shells C */
/*                                    and D */
/*                    Qx           =  will hold current MKL coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers Q=C+D */
/*                    QCx          =  will hold current MKL coordinate */
/*                                    x=X,Y,Z differences Q-C between */
/*                                    centers Q and C */
/*                    QINVHF       =  will hold current MKL values of */
/*                                    1/(2*Q), where Q are the exponent */
/*                                    sums for contraction shells C */
/*                                    and D */
/*                    SCALEQ       =  will hold current MKL values of */
/*                                    scaling factors related to point Q */
/*                    RTS          =  will hold all current MGQIJKL */
/*                                    quadrature roots */
/*                    WTS          =  will hold all current MGQIJKL */
/*                                    quadrature weights */
/*                    GQSCR        =  will be used as scratch space */
/*                                    for determining the quadrature */
/*                                    roots and weights */
/*                    TVAL         =  will hold current MIJKL values */
/*                                    of T-exponents defining the Rys */
/*                                    weight functions */
/*                    PQPINV       =  will hold current MIJKL values */
/*                                    of 1/(P+Q), i.e. the inverses */
/*                                    of all total exponent sums */
/*                    SCALEPQ      =  will hold current distinct MIJKL */
/*                                    (expanded to MGQIJKL) values of */
/*                                    the overal scaling factors for */
/*                                    the integrals */
/*                    Bxx          =  will hold the current MGQIJKL */
/*                                    coordinate independent VRR */
/*                                    B-coefficients (xx=00,01,10) */
/*                    C00x         =  will hold the current MGQIJKL */
/*                                    VRR C-coefficients (individual */
/*                                    cartesian components x=X,Y,Z) for */
/*                                    shell expansion on center P */
/*                    D00x         =  will hold the current MGQIJKL */
/*                                    VRR D-coefficients (individual */
/*                                    cartesian components x=X,Y,Z) for */
/*                                    shell expansion on center Q */
/*                    INT2Dx       =  will hold all current 2D integrals */
/*                                    for each cartesian component */
/*                                    (x = X,Y,Z) */
/*                  Output: */
/*                    BATCH        =  current batch of primitive */
/*                                    cartesian [E0|F0] integrals */
/* ------------------------------------------------------------------------ */
int erd__e0f0_pcgto_block (int nij, int nkl,
                           int ngqp, int nmom,
                           int nxyzet, int nxyzft,
                           int nxyzp, int nxyzq,
                           int shella, int shellp,
                           int shellc, int shellq,
                           double xa, double ya, double za,
                           double xb, double yb, double zb,
                           double xc, double yc, double zc,
                           double xd, double yd, double zd,
                           double *alphaa, double *alphab,
                           double *alphac, double *alphad,
                           double *cca, double *ccb,
                           double *ccc, double *ccd,
                           double *ftable, int mgrid, int ngrid,
                           double tmax, double tstep, double tvstep,
                           int *prima, int *primb,
                           int *primc, int *primd,
                           double *norma, double *normb,
                           double *normc, double *normd,
                           double *rhoab, double *rhocd,
                           double *p, double *px,
                           double *py, double *pz,
                           double *pax, double *pay, double *paz,
                           double *pinvhf, double *scalep, double *q,
                           double *qx, double *qy, double *qz,
                           double *qcx, double *qcy, double *qcz,
                           double *qinvhf, double *scaleq,
                           double *rts, double *wts,
                           double *gqscr, double *tval,
                           double *pqpinv, double *scalepq,
                           double *b00, double *b01, double *b10,
                           double *c00x, double *c00y, double *c00z,
                           double *d00x, double *d00y, double *d00z,
                           double *int2dx, double *int2dy,
                           double *int2dz, double *batch)
{
    int i, j, k, l, m, n;
    int ij, kl, g000, g010, g020, g030, g040, g050, g060;
    double pqx, pqy, pqz;
    double expa, expb, expc, expd, pval, qval;
    double pinv, qinv, pxval, pyval, pzval;
    int case2d;
    double pscale, invers, pqmult, pqplus;
    int nijkl;
    double abx;
    double aby;
    double abz;
    double cdx;
    double cdy;
    double cdz;
    int mgqijkl;
    uint64_t start_clock, end_clock; 
    int tid = omp_get_thread_num();
    
/*            ...predetermine 2D integral case. This is done in */
/*               order to distinguish the P- and Q-shell combinations */
/*               for efficient evaluation of the 2D integrals and */
/*               their VRR coefficients. The cases distinguished */
/*               are summarized in the following table, indicating */
/*               the value of CASE2D: */
/*                                Q-shell */
/*                              s    p   >p */
/*                            --------------- */
/*                           | */
/*                         s |  1    4    7 */
/*                           | */
/*                P-shell  p |  2    5    8 */
/*                           | */
/*                        >p |  3    6    9 */
/*                           | */
    --batch;
    --int2dz;
    --int2dy;
    --int2dx;
    --scalep;
    --pinvhf;
    --paz;
    --pay;
    --pax;
    --pz;
    --py;
    --px;
    --p;
    --primb;
    --prima;
    --scaleq;
    --qinvhf;
    --qcz;
    --qcy;
    --qcx;
    --qz;
    --qy;
    --qx;
    --q;
    --primd;
    --primc;
    --pqpinv;
    --tval;
    --rhoab;
    --rhocd;
    --gqscr;
    --d00z;
    --d00y;
    --d00x;
    --c00z;
    --c00y;
    --c00x;
    --b10;
    --b01;
    --b00;
    --scalepq;
    --wts;
    --rts;
    
/*             ...predetermine in 'K2' loops the quantities associated */
/*                with the A,B-part and C,D-part. Set the atom equality */
/*                case CASEAT here to exploit simplifications due to */
/*                center equalities later on: */
/*                     CASEAT = 1  -->    atomic (AA|AA) integrals */
/*                            = 2  -->  2-center (AA|CC) integrals */
/*                            = 3  -->  3-center (AB|CC) integrals */
/*                            = 4  -->  3-center (AA|CD) integrals */
/*                            = 5  -->  4-center (AB|CD) integrals */
    case2d = MIN(2, shellq) * 3 + MIN(2, shellp) + 1;
    nijkl = nij * nkl;
    mgqijkl = ngqp * nijkl;
    abx = xa - xb;
    aby = ya - yb;
    abz = za - zb;
    cdx = xc - xd;
    cdy = yc - yd;
    cdz = zc - zd; 
        
    for (ij = 1; ij <= nij; ++ij)
    {
        i = prima[ij];
        j = primb[ij];
        expa = alphaa[i];
        expb = alphab[j];
        pval = expa + expb;
        pinv = 1. / pval;
        p[ij] = pval;
        pval = -expb * pinv;
        pax[ij] = pval * abx;
        pay[ij] = pval * aby;
        paz[ij] = pval * abz;
        px[ij] = pax[ij] + xa;
        py[ij] = pay[ij] + ya;
        pz[ij] = paz[ij] + za;
        pinvhf[ij] = pinv * 0.5;
        scalep[ij] = norma[i] * normb[j] * rhoab[ij];
        scalep[ij] *= cca[i] * ccb[j];
    }

    for (kl = 1; kl <= nkl; ++kl)
    {
        k = primc[kl];
        l = primd[kl];
        expc = alphac[k];
        expd = alphad[l];
        qval = expc + expd;
        qinv = 1.0 / qval;
        q[kl] = qval;
        qval = -expd * qinv;
        qcx[kl] = qval * cdx;
        qcy[kl] = qval * cdy;
        qcz[kl] = qval * cdz;
        qx[kl] = qcx[kl] + xc;
        qy[kl] = qcy[kl] + yc;
        qz[kl] = qcz[kl] + zc;
        qinvhf[kl] = qinv * 0.5;
        scaleq[kl] = normc[k] * normd[l] * rhocd[kl];
        scaleq[kl] *= ccc[k] * ccd[l];
    }

/*             ...the 'K4' loop over all ij- and kl-exponent pairs */
/*                in present ij and kl block to calculate all T's */
/*                and scaling factors for the cases: */
/*                     CASEAT = 1  -->    atomic (AA|AA) integrals */
/*                            = 2  -->  2-center (AA|CC) integrals */
/*                            = 3  -->  3-center (AB|CC) integrals */
/*                            = 4  -->  3-center (AA|CD) integrals */
/*                            = 5  -->  4-center (AB|CD) integrals */
/*                4-center (AB|CD) integrals are checked first */
/*                (most common occurence in large systems). */
    m = 1;
    for (ij = 1; ij <= nij; ++ij)
    {
        pval = p[ij];
        pxval = px[ij];
        pyval = py[ij];
        pzval = pz[ij];
        pscale = scalep[ij];
        for (kl = 1; kl <= nkl; ++kl)
        {
            qval = q[kl];
            pqmult = pval * qval;
            pqplus = pval + qval;
            invers = 1. / pqplus;
            pqx = pxval - qx[kl];
            pqy = pyval - qy[kl];
            pqz = pzval - qz[kl];
            tval[m] = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult * invers;
            pqpinv[m] = invers;
            scalepq[m] = pscale * scaleq[kl] / (pqmult * sqrt (pqplus));
            ++m;
        }
    }

/*             ...if necessary, expand the scaling array size from */
/*                MIJKL to MGQIJKL starting from the last elements. */
    if (ngqp > 1)
    {
        n = mgqijkl + 1;
        for (m = nijkl; m >= 1; --m)
        {
            for (i = 1; i <= ngqp; ++i)
            {
                scalepq[n - i] = scalepq[m];
            }
            n -= ngqp;
        }
    }

/*             ...determine memory allocation offsets for the scratch */
/*                arrays used to calculate the quadrature roots + */
/*                weights: */
/*                   G000 = offset for A coefficients (Jacobi/Laguerre) */
/*                   G010 = offset for B coefficients (Jacobi/Laguerre) */
/*                   G020 = offset for moments (Jacobi/Laguerre) */
/*                   G030 = offset for diagonals of symmetric termat */
/*                   G040 = offset for offdiagonals of symmetric termat */
/*                   G050 = offset for first row intermediates during */
/*                          evaluation of symmetric termat */
/*                   G060 = offset for second row intermediates during */
/*                          evaluation of symmetric termat */
    g000 = 1;
    g010 = g000 + nmom;
    g020 = g010 + nmom - 1;
    g030 = g020 + nmom;
    g040 = g030 + ngqp;
    g050 = g040 + ngqp;
    g060 = g050 + nmom;

/*             ...calculate all roots and weights. Array B00 is passed */
/*                as a scratch array. */
    start_clock = __rdtsc();
    erd__rys_roots_weights_ (&nijkl, &mgqijkl, &ngqp, &nmom, &tval[1], &b00[1],
                             ftable, &mgrid, &ngrid,
                             &tmax, &tstep, &tvstep, &gqscr[g000], &gqscr[g010],
                             &gqscr[g020], &gqscr[g030], &gqscr[g040],
                             &gqscr[g050], &gqscr[g060], &rts[1], &wts[1]);
    end_clock = __rdtsc();
    erd__rys_roots_weights_ticks[tid] += (end_clock - start_clock);

/*             ...perform the following steps: */
/*                1) generate all VRR coefficients. */
/*                2) construct all 2D PQ x,y,z integrals using all the */
/*                   weights and all the generated VRR coefficients for */
/*                   all exponent quadruples. */
/*                3) assemble the complete [E0|F0] batch for all ij and */
/*                   kl pairs using the 2D integrals. Arrays B00 and B01 */
/*                   are passed as scratch arrays. */
/*                The last step 3) is the most compute intensive and */
/*                separate routines are provided depending on presence */
/*                of s-shells. The gainings are in the innermost loops */
/*                of these routines, which are considerably simplified */
/*                for the special s-shell cases. Note, that the case */
/*                in which both P- and Q-shells are s-shells cannot */
/*                arise, as this case is dealt with in separate routines. */
    start_clock = __rdtsc();
    erd__2d_coefficients (nij, nkl, ngqp, &p[1], &q[1],
                          &px[1], &py[1], &pz[1], &qx[1], &qy[1], &qz[1],
                          &pax[1], &pay[1], &paz[1], &qcx[1], &qcy[1], &qcz[1],
                          &pinvhf[1], &qinvhf[1], &pqpinv[1], &rts[1],
                          case2d, &b00[1], &b01[1], &b10[1],
                          &c00x[1], &c00y[1], &c00z[1],
                          &d00x[1], &d00y[1], &d00z[1]);
    end_clock = __rdtsc();
    erd__2d_coefficients_ticks[tid] += (end_clock - start_clock);

    int mgqijkl_aligned = ((mgqijkl + SIMD_WIDTH - 1)/SIMD_WIDTH) * SIMD_WIDTH;
    int int2d_size_aligned = mgqijkl_aligned * (shellp + 1) * (shellq + 1);
    double rts_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double wts_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double b00_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double b01_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double b10_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double c00x_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double c00y_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double c00z_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double d00x_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double d00y_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double d00z_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double scalepq_aligned[mgqijkl_aligned]  __attribute__((aligned(64)));
    double int2dx_aligned[int2d_size_aligned]  __attribute__((aligned(64)));
    double int2dy_aligned[int2d_size_aligned]  __attribute__((aligned(64)));
    double int2dz_aligned[int2d_size_aligned]  __attribute__((aligned(64)));
    double batch_aligned[nxyzet * nxyzft]  __attribute__((aligned(64)));
 
    memset(rts_aligned, 0, mgqijkl_aligned*sizeof(double));
    memset(wts_aligned, 0, mgqijkl_aligned*sizeof(double));
    memset(b00_aligned, 0, mgqijkl_aligned*sizeof(double));
    memset(b01_aligned, 0, mgqijkl_aligned*sizeof(double));
    memset(b10_aligned, 0, mgqijkl_aligned*sizeof(double));
    memset(c00x_aligned, 0, mgqijkl_aligned*sizeof(double));
    memset(c00y_aligned, 0, mgqijkl_aligned*sizeof(double));
    memset(c00z_aligned, 0, mgqijkl_aligned*sizeof(double));
    memset(d00x_aligned, 0, mgqijkl_aligned*sizeof(double));
    memset(d00y_aligned, 0, mgqijkl_aligned*sizeof(double));
    memset(d00z_aligned, 0, mgqijkl_aligned*sizeof(double));
    memset(scalepq_aligned, 0, mgqijkl_aligned*sizeof(double));

/*    for(j = 0; j < nijkl; j++)
    {
        for(k = 0; k < ngqp; k++)
        {
            rts_aligned[k * nijkl_aligned + j] = rts[1 + j * ngqp + k];
            wts_aligned[k * nijkl_aligned + j] = wts[1 + j * ngqp + k];
            b00_aligned[k * nijkl_aligned + j] = b00[1 + j * ngqp + k];
            b01_aligned[k * nijkl_aligned + j] = b01[1 + j * ngqp + k];
            b10_aligned[k * nijkl_aligned + j] = b10[1 + j * ngqp + k];
            c00x_aligned[k * nijkl_aligned + j] = c00x[1 + j * ngqp + k];
            c00y_aligned[k * nijkl_aligned + j] = c00y[1 + j * ngqp + k];
            c00z_aligned[k * nijkl_aligned + j] = c00z[1 + j * ngqp + k];
            d00x_aligned[k * nijkl_aligned + j] = d00x[1 + j * ngqp + k];
            d00y_aligned[k * nijkl_aligned + j] = d00y[1 + j * ngqp + k];
            d00z_aligned[k * nijkl_aligned + j] = d00z[1 + j * ngqp + k];
            scalepq_aligned[k * nijkl_aligned + j] = scalepq[1 + j * ngqp + k];
        }
    }
*/
    memcpy(rts_aligned, rts + 1, mgqijkl*sizeof(double));
    memcpy(wts_aligned, wts + 1, mgqijkl*sizeof(double));
    memcpy(b00_aligned, b00 + 1, mgqijkl*sizeof(double));
    memcpy(b01_aligned, b01 + 1, mgqijkl*sizeof(double));
    memcpy(b10_aligned, b10 + 1, mgqijkl*sizeof(double));
    memcpy(c00x_aligned, c00x + 1, mgqijkl*sizeof(double));
    memcpy(c00y_aligned, c00y + 1, mgqijkl*sizeof(double));
    memcpy(c00z_aligned, c00z + 1, mgqijkl*sizeof(double));
    memcpy(d00x_aligned, d00x + 1, mgqijkl*sizeof(double));
    memcpy(d00y_aligned, d00y + 1, mgqijkl*sizeof(double));
    memcpy(d00z_aligned, d00z + 1, mgqijkl*sizeof(double));
    memcpy(scalepq_aligned, scalepq + 1, mgqijkl*sizeof(double));


    start_clock = __rdtsc();
    /*erd__2d_pq_integrals (shellp, shellq, mgqijkl, &wts[1],
                          &b00[1], &b01[1], &b10[1],
                          &c00x[1], &c00y[1], &c00z[1], &d00x[1],
                          &d00y[1], &d00z[1], case2d,
                          &int2dx[1], &int2dy[1], &int2dz[1]);
    */
    erd__2d_pq_integrals (shellp, shellq, mgqijkl_aligned, wts_aligned,
                          b00_aligned, b01_aligned, b10_aligned,
                          c00x_aligned, c00y_aligned, c00z_aligned, d00x_aligned,
                          d00y_aligned, d00z_aligned, case2d,
                          int2dx_aligned, int2dy_aligned, int2dz_aligned);

    end_clock = __rdtsc();
    erd__2d_pq_integrals_ticks[tid] += (end_clock - start_clock);

    start_clock = __rdtsc();
    /*if (shellq == 0)
    {
        erd__int2d_to_e000 (shella, shellp, ngqp, nijkl, mgqijkl,
                            nxyzet, nxyzp,
                            &int2dx[1], &int2dy[1], &int2dz[1],
                            &b00[1], &b01[1], &scalepq[1], &batch[1]);
    }
    else if (shellp == 0)
    {
        erd__int2d_to_e000 (shellc, shellq, ngqp, nijkl, mgqijkl,
                            nxyzft, nxyzq,
                            &int2dx[1], &int2dy[1], &int2dz[1],
                            &b00[1], &b01[1], &scalepq[1], &batch[1]);
    }
    else
    {
        erd__int2d_to_e0f0 (shella, shellp, shellc, shellq, ngqp,
                            nijkl, mgqijkl, nxyzet, nxyzft, nxyzp, nxyzq,
                            &int2dx[1], &int2dy[1], &int2dz[1],
                            &b00[1], &b01[1], &scalepq[1], &batch[1]);
    }
    */
    erd__int2d_to_e0f0 (shella, shellp, shellc, shellq,
            mgqijkl_aligned, nxyzet, nxyzft, nxyzp, nxyzq,
            int2dx_aligned, int2dy_aligned, int2dz_aligned,
            b00_aligned, b01_aligned, scalepq_aligned, batch_aligned);
    end_clock = __rdtsc();
    erd__int2d_to_e0f0_ticks[tid] += (end_clock - start_clock);
    memcpy(&batch[1], batch_aligned, nxyzet * nxyzft * sizeof(double)); 

    return 0;
}
