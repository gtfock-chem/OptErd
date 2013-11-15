#ifndef __ERD_H__
#define __ERD_H__


#include <yepPredefines.h>


#define MAX(a,b)    ((a) < (b) ? (b) : (a))
#define MIN(a,b)    ((a) > (b) ? (b) : (a))
#define PREFACT     9.027033336764101

	
int erd__map_ijkl_to_ikjl (int ni, int nj, int nk, int nl,
                           double *x, double *y);

int erd__prepare_ctr (int ncsize, int nij, int nkl,
                      int npgtoa, int npgtob,
                      int npgtoc, int npgtod,
                      int shella, int shellb,
                      int shellc, int shelld,
                      double *alphaa, double *alphab,
                      double *alphac, double *alphad,
                      double prefact, double spnorm,
                      int equalab, int equalcd,
                      int blocked, double *rho,
                      double *norma, double *normb,
                      double *normc, double *normd,
                      double *rhoab, double *rhocd, double *cbatch);

int erd__ctr_1st_half (int n, int npmax, int npmin,
                       int mij, int nrs, int nblock,
                       int ncr, int ncs, int npr, int nps,
                       double *ccr, double *ccs,
                       int *ccbegr, int *ccbegs, int *ccendr,
                       int *ccends, int *primr, int *prims,
                       int equalrs, int swaprs, int *pused,
                       int *psave, int *ppair,
                       double *x, double *w, double *y);

int erd__ctr_2nd_half_new (int n, int npmax, int npmin,
                           int mkl, int ntu, int nblock,
                           int nct, int ncu, int npt, int npu,
                           double *cct, double *ccu,
                           int *ccbegt, int *ccbegu,
                           int *ccendt, int *ccendu,
                           int *primt, int *primu,
                           int equaltu, int swaptu,
                           int *pused, int *psave, int *ppair,
                           double *x, double *w, double *y);

int erd__transpose_batch (int nrow, int ncol, double *batch, double *obatch);

int erd__1111_def_blocks (int zmax, int npgto1, int npgto2,
                          int npgto3, int npgto4,
                          int nij, int nkl, int nrs, int ntu,
                          int nrstu, int nxyzt,
                          int l1cache, int nctrow, int memory,
                          int *nijblk, int *nklblk,
                          int *npsize, int *ncsize, int *nwsize,
                          int *mxprim, int *mnprim,
                          int *zcbatch, int *zpbatch, int *zwork,
                          int *znorm1, int *znorm2,
                          int *znorm3, int *znorm4,
                          int *zrho12, int *zrho34,
                          int *zp, int *zpx, int *zpy, int *zpz,
                          int *zscpk2, int *zq, int *zqx,
                          int *zqy, int *zqz, int *zscqk2);

int erd__pppp_pcgto_block (int nbatch, int atomic, int atom12, int atom34,
                           int mij, int mkl,
                           int nij, int nijbeg, int nijend,
                           int nkl, int nklbeg, int nklend,
                           int npgto1, int npgto2,
                           int npgto3, int npgto4,
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
                           double *scaleq, double *batch);

int erd__sppp_pcgto_block (int nbatch, int atomic, int atom12, int atom34,
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
                           double *scaleq, double *batch);

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
                           double *scaleq, double *batch);

int erd__sssp_pcgto_block (int nbatch, int atomic, int atom12, int atom34,
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
                           double *scaleq, double *batch);

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
                           double *scaleq, double *batch);

int erd__ctr_4index_block (int npsize, int ncsize, int nwsize,
                           int nxyzt, int mijkl,
                           int mij, int mkl, int nrs, int ntu,
                           int npr, int nps, int npt, int npu,
                           int ncr, int ncs, int nct, int ncu,
                           int mxprim, int mnprim,
                           double *ccr, double *ccs,
                           double *cct, double *ccu,
                           int *ccbegr, int *ccbegs,
                           int *ccbegt, int *ccbegu,
                           int *ccendr, int *ccends,
                           int *ccendt, int *ccendu,
                           int *primr, int *prims,
                           int *primt, int *primu,
                           int l1cache, int tile, int nctrow,
                           int equalrs, int equaltu,
                           int swaprs, int swaptu,
                           int ptrans, int blocked,
                           int *pused, int *psave,
                           int *ppair, double *pbatch,
                           double *work, double *cbatch);

int erd__normalize_cartesian_ (int *, int *, int *, double *, double *);

int erd__spherical_transform_ (int *, int *, int *,
                               int *, int *, int *,
                               double *, double *, double *);

int erd__e0f0_def_blocks (int zmax, int npgtoa, int npgtob,
                          int npgtoc, int npgtod,
                          int shellp, int shellq,
                          int nij, int nkl, int nrs, int ntu,
                          int nrstu, int ngqp, int ngqscr,
                          int nxyzt, int l1cache, int nctrow,
                          int memory,
                          int *nijblk, int *nklblk,
                          int *npsize, int *ncsize, int *nwsize,
                          int *nint2d, int *mxprim, int *mnprim,
                          int *zcbatch, int *zpbatch, int *zwork,
                          int *znorma, int *znormb, int *znormc,
                          int *znormd, int *zrhoab, int *zrhocd,
                          int *zp, int *zpx, int *zpy,
                          int *zpz, int *zpax, int *zpay,
                          int *zpaz, int *zpinvhf, int *zscpk2,
                          int *zq, int *zqx, int *zqy,
                          int *zqz, int *zqcx, int *zqcy,
                          int *zqcz, int *zqinvhf, int *zscqk2,
                          int *zrts, int *zwts, int *zgqscr,
                          int *ztval, int *zpqpinv, int *zscpqk4,
                          int *zb00, int *zb01, int *zb10,
                          int *zc00x, int *zc00y, int *zc00z,
                          int *zd00x, int *zd00y, int *zd00z,
                          int *zint2dx, int *zint2dy, int *zint2dz);

int erd__set_abcd (int ncgto1, int ncgto2,
                   int ncgto3, int ncgto4,
                   int npgto1, int npgto2,
                   int npgto3, int npgto4,
                   int shell1, int shell2,
                   int shell3, int shell4,
                   double x1, double y1, double z1,
                   double x2, double y2, double z2,
                   double x3, double y3, double z3,
                   double x4, double y4, double z4,
                   double *exp1, double *exp2,
                   double *exp3, double *exp4,
                   double *cc1, double *cc2,
                   double *cc3, double *cc4,
                   int spheric, int *ncgtoa, int *ncgtob,
                   int *ncgtoc, int *ncgtod,
                   int *npgtoa, int *npgtob,
                   int *npgtoc, int *npgtod,
                   int *shella, int *shellb,
                   int *shellc, int *shelld,
                   int *shellp, int *shellq, int *shellt,
                   int *mxshell,
                   double *xa, double *ya, double *za,
                   double *xb, double *yb, double *zb,
                   double *xc, double *yc, double *zc,
                   double *xd, double *yd, double *zd,
                   int *atomic, int *atomab, int *atomcd,
                   int *equalab, int *equalcd,
                   double *abx, double *aby, double *abz,
                   double *cdx, double *cdy, double *cdz,
                   int *nabcoor, int *ncdcoor,
                   double *rnabsq, double *rncdsq,
                   double *spnorm, int *nxyza, int *nxyzb,
                   int *nxyzc, int *nxyzd,
                   int *nxyzet, int *nxyzft,
                   int *nxyzp, int *nxyzq,
                   int *nrya, int *nryb,
                   int *nryc, int *nryd,
                   int *indexa, int *indexb,
                   int *indexc, int *indexd,
                   int *swap12, int *swap34,
                   int *swaprs, int *swaptu,
                   int *tr1234, int *lexpa, int *lexpb,
                   int *lexpc, int *lexpd,
                   int *lcca, int *lccb, int *lccc, int *lccd,
                   int *lccsega, int *lccsegb,
                   int *lccsegc, int *lccsegd,
                   int *nxyzhrr, int *ncolhrr, int *nrothrr, int *empty);

/*******************************************************************/


#endif /* __ERD_H__ */