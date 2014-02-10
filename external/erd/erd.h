#ifndef __ERD_H__
#define __ERD_H__

#include <stdint.h>
#include <math.h>
#include <omp.h>
#include <yepPredefines.h>

#ifdef __ERD_PROFILE__
#include "erd_profile.h"
#endif

#define MAX(a,b)    ((a) < (b) ? (b) : (a))
#define MIN(a,b)    ((a) > (b) ? (b) : (a))
#define PREFACT     9.027033336764101
#if defined (__MIC__)
#define SIMDW      8
#elif defined (__AVX__)
#define SIMDW      4
#elif defined (__SSE__)
#define SIMDW      2
#else
#define SIMDW      8
#endif

#define PAD_LEN(N)  ((N+SIMDW-1)/SIMDW * SIMDW )
#define PAD_LEN2(N) ((N+SIMDW*2-1)/(SIMDW*2) * SIMDW*2 )

/*******************************************************************/
// C functions

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

int erd__move_ry (
		  int nindex, int notmove, int move, int nry,
                  int index, double *x,
                  int *ixoff, double *y);

int erd__set_ij_kl_pairs (int npgtoa, int npgtob, int npgtoc, int npgtod,
                          double xa, double ya, double za,
                          double xb, double yb, double zb,
                          double xc, double yc, double zc,
                          double xd, double yd, double zd,
                          double rnabsq, double rncdsq, double prefact,
                          double *YEP_RESTRICT alphaa,
                          double *YEP_RESTRICT alphab,
                          double *YEP_RESTRICT alphac,
                          double *YEP_RESTRICT alphad,
                          int screen, int *YEP_RESTRICT empty,
                          int *YEP_RESTRICT nij_ptr,
                          int *YEP_RESTRICT nkl_ptr, int *YEP_RESTRICT prima,
                          int *YEP_RESTRICT primb, int *YEP_RESTRICT primc,
                          int *YEP_RESTRICT primd, double *YEP_RESTRICT rho);
	
int erd__map_ijkl_to_ikjl (int ni, int nj, int nk, int nl,
                           double *x, double *y);

int erd__prepare_ctr (int npgtoa, int npgtob,
                      int npgtoc, int npgtod,
                      int shella, int shellb,
                      int shellc, int shelld,
                      double *alphaa, double *alphab,
                      double *alphac, double *alphad,
                      double spnorm,
                      double *norma, double *normb,
                      double *normc, double *normd);

int erd__ctr_1st_half (int n, int mij, double *ccr, double *ccs,
                       int *primr, int *prims,
                       int equalrs, double *x, double *y);

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
                          int nij, int nkl, int nxyzt,
                          int memory, int *zcbatch, 
                          int *znorm, int *zrho12, int *zrho34,
                          int *zp, int *zpx, int *zpy, int *zpz,
                          int *zscpk2, int *zq, int *zqx,
                          int *zqy, int *zqz, int *zscqk2);

void erd__pppp_pcgto_block (int nij, int nkl,
                           double x1, double y1, double z1,
                           double x2, double y2, double z2,
                           double x3, double y3, double z3,
                           double x4, double y4, double z4,
                           double *alpha1, double *alpha2,
                           double *alpha3, double *alpha4,
                           double *cc1, double *cc2,
                           double *cc3, double *cc4,
                           int *prim1, int *prim2,
                           int *prim3, int *prim4,
                           double *norm1, double *norm2,
                           double *norm3, double *norm4,
                           double *rho12, double *rho34,
                           double *p, double *px,
                           double *py, double *pz, double *scalep,
                           double *q, double *qx,
                           double *qy, double *qz,
                           double *scaleq, double *cbatch);

void erd__sppp_pcgto_block (int nij, int nkl,
                           int shell1, int shell3, int shellp,
                           double x1, double y1, double z1,
                           double x2, double y2, double z2,
                           double x3, double y3, double z3,
                           double x4, double y4, double z4,
                           double *alpha1, double *alpha2,
                           double *alpha3, double *alpha4,
                           double *cc1, double *cc2,
                           double *cc3, double *cc4,
                           int *prim1, int *prim2,
                           int *prim3, int *prim4,
                           double *norm1, double *norm2,
                           double *norm3, double *norm4,
                           double *rho12, double *rho34,
                           double *p, double *px,
                           double *py, double *pz, double *scalep,
                           double *q, double *qx,
                           double *qy, double *qz,
                           double *scaleq, double *cbatch);

void erd__sspp_pcgto_block (int nij, int nkl,
                           int shell1, int shell3, int shellp,
                           double x1, double y1, double z1,
                           double x2, double y2, double z2,
                           double x3, double y3, double z3,
                           double x4, double y4, double z4,
                           double *alpha1, double *alpha2,
                           double *alpha3, double *alpha4,
                           double *cc1, double *cc2,
                           double *cc3, double *cc4,
                           int *prim1, int *prim2,
                           int *prim3, int *prim4,
                           double *norm1, double *norm2,
                           double *norm3, double *norm4,
                           double *rho12, double *rho34,
                           double *p, double *px,
                           double *py, double *pz, double *scalep,
                           double *q, double *qx,
                           double *qy, double *qz,
                           double *scaleq, double *cbatch);

void erd__sssp_pcgto_block(int nij, int nkl,
                           int shell1, int shell3, int shellp,
                           double x1, double y1, double z1,
                           double x2, double y2, double z2,
                           double x3, double y3, double z3,
                           double x4, double y4, double z4,
                           double *alpha1, double *alpha2,
                           double *alpha3, double *alpha4,
                           double *cc1, double *cc2,
                           double *cc3, double *cc4,
                           int *prim1, int *prim2,
                           int *prim3, int *prim4,
                           double *norm1, double *norm2,
                           double *norm3, double *norm4,
                           double *rho12, double *rho34,
                           double *p, double *px,
                           double *py, double *pz, double *scalep,
                           double *q, double *qx,
                           double *qy, double *qz,
                           double *scaleq, double *cbatch);

void erd__ssss_pcgto_block (int nij, int nkl,
                           double x1, double y1, double z1,
                           double x2, double y2, double z2,
                           double x3, double y3, double z3,
                           double x4, double y4, double z4,
                           double *alpha1, double *alpha2,
                           double *alpha3, double *alpha4,
                           double *cc1, double *cc2,
                           double *cc3, double *cc4,
                           int *prim1, int *prim2,
                           int *prim3, int *prim4,
                           double *norm1, double *norm2,
                           double *norm3, double *norm4,
                           double *rho12, double *rho34,
                           double *p, double *px, double *py, double *pz,
                           double *scalep, double *q, double *qx,
                           double *qy, double *qz,
                           double *scaleq, double *cbatch);

int erd__hrr_transform (int m, int nrow,
                        int nxyza, int nxyzb,
                        int *lrow, int *row,
                        double *rot, double *x, double *y);

int erd__xyz_to_ry_abcd (int nxyza, int nxyzb, int nxyzc, int nxyzd,
                         int nrya, int nryb, int nryc, int nryd,
                         int shella, int shellb,
                         int shellc, int shelld,
                         int istart, int zstart,
                         int *nrowa, int *nrowb,
                         int *nrowc, int *nrowd,
                         int *nrota, int *nrotb,
                         int *nrotc, int *nrotd,
                         int *z00a, int *z00b, int *z00c, int *z00d,
                         int *i0a1, int *i0b1, int *i0c1, int *i0d1,
                         int *i0a2, int *i0b2, int *i0c2, int *i0d2,
                         int *iused, int *zused,
                         int *icore, double *zcore);

int erd__xyz_to_ry_matrix (int nxyz, int nrowmx, int l,
                           double *temp, int *nrow,
                           int *row, double *tmat);

int erd__spherical_transform (int m, int nrow, int nry,
                              int *lrow, int *row, double *rot,
                              double *x, double *y);

int erd__hrr_step (int nabo, int mrowin,
                   int mrowout, int nxyzx,
                   int nxyza, int nxyzb, int nxyzao,
                   int shellx, int shellp, int shellb,
                   double abx, double aby, double abz,
                   int *cpair, int *nrowin, int *rowin, double *win,
                   int *nrowout, int *rowout, double *wout);

int erd__hrr_matrix (int nrothrr, int ncolhrr,
                     int nxyzet, int nxyza, int nxyzp,
                     int shella, int shellb, int shellp,
                     int nabcoor, double abx, double aby, double abz,
                     int *work, int *in1, int *in2,
                     int *nrowout, int *nrow,
                     int *row, double *t);

double erd__dsqmin_line_segments (double xp0, double yp0,
                                  double zp0, double xp1,
                                  double yp1, double zp1,
                                  double xq0, double yq0,
                                  double zq0, double xq1,
                                  double yq1, double zq1);

int erd__set_abcd (int npgto1, int npgto2, int npgto3, int npgto4,
                   int shell1, int shell2, int shell3, int shell4,
                   double x1, double y1, double z1,
                   double x2, double y2, double z2, 
                   double x3, double y3, double z3,
                   double x4, double y4, double z4, int spheric,
                   int *npgtoa, int *npgtob, int *npgtoc, int *npgtod,
                   int *shella, int *shellb, int *shellc, int *shelld,
                   double *xa, double *ya, double *za,
                   double *xb, double *yb, double *zb,
                   double *xc, double *yc, double *zc,
                   double *xd, double *yd, double *zd,
                   int *nxyza, int *nxyzb, int *nxyzc, int *nxyzd,
                   int *nxyzet, int *nxyzft,
                   int *nrya, int *nryb, int *nryc, int *nryd,
                   int *nabcoor, int *ncdcoor,
                   int *ncolhrr, int *nrothrr,
                   int *nxyzhrr, int *empty, int *tr1234);                  

int erd__normalize_cartesian (int m, int l, double *norm, double *batch);

int erd__cartesian_norms (int l, double *norm);

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
                           double *int2dz, double *batch);

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
                          double *d00x, double *d00y, double *d00z);

int erd__2d_pq_integrals (int shellp, int shellq, int ngqexq,
                          double *wts, double *b00, double *b01, double *b10,
                          double *c00x, double *c00y, double *c00z,
                          double *d00x, double *d00y, double *d00z,
                          int case2d, double *int2dx,
                          double *int2dy, double *int2dz);

int erd__int2d_to_e0f0 (int shella, int shellp, int shellc, int shellq,
                        int ngqexq,
                        int nxyzet, int nxyzft, int nxyzp, int nxyzq,
                        double *int2dx, double *int2dy, double *int2dz,
                        double *scale, double *batch);

int erd__int2d_to_e000 (int shella, int shellp, int ngqp, int nexq, int ngqexq,
                        int nxyzet, int nxyzp,
                        double *int2dx, double *int2dy, double *int2dz,
                        double *temp1, double *temp2,
                        double *scale, double *batch);

int erd__e0f0_def_blocks (int zmax, int npgtoa, int npgtob,
                          int npgtoc, int npgtod,
                          int shellp, int shellq,
                          int nij, int nkl, int ngqp, int ngqscr,
                          int nxyzt, int memory, int *nint2d, int *zcbatch,
                          int *znorm, int *zrhoab, int *zrhocd,
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
                          
int erd__rys_roots_weights_ (int * nt, int * ngqp, int * nmom, double * tval,
                             double * ryszero,
                             double * a, double * b, double * mom,
                             double * dia, double * off,
                             double * row1, double * row2,
                             double * rts, double * wts);


int erd__rys_1_roots_weights_ (int * nt, double * tval,
                               double * rts, double * wts);

int erd__rys_2_roots_weights_ (int * nt, double * tval, double * rts,
                               double * wts);

int erd__rys_3_roots_weights_ (int * nt, double * tval, double * rts,
                               double * wts);

int erd__rys_4_roots_weights_ (int * nt, double * tval, double * rts,
                               double * wts);

int erd__rys_5_roots_weights_ (int * nt, double * tval, double * rts,
                               double * wts);

int erd__rys_x_roots_weights_ (int * nt, int * ngqp, int * nmom, double * tval,
                               double * ryszero, double * a,
                               double * b, double * mom,
                               double * dia, double * off,
                               double * row1, double * row2,
                               double * rts, double * wts);

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif


#endif /* __ERD_H__ */
