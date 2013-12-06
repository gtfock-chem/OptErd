#ifndef __ERD_H__
#define __ERD_H__


#include <yepPredefines.h>


#define MAX(a,b)    ((a) < (b) ? (b) : (a))
#define MIN(a,b)    ((a) > (b) ? (b) : (a))
#define PREFACT     9.027033336764101


/*******************************************************************/
// C functions

int erd__move_ry (int nindex, int notmove, int move, int nry,
                  int index, double *x,
                  int *ixoff, double *y);

int erd__set_ij_kl_pairs (int npgtoa, int npgtob, int npgtoc, int npgtod,
                          int atomab, int atomcd, int equalab, int equalcd,
                          double xa, double ya, double za,
                          double xb, double yb, double zb,
                          double xc, double yc, double zc,
                          double xd, double yd, double zd,
                          double rnabsq, double rncdsq, double prefact,
                          double *YEP_RESTRICT alphaa,
                          double *YEP_RESTRICT alphab,
                          double *YEP_RESTRICT alphac,
                          double *YEP_RESTRICT alphad,
                          double *YEP_RESTRICT ftable, int mgrid, int ngrid,
                          double tmax, double tstep, double tvstep, int screen,
                          int *YEP_RESTRICT empty, int *YEP_RESTRICT nij_ptr,
                          int *YEP_RESTRICT nkl_ptr, int *YEP_RESTRICT prima,
                          int *YEP_RESTRICT primb, int *YEP_RESTRICT primc,
                          int *YEP_RESTRICT primd, double *YEP_RESTRICT rho);
	
int erd__map_ijkl_to_ikjl (int ni, int nj, int nk, int nl,
                           double *x, double *y);

int erd__prepare_ctr (int ncsize, int nij, int nkl,
                      int npgtoa, int npgtob,
                      int npgtoc, int npgtod,
                      int shella, int shellb,
                      int shellc, int shelld,
                      double *alphaa, double *alphab,
                      double *alphac, double *alphad, double spnorm,
                      int equalab, int equalcd, double *rho,
                      double *norma, double *normb,
                      double *normc, double *normd,
                      double *rhoab, double *rhocd, double *cbatch);

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
                          int memory, int *ncsize,
                          int *zcbatch, int *zpbatch, int *zwork, 
                          int *znorm1, int *znorm2,
                          int *znorm3, int *znorm4,
                          int *zrho12, int *zrho34,
                          int *zp, int *zpx, int *zpy, int *zpz,
                          int *zscpk2, int *zq, int *zqx,
                          int *zqy, int *zqz, int *zscqk2);

int erd__pppp_pcgto_block (int nij, int nkl,
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
                           double *scaleq, double *batch);

int erd__sppp_pcgto_block (int atomic, int nij, int nkl,
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
                           double *scaleq, double *batch);

int erd__sspp_pcgto_block (int nij, int nkl,
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
                           double *scaleq, double *batch);

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
                           double *scaleq, double *batch);

int erd__ssss_pcgto_block (int nij, int nkl,
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
                           double *p, double *px, double *py, double *pz,
                           double *scalep, double *q, double *qx,
                           double *qy, double *qz,
                           double *scaleq, double *batch);

int erd__ctr_4index_block (int nxyzt, int mij, int mkl,
                           double *ccr, double *ccs,
                           double *cct, double *ccu,
                           int *primr, int *prims,
                           int *primt, int *primu,
                           int equalrs, int equaltu, int ptrans,
                           double *pbatch, double *work, double *cbatch);

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

int erd__xyz_to_ry_matrix (int nxyz, int nry,
                           int nrowmx, int l,
                           double *temp, int *nrow,
                           int *row, double *tmat);

int erd__spherical_transform (int m, int nrow, int nry,
                              int *lrow, int *row, double *rot,
                              double *x, double *y);

int erd__hrr_step (int nab, int nabo, int mrowin,
                   int mrowout, int nxyzx, int nxyzp,
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

int erd__memory_1111_csgto (int npgto1, int npgto2,
                            int npgto3, int npgto4,
                            int shell1, int shell2,
                            int shell3, int shell4,
                            double x1, double y1, double z1,
                            double x2, double y2, double z2,
                            double x3, double y3, double z3,
                            double x4, double y4, double z4,
                            double *alpha, double *cc,
                            int *imin, int *iopt, int *zmin, int *zopt);

int erd__set_abcd (int npgto1, int npgto2, int npgto3, int npgto4,
                   int shell1, int shell2, int shell3, int shell4,
                   double x1, double y1, double z1,
                   double x2, double y2, double z2, 
                   double x3, double y3, double z3,
                   double x4, double y4, double z4,
                   double *exp1, double *exp2,
                   double *exp3, double *exp4,
                   double *cc1, double *cc2,
                   double *cc3, double *cc4, int spheric,
                   int *npgtoa, int *npgtob,
                   int *npgtoc, int *npgtod,
                   int *shella, int *shellb, int *shellc, int *shelld,
                   int *shellp, int *shellq, int *shellt, int *mxshell,
                   double *xa, double *ya, double *za,
                   double *xb, double *yb, double *zb,
                   double *xc, double *yc, double *zc,
                   double *xd, double *yd, double *zd,
                   int *atomic, int *atomab, int *atomcd,
                   int *equalab, int *equalcd,
                   double *abx, double *aby,
                   double *abz, double *cdx, double *cdy,
                   double *cdz, int *nabcoor, int *ncdcoor,
                   double *rnabsq, double *rncdsq,
                   double *spnorm, int *nxyza, int *nxyzb,
                   int *nxyzc, int *nxyzd, int *nxyzet,
                   int *nxyzft, int *nxyzp, int *nxyzq,
                   int *nrya, int *nryb, int *nryc,
                   int *nryd, int *indexa, int *indexb,
                   int *indexc, int *indexd, int *swap12,
                   int *swap34, int *swaprs, int *swaptu,
                   int *tr1234, int *lexpa, int *lexpb,
                   int *lexpc, int *lexpd, int *lcca,
                   int *lccb, int *lccc, int *lccd,
                   int *lccsega, int *lccsegb, int *lccsegc,
                   int *lccsegd, int *nxyzhrr, int *ncolhrr,
                   int *nrothrr, int *empty);

int erd__e0f0_def_blocks (int zmax, int npgtoa, int npgtob,
                          int npgtoc, int npgtod,
                          int shellp, int shellq,
                          int nij, int nkl, int ngqp, int ngqscr,
                          int nxyzt, int memory,
                          int *npsize, int *ncsize,
                          int *nint2d, int *zcbatch,
                          int *zpbatch, int *zwork,
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
                          
int erd__memory_csgto (int npgto1, int npgto2,
                       int npgto3, int npgto4,
                       int shell1, int shell2, int shell3, int shell4,
                       double x1, double y1, double z1,
                       double x2, double y2, double z2,
                       double x3, double y3, double z3,
                       double x4, double y4, double z4,
                       double *alpha, double *cc, int spheric,
                       int *imin, int *iopt,
                       int *zmin, int *zopt);

int erd__normalize_cartesian (int m, int l, double *norm, double *batch);

int erd__cartesian_norms (int l, double *norm);


/*******************************************************************/
// Fortran functions

int erd__e0f0_pcgto_block_ (int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, int *, int *, double *, 
	    double *, double *, int *, int *, int *, 
	    int *, double *, double *, double *, double *,
	     double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *);


#endif /* __ERD_H__ */
