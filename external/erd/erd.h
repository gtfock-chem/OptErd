#ifndef __ERD_H__
#define __ERD_H__

#include <stdint.h>
#include <stdbool.h>
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

void erd__move_ry(uint32_t nindex, uint32_t notmove, uint32_t move, uint32_t nry,
    uint32_t index, const double x[restrict static notmove * move * nry],
    uint32_t ixoff[restrict static nindex], double y[restrict static notmove * move * nry]);

void erd__set_ij_kl_pairs(
    uint32_t npgtoa, uint32_t npgtob, uint32_t npgtoc, uint32_t npgtod,
    double xa, double ya, double za,
    double xb, double yb, double zb,
    double xc, double yc, double zc,
    double xd, double yd, double zd,
    double rnabsq, double rncdsq, double prefact,
    const double alphaa[restrict static npgtoa],
    const double alphab[restrict static npgtob],
    const double alphac[restrict static npgtoc],
    const double alphad[restrict static npgtod],
    uint32_t nij_ptr[restrict static 1], uint32_t nkl_ptr[restrict static 1],
    uint32_t prima[restrict static npgtoa*npgtob], uint32_t primb[restrict static npgtoa*npgtob], uint32_t primc[restrict static npgtoc*npgtod], uint32_t primd[restrict static npgtoc*npgtod],
    double rhoab[restrict static npgtoa*npgtob],
    double rhocd[restrict static npgtoc*npgtod]);

int erd__transpose_batch(int nrow, int ncol, double *batch, double *obatch);

void erd__1111_def_blocks(uint32_t zmax, uint32_t npgto1, uint32_t npgto2,
                          uint32_t npgto3, uint32_t npgto4,
                          uint32_t nij, uint32_t nkl, uint32_t nxyzt,
                          bool memory, uint32_t *zcbatch, 
                          uint32_t *znorm, uint32_t *zrho12, uint32_t *zrho34,
                          uint32_t *zp, uint32_t *zpx, uint32_t *zpy, uint32_t *zpz,
                          uint32_t *zscpk2, uint32_t *zq, uint32_t *zqx,
                          uint32_t *zqy, uint32_t *zqz, uint32_t *zscqk2);

void erd__e0f0_def_blocks(int zmax, int npgtoa, int npgtob,
                          int npgtoc, int npgtod,
                          int shellp, int shellq,
                          int nij, int nkl, int ngqp, int ngqscr,
                          int nxyzt, int memory, int *nint2d, int *zcbatch);

void erd__pppp_pcgto_block(uint32_t nij, uint32_t nkl,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double alpha1[restrict static nij], const double alpha2[restrict static nij], const double alpha3[restrict static nkl], const double alpha4[restrict static nkl],
    const double cc1[restrict static nij], const double cc2[restrict static nij], const double cc3[restrict static nkl], const double cc4[restrict static nkl],
    const uint32_t prim1[restrict static nij], const uint32_t prim2[restrict static nij], const uint32_t prim3[restrict static nkl], const uint32_t prim4[restrict static nkl],
    const double norm1[restrict static nij], const double norm2[restrict static nij], const double norm3[restrict static nkl], const double norm4[restrict static nkl],
    const double rho12[restrict static nij], const double rho34[restrict static nkl],
    double *restrict cbatch);

void erd__pppp_pcgto_block(uint32_t nij, uint32_t nkl,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double alpha1[restrict static nij], const double alpha2[restrict static nij], const double alpha3[restrict static nkl], const double alpha4[restrict static nkl],
    const double cc1[restrict static nij], const double cc2[restrict static nij], const double cc3[restrict static nkl], const double cc4[restrict static nkl],
    const uint32_t prim1[restrict static nij], const uint32_t prim2[restrict static nij], const uint32_t prim3[restrict static nkl], const uint32_t prim4[restrict static nkl],
    const double norm1[restrict static nij], const double norm2[restrict static nij], const double norm3[restrict static nkl], const double norm4[restrict static nkl],
    const double rho12[restrict static nij], const double rho34[restrict static nkl],
    double cbatch[restrict static 81]);

void erd__sppp_pcgto_block(uint32_t nij, uint32_t nkl,
    uint32_t shell1, uint32_t shell3, uint32_t shellp,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double alpha1[restrict static nij], const double alpha2[restrict static nij], const double alpha3[restrict static nkl], const double alpha4[restrict static nkl],
    const double cc1[restrict static nij], const double cc2[restrict static nij], const double cc3[restrict static nkl], const double cc4[restrict static nkl],
    const uint32_t prim1[restrict static nij], const uint32_t prim2[restrict static nij], const uint32_t prim3[restrict static nkl], const uint32_t prim4[restrict static nkl],
    const double norm1[restrict static nij], const double norm2[restrict static nij], const double norm3[restrict static nkl], const double norm4[restrict static nkl],
    const double rho12[restrict static nij], const double rho34[restrict static nkl],
    double cbatch[restrict static 27]);

void erd__sspp_pcgto_block(uint32_t nij, uint32_t nkl,
    uint32_t shell1, uint32_t shell3, uint32_t shellp,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double alpha1[restrict static nij], const double alpha2[restrict static nij], const double alpha3[restrict static nkl], const double alpha4[restrict static nkl],
    const double cc1[restrict static nij], const double cc2[restrict static nij], const double cc3[restrict static nkl], const double cc4[restrict static nkl],
    const uint32_t prim1[restrict static nij], const uint32_t prim2[restrict static nij], const uint32_t prim3[restrict static nkl], const uint32_t prim4[restrict static nkl],
    const double norm1[restrict static nij], const double norm2[restrict static nij], const double norm3[restrict static nkl], const double norm4[restrict static nkl],
    const double rho12[restrict static nij], const double rho34[restrict static nkl],
    double cbatch[restrict static 9]);

void erd__sssp_pcgto_block(uint32_t nij, uint32_t nkl,
    uint32_t shell1, uint32_t shell3, uint32_t shellp,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double alpha1[restrict static nij], const double alpha2[restrict static nij], const double alpha3[restrict static nkl], const double alpha4[restrict static nkl],
    const double cc1[restrict static nij], const double cc2[restrict static nij], const double cc3[restrict static nkl], const double cc4[restrict static nkl],
    const uint32_t prim1[restrict static nij], const uint32_t prim2[restrict static nij], const uint32_t prim3[restrict static nkl], const uint32_t prim4[restrict static nkl],
    const double norm1[restrict static nij], const double norm2[restrict static nij], const double norm3[restrict static nkl], const double norm4[restrict static nkl],
    const double rho12[restrict static nij], const double rho34[restrict static nkl],
    double cbatch[restrict static 3]);

void erd__ssss_pcgto_block(uint32_t nij, uint32_t nkl,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double alpha1[restrict static nij], const double alpha2[restrict static nij], const double alpha3[restrict static nkl], const double alpha4[restrict static nkl],
    const double cc1[restrict static nij], const double cc2[restrict static nij], const double cc3[restrict static nkl], const double cc4[restrict static nkl],
    const uint32_t prim1[restrict static nij], const uint32_t prim2[restrict static nij], const uint32_t prim3[restrict static nkl], const uint32_t prim4[restrict static nkl],
    const double norm1[restrict static nij], const double norm2[restrict static nij], const double norm3[restrict static nkl], const double norm4[restrict static nkl],
    const double rho12[restrict static nij], const double rho34[restrict static nkl],
    double cbatch[restrict static 1]);

void erd__xyz_to_ry_abcd(uint32_t nxyza, uint32_t nxyzb, uint32_t nxyzc, uint32_t nxyzd,
    uint32_t nrya, uint32_t nryb, uint32_t nryc, uint32_t nryd,
    uint32_t shella, uint32_t shellb, uint32_t shellc, uint32_t shelld,
    uint32_t nrowa[restrict static 1], uint32_t nrowb[restrict static 1], uint32_t nrowc[restrict static 1], uint32_t nrowd[restrict static 1],
    uint32_t nrota[restrict static 1], uint32_t nrotb[restrict static 1], uint32_t nrotc[restrict static 1], uint32_t nrotd[restrict static 1],
    uint32_t z00a[restrict static 1], uint32_t z00b[restrict static 1], uint32_t z00c[restrict static 1], uint32_t z00d[restrict static 1],
    uint32_t i0a1[restrict static 1], uint32_t i0b1[restrict static 1], uint32_t i0c1[restrict static 1], uint32_t i0d1[restrict static 1],
    uint32_t i0a2[restrict static 1], uint32_t i0b2[restrict static 1], uint32_t i0c2[restrict static 1], uint32_t i0d2[restrict static 1],
    uint32_t icore[restrict static 1], double zcore[restrict static 1]);

void erd__xyz_to_ry_matrix(uint32_t nxyz, uint32_t nrowmx, uint32_t l, double *temp, uint32_t *nrow, uint32_t *row, double *tmat);

void erd__spherical_transform(uint32_t m, uint32_t nrow, uint32_t nry, uint32_t lrow[restrict static nry], uint32_t row[restrict], const double rot[restrict], const double x[restrict], double y[restrict]);

void erd__hrr_step(uint32_t nabo, uint32_t mrowin,
    uint32_t mrowout, uint32_t nxyzx,
    uint32_t nxyza, uint32_t nxyzb, uint32_t nxyzao,
    uint32_t shellx, uint32_t shellp, uint32_t shellb,
    double abx, double aby, double abz,
    const uint32_t nrowin[restrict], const uint32_t rowin[restrict], const double win[restrict],
    uint32_t nrowout[restrict static nabo], uint32_t rowout[restrict], double wout[restrict]);

void erd__hrr_matrix(uint32_t nrothrr, uint32_t ncolhrr,
    uint32_t nxyzet, uint32_t nxyza, uint32_t nxyzp,
    uint32_t shella, uint32_t shellb, uint32_t shellp,
    uint32_t nabcoor, double abx, double aby, double abz,
    uint32_t *in1_ptr, uint32_t *in2_ptr,
    uint32_t *nrowout_ptr, uint32_t *nrow,
    uint32_t *row, double *t);

void erd__hrr_transform(uint32_t m, uint32_t nrow,
    uint32_t nxyza, uint32_t nxyzb,
    const uint32_t lrow[restrict static nxyzb], const uint32_t row[restrict],
    const double rot[restrict], const double x[restrict], double y[restrict]);

double erd__dsqmin_line_segments(double xp0, double yp0,
    double zp0, double xp1,
    double yp1, double zp1,
    double xq0, double yq0,
    double zq0, double xq1,
    double yq1, double zq1);

void erd__set_abcd(
    uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
    uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
    bool atomic,
    double x1, double y1, double z1,
    double x2, double y2, double z2, 
    double x3, double y3, double z3,
    double x4, double y4, double z4, bool spheric,
    uint32_t *restrict npgtoa_ptr, uint32_t *restrict npgtob_ptr, uint32_t *restrict npgtoc_ptr, uint32_t *restrict npgtod_ptr,
    uint32_t *restrict shella_ptr, uint32_t *restrict shellb_ptr, uint32_t *restrict shellc_ptr, uint32_t *restrict shelld_ptr,
    double *restrict xa_ptr, double *restrict ya_ptr, double *restrict za_ptr,
    double *restrict xb_ptr, double *restrict yb_ptr, double *restrict zb_ptr,
    double *restrict xc_ptr, double *restrict yc_ptr, double *restrict zc_ptr,
    double *restrict xd_ptr, double *restrict yd_ptr, double *restrict zd_ptr,
    uint32_t *restrict nxyza_ptr, uint32_t *restrict nxyzb_ptr, uint32_t *restrict nxyzc_ptr, uint32_t *restrict nxyzd_ptr,
    uint32_t *restrict nxyzet_ptr, uint32_t *restrict nxyzft_ptr,
    uint32_t *restrict nrya_ptr, uint32_t *restrict nryb_ptr, uint32_t *restrict nryc_ptr, uint32_t *restrict nryd_ptr,
    uint32_t *restrict nabcoor_ptr, uint32_t *restrict ncdcoor_ptr,
    uint32_t *restrict ncolhrr, uint32_t *restrict nrothrr,
    uint32_t *restrict nxyzhrr, bool *restrict empty, bool *restrict tr1234);  

void erd__normalize_cartesian(uint32_t m, uint32_t l, const double norm[restrict static l+1], double batch[restrict]);

void erd__cartesian_norms(uint32_t length, double norm[restrict static length+1]);

void erd__e0f0_pcgto_block(uint32_t nij, uint32_t nkl,
    uint32_t ngqp, uint32_t nmom,
    uint32_t nxyzet, uint32_t nxyzft,
    uint32_t nxyzp, uint32_t nxyzq,
    uint32_t shella, uint32_t shellp,
    uint32_t shellc, uint32_t shellq,
    double xa, double ya, double za,
    double xb, double yb, double zb,
    double xc, double yc, double zc,
    double xd, double yd, double zd,
    const double alphaa[restrict static nij], const double alphab[restrict static nij], const double alphac[restrict static nkl], const double alphad[restrict static nkl],
    const double cca[restrict static nij], const double ccb[restrict static nij], const double ccc[restrict static nkl], const double ccd[restrict static nkl],
    int **vrrtab, int ldvrrtab,
    const uint32_t prima[restrict static nij], const uint32_t primb[restrict static nij], const uint32_t primc[restrict static nkl], const uint32_t primd[restrict static nkl],
    const double norma[restrict static nij], const double normb[restrict static nij], const double normc[restrict static nkl], const double normd[restrict static nkl],
    const double rhoab[restrict static nij], const double rhocd[restrict static nkl],
    double batch[restrict static 1]);

void erd__2d_coefficients(uint32_t mij, uint32_t mkl, uint32_t ngqp,
    const double *restrict p, const double *restrict q,
    const double *restrict px, const double *restrict py, const double *restrict pz,
    const double *restrict qx, const double *restrict qy, const double *restrict qz,
    double  xa, double  ya, double  za,
    double  xc, double  yc, double  zc,
    const double *restrict pinvhf, const double *restrict qinvhf, const double *restrict pqpinv,
    const double *restrict rts, uint32_t case2d,
    double *restrict b00, double *restrict b01, double *restrict b10,
    double *restrict c00x, double *restrict c00y, double *restrict c00z,
    double *restrict d00x, double *restrict d00y, double *restrict d00z);

int erd__2d_pq_integrals(int shellp, int shellq,
                        int ngqexq, double *b00,
                        double *b01, double *b10, double *c00x,
                        double *c00y, double *c00z,
                        double *d00x, double *d00y,
                        double *d00z, int case2d,
                        double *int2dx, double *int2dy,
                        double *int2dz);

int erd__int2d_to_e0f0 (int shella, int shellp, int shellc, int shellq,
                        int ngqexq, int nxyzet, int nxyzft,
                        double *int2dx, double *int2dy, double *int2dz,
                        int **vrrtab, double *batch);

int erd__int2d_to_e000 (int shella, int shellp, int ngqp, int nexq, int ngqexq,
                        int nxyzet, int nxyzp,
                        double *int2dx, double *int2dy, double *int2dz,
                        double *temp1, double *temp2,
                        double *scale, double *batch);

void erd__rys_roots_weights(uint32_t nt, uint32_t ngqp, uint32_t nmom,
                            const double tval[restrict],
                            double rts[restrict], double wts[restrict]);

void erd__rys_1_roots_weights(int nt, const double tval[restrict], double rts[restrict], double wts[restrict]);

void erd__rys_2_roots_weights(int nt, const double tval[restrict], double rts[restrict], double wts[restrict]);

void erd__rys_3_roots_weights(int nt, const double tval[restrict], double rts[restrict], double wts[restrict]);

void erd__rys_4_roots_weights(int nt, const double tval[restrict], double rts[restrict], double wts[restrict]);

void erd__rys_5_roots_weights(int nt, const double tval[restrict], double rts[restrict], double wts[restrict]);

void erd__rys_x_roots_weights(int nt, int ntgqp, int ngqp,
    int nmom, const double tval[restrict],
    const double ryszero[restrict],
    double rts[restrict], double wts[restrict]);

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif


#endif /* __ERD_H__ */
