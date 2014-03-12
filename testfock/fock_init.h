#ifndef __FOCK_INIT_H__
#define __FOCK_INIT_H__


#include "CInt.h"


#define TOLSRC 1e-10
#define ALLOC alloc_if(1) free_if(0)
#define REUSE alloc_if(0) free_if(0)

#define ALIGNED_8(size) (((size + 7)/8)*8)

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
ERD_t erd_mic;
double *F1_mic;
double *F2_mic;
double *F3_mic;
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

void schwartz_screening (BasisSet_t basis, int **shellptrOut,
                         int **shellidOut, int **shellridOut,
                         double **shellvalueOut, int *nnzOut);

void create_buffers (int nshells, int nnz,
                     int *shellptr, int *shellid, int *f_startind,
                     int startM, int endM, int startP, int endP,
                     int **rowposOut, int **colposOut,
                     int **rowptrOut, int **colptrOut,
                     int *rowfuncs, int *colfuncs,
                     int *rowsize, int *colsize);

int MIC_init_devices();

void MIC_copy_buffers (int num_devices, int nshells, int nnz,
                       int *shellptr, int *shellid, int *shellrid, double *shellvalue,
                       int *f_startind, int startM, int endM, int startP, int endP,
                       int *rowpos, int *colpos,
                       int *rowptr, int *colptr,
                       int rowfuncs, int colfuncs,
                       int rowsize, int colsize);

void MIC_create_matrices(int num_devices,
                         double *D1, double *D2, double *D3,
                         double *F1_hetero, double *F2_hetero,
                         double *F3_hetero,
                         int sizeD1, int sizeD2, int sizeD3,
                         int nthreads_mic);

void MIC_copy_D_matrices(int num_devices,
                         double *D1, double *D2, double *D3,
                         int sizeD1, int sizeD2, int sizeD3);

void MIC_reset_F_matrices(int num_devices,
                          int sizeD1, int sizeD2, int sizeD3,
                          int nthreads_mic);

void compute_task (BasisSet_t basis, ERD_t erd,
                   int *shellptr, double *shellvalue,
                   int *shellid, int *shellrid, int *f_startind,
                   int *rowpos, int *colpos, int *rowptr, int *colptr,
                   double tolscr2, int startrow, int startcol,
                   int startM, int endM, int startP, int endP,
                   double *D1, double *D2, double *D3,
                   double *F1, double *F2, double *F3,
                   int ldX1, int ldX2, int ldX3,
                   int sizeX1, int sizeX2, int sizeX3,
                   double *totalcalls, double *totalnintls);

void compute_task_hetero (int num_devices,
                   BasisSet_t basis, ERD_t erd,
                   int *shellptr, double *shellvalue,
                   int *shellid, int *shellrid, int *f_startind,
                   int *rowpos, int *colpos, int *rowptr, int *colptr,
                   double tolscr2, int startrow, int startcol,
                   int startM, int endM, int startP, int endP,
                   double *D1, double *D2, double *D3,
                   double *F1, double *F2, double *F3,
                   double *F1_hetero, double *F2_hetero, double *F3_hetero,
                   int ldX1, int ldX2, int ldX3,
                   int sizeX1, int sizeX2, int sizeX3, double mic_fraction,
                   double *totalcalls, double *totalnintls);

void reduce_F(int num_devices,
              double *F1, double *F2, double *F3,
              double *F1_hetero, double *F2_hetero, double *F3_hetero,
              int sizeD1, int sizeD2, int sizeD3,
              int nthreads, int nthreads_mic);

void reduce_F_CPU_MIC(int num_devices,
                      double *F1, double *F2, double *F3,
                      double *F1_hetero, double *F2_hetero, double *F3_hetero,
                      int sizeD1, int sizeD2, int sizeD3);
#endif /* __FOCK_INIT_H__ */
