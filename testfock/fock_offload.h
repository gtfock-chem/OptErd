#ifndef __FOCK_INIT_H__
#define __FOCK_INIT_H__


#include "CInt.h"


#define TOLSRC 1e-10
#define ALLOC alloc_if(1) free_if(0)
#define REUSE alloc_if(0) free_if(0)

#define ALIGNED_8(size) ((((size) + 7)/8)*8)
#define NUM_F_COPIES_MIC(n) (((n)+3) >> 2)
#define MY_F_COPY_MIC(tid) ((tid) >> 2)

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
double *F1;
double *F2;
double *F3;
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

int MIC_init_devices ();

void MIC_copy_buffers (int num_devices, int nshells, int nnz,
                       int *shellptr, int *shellid, int *shellrid,
                       double *shellvalue, int *f_startind, int *rowpos,
                       int *colpos, int *rowptr, int *colptr);

void MIC_create_matrices (int num_devices,
                          double *D1, double *D2, double *D3,
                          double *F1_mic, double *F2_mic,
                          double *F3_mic,
                          int sizeD1, int sizeD2, int sizeD3,
                          int nthreads_mic);

void copy_double_array_CPU_to_MIC (int num_devices, double *A, int size);

void MIC_copy_D_matrices (int num_devices,
                          double *D1, double *D2, double *D3,
                          int sizeD1, int sizeD2, int sizeD3);

void MIC_reset_F1_matrix (int num_devices, int sizeD1, int nthreads_mic);

void MIC_reset_F2_matrix (int num_devices, int sizeD2, int nthreads_mic);

void MIC_reset_F3_matrix (int num_devices, int sizeD3, int nthreads_mic);

void MIC_reset_F_matrices (int num_devices,
                           int sizeD1, int sizeD2, int sizeD3,
                           int nthreads_mic);

void reset_F_matrices (int num_devices,
                       int sizeD1, int sizeD2, int sizeD3,
                       int nthreads, int nthreads_mic, int toOffload);

void compute_task (int num_devices,
                   BasisSet_t basis, ERD_t erd,
                   int *shellptr, double *shellvalue,
                   int *shellid, int *shellrid, int *f_startind,
                   int *rowpos, int *colpos, int *rowptr, int *colptr,
                   double tolscr2, int startrow, int startcol,
                   int startM, int endM, int startP, int endP,
                   double *D1, double *D2, double *D3,
                   int ldX1, int ldX2, int ldX3,
                   int sizeX1, int sizeX2, int sizeX3, double mic_fraction,
                   double *totalcalls, double *totalnintls, int toOffload);

void reduce_F_on_individual_devices (int num_devices,
                                     int sizeD1, int sizeD2, int sizeD3,
                                     int nthreads, int nthreads_mic,
                                     int toOffload);

void copy_F_MIC_to_CPU (int num_devices,
                        double *F1_mic, double *F2_mic, double *F3_mic,
                        int sizeD1, int sizeD2, int sizeD3, int *finish_tag);

void wait_for_MIC_to_CPU_copy (int num_devices, int *finish_tag);

void reduce_F_across_devices (int num_devices,
                              double *F1_mic, double *F2_mic, double *F3_mic,
                              int sizeD1, int sizeD2, int sizeD3);
#endif /* __FOCK_INIT_H__ */
