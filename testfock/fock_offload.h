#ifndef __FOCK_OFFLOAD_H__
#define __FOCK_OFFLOAD_H__


#include "CInt.h"

#define ALLOC alloc_if(1) free_if(0)
#define REUSE alloc_if(0) free_if(0)
#define ONCE  alloc_if(1) free_if(1)
#define FREE  alloc_if(0) free_if(1)

#define ALIGNED_8(size) ((((size) + 7)/8)*8)
#define NUM_F_COPIES_MIC(n) (((n)+3) >> 2)
#define MY_F_COPY_MIC(tid) ((tid) >> 2)


typedef struct _pfock_mic_t
{
    double **F1;
    double **F2;
    double **F3;
    int sizeD1;
    int sizeD2;
    int sizeD3;
    int num_F;
    int nthreads;
} pfock_mic_t;


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
extern pfock_mic_t pfock_mic;
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif


void offload_reset_F (int num_devices);


void compute_task (int num_devices,
                   BasisSet_t basis, ERD_t erd,
                   int *shellptr, double *shellvalue,
                   int *shellid, int *shellrid, int *f_startind,
                   int *rowpos, int *colpos, int *rowptr, int *colptr,
                   double tolscr2, int startrow, int startcol,
                   int startM, int endM, int startP, int endP,
                   double *D1, double *D2, double *D3,
                   double *F1, double *F2, double *F3,
                   int ldX1, int ldX2, int ldX3,
                   int sizeX1, int sizeX2, int sizeX3, double mic_fraction,
                   double *totalcalls, double *totalnintls,
                   int toOffload);


void offload_reduce_mic (int num_devices,
                         double *F1_offload,
                         double *F2_offload,
                         double *F3_offload,
                         int sizeD1, int sizeD2, int sizeD3);


void offload_reduce (int num_devices,
                     double *F1, double *F2, double *F3,
                     double *F1_offload, double *F2_offload, double *F3_offload,
                     int sizeD1, int sizeD2, int sizeD3);


void offload_init (int num_devices, int nshells, int nnz,
                   int *shellptr, int *shellid, 
                   int *shellrid, double *shellvalue,
                   int *f_startind,
                   int *rowpos, int *colpos,
                   int *rowptr, int *colptr,
                   double *D1, double *D2, double *D3,
                   double *F1_offload, double *F2_offload, double *F3_offload,
                   int sizeD1, int sizeD2, int sizeD3,
                   int nthreads_mic);


void offload_copy_D (int num_devices, double *D, int sizeD);


void offload_wait_mic (int num_devices);


#endif /* __FOCK_OFFLOAD_H__ */
