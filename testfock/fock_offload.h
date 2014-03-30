#ifndef __FOCK_OFFLOAD_H__
#define __FOCK_OFFLOAD_H__


#include "CInt.h"


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



#pragma offload_attribute(push, target(mic))

extern pfock_mic_t pfock_mic;

#pragma offload_attribute(pop)


void offload_fock_task (int num_devices,
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
                        double *totalcalls, double *totalnintls);

void offload_init (int nshells, int nnz,
                   int *shellptr, int *shellid, 
                   int *shellrid, double *shellvalue,
                   int *f_startind,
                   int *rowpos, int *colpos,
                   int *rowptr, int *colptr,
                   double *D1, double *D2, double *D3,
                   double *VD1, double *VD2, double *VD3,
                   double **_F1_offload,
                   double **_F2_offload,
                   double **_F3_offload,
                   int sizeD1, int sizeD2, int sizeD3,
                   int *_mic_numdevs,
                   int *_nthreads_mic);
                   
void offload_deinit (int mic_numdevs,
                     int *shellptr, int *shellid, 
                     int *shellrid, double *shellvalue,
                     int *f_startind,
                     int *rowpos, int *colpos,
                     int *rowptr, int *colptr,
                     double *D1, double *D2, double *D3,
                     double *VD1, double *VD2, double *VD3,
                     double *F1_offload,
                     double *F2_offload,
                     double *F3_offload);

void offload_reset_F (int num_devices);

void offload_reduce_mic (int num_devices,
                         double *F1_offload,
                         double *F2_offload,
                         double *F3_offload,
                         int sizeD1, int sizeD2, int sizeD3);


void offload_reduce (int num_devices,
                     double *F1, double *F2, double *F3,
                     double *F1_offload, double *F2_offload, double *F3_offload,
                     int sizeD1, int sizeD2, int sizeD3);

void offload_copy_D (int num_devices, double *D, int sizeD);


void offload_wait_mic (int num_devices);


#endif /* __FOCK_OFFLOAD_H__ */
