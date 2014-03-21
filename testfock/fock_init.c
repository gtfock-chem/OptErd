#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <assert.h>

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include "CInt.h"
#include "fock_init.h"


void schwartz_screening (BasisSet_t basis, int **shellptrOut,
                         int **shellidOut, int **shellridOut,
                         double **shellvalueOut, int *nnzOut)
{
    ERD_t erd;
    CInt_createERD (basis, &erd, 1);
    const int nshells = CInt_getNumShells (basis);

    double *vpairs = (double *) malloc (sizeof (double) * nshells * nshells);
    assert (vpairs != NULL);

    double allmax = 0.0;
    for (int M = 0; M < nshells; M++)
    {
        const int dimM = CInt_getShellDim (basis, M);
        for (int N = 0; N < nshells; N++)
        {
            const int dimN = CInt_getShellDim (basis, N);

            int nints;
            double *integrals;
            CInt_computeShellQuartet (basis, erd, 0, M, N, M, N, &integrals,
                                      &nints);

            double mvalue = 0.0;
            if (nints != 0)
            {
                for (int iM = 0; iM < dimM; iM++)
                {
                    for (int iN = 0; iN < dimN; iN++)
                    {
                        const int index =
                            iM * (dimN * dimM * dimN + dimN) +
                            iN * (dimM * dimN + 1);
                        if (mvalue < fabs (integrals[index]))
                        {
                            mvalue = fabs (integrals[index]);
                        }
                    }
                }
            }
            vpairs[M * nshells + N] = mvalue;
            if (mvalue > allmax)
            {
                allmax = mvalue;
            }
        }
    }

    // init shellptr
    int nnz = 0;
    const double eta = TOLSRC * TOLSRC / allmax;
    int *shellptr = (int *) malloc (sizeof (int) * (nshells + 1));
    assert (shellptr != NULL);
    memset (shellptr, 0, sizeof (int) * (nshells + 1));
    for (int M = 0; M < nshells; M++)
    {
        for (int N = 0; N < nshells; N++)
        {
            double mvalue = vpairs[M * nshells + N];
            if (mvalue > eta)
            {
                if (M > N && (M + N) % 2 == 1 || M < N && (M + N) % 2 == 0)
                {
                    continue;
                }
                else
                {
                    nnz++;
                }
            }
        }
        shellptr[M + 1] = nnz;
    }

    double *shellvalue = (double *) malloc (sizeof (double) * nnz);
    int *shellid = (int *) malloc (sizeof (int) * nnz);
    int *shellrid = (int *) malloc (sizeof (int) * nnz);
    assert ((shellvalue != NULL) && (shellid != NULL) && (shellrid != NULL));
    nnz = 0;
    for (int M = 0; M < nshells; M++)
    {
        for (int N = 0; N < nshells; N++)
        {
            const double mvalue = vpairs[M * nshells + N];
            if (mvalue > eta)
            {
                if (M > N && (M + N) % 2 == 1 || M < N && (M + N) % 2 == 0)
                {
                    continue;
                }
                if (M == N)
                {
                    shellvalue[nnz] = mvalue;                       
                }
                else
                {
                    shellvalue[nnz] = -mvalue;
                }
                shellid[nnz] = N;
                shellrid[nnz] = M;
                nnz++;
            }
        }
    }
    *nnzOut = nnz;
    free (vpairs);
    CInt_destroyERD (erd);

    *shellidOut = shellid;
    *shellridOut = shellrid;
    *shellptrOut = shellptr;
    *shellvalueOut = shellvalue;
}


static void compute_FD_ptr (int nshells, int *shellptr, int *shellid,
                            int *f_startind, int startM, int endM,
                            int *ptrrow, int *rowsize)
{
    int A;
    int B;
    int i;
    int start;
    int end;

    for (A = 0; A < nshells; A++)
    {
        ptrrow[A] = -1;
    }    
    // init row pointers
    for (A = startM; A <= endM; A++)
    {
        start = shellptr[A];
        end = shellptr[A + 1]; 
        for (i = start; i < end; i++)
        {
            B = shellid[i];
            ptrrow[B] = 1;
        }
    }
    
    *rowsize = 0;
    for (A = 0; A < nshells; A++)
    {
        if (ptrrow[A] == 1)
        {
            ptrrow[A] = *rowsize;           
            *rowsize += f_startind[A + 1] - f_startind[A];
        }
    }
}


void create_buffers (int nshells, int nnz,
                     int *shellptr, int *shellid, int *f_startind,
                     int startM, int endM, int startP, int endP,
                     int **rowposOut, int **colposOut,
                     int **rowptrOut, int **colptrOut,
                     int *rowfuncs, int *colfuncs,
                     int *rowsize, int *colsize)
{
    int *ptrrow;
    int *ptrcol;
    int *rowpos;
    int *colpos;
    int *rowptr;
    int *colptr;
    int j;
    int k;
    int sh;
    int count;
    
    ptrrow = (int *)malloc (sizeof(int) * nshells);
    ptrcol = (int *)malloc (sizeof(int) * nshells);
    assert (NULL != ptrrow && NULL != ptrcol);

    // compute rowptr/pos and colptr/pos
    rowpos = (int *)malloc (sizeof(int) * nshells);
    colpos = (int *)malloc (sizeof(int) * nshells);
    rowptr = (int *)malloc (sizeof(int) * nnz);
    colptr = (int *)malloc (sizeof(int) * nnz);
    assert (rowpos != NULL &&
            colpos != NULL &&
            rowptr != NULL &&
            colptr != NULL);
    
    
    compute_FD_ptr (nshells, shellptr, shellid, f_startind, startM, endM,
                    ptrrow, rowsize);
    count = 0;
    for (j = startM; j <= endM; j++)
    {
        rowpos[j] = ptrrow[j];
        for (k = shellptr[j]; k < shellptr[j+1]; k++)
        {
            sh = shellid[k];
            rowptr[count] = ptrrow[sh];
            count++;
        }
    }

    compute_FD_ptr (nshells, shellptr, shellid, f_startind,
                    startP, endP,
                    ptrcol, colsize);
    count= 0;
    for (j = startP; j <= endP; j++)
    {
       colpos[j] = ptrcol[j];
       for (k = shellptr[j]; k < shellptr[j+1]; k++)
       {
           sh = shellid[k];
           colptr[count] = ptrcol[sh];
           count++;
        }
    }
    
    free (ptrrow);
    free (ptrcol);

    *rowfuncs = f_startind[endM + 1] - f_startind[startM];
    *colfuncs = f_startind[endP + 1] - f_startind[startP];    
    *rowposOut = rowpos;
    *colposOut = colpos;
    *rowptrOut = rowptr;
    *colptrOut = colptr;

    printf ("FD size (%d %d %d %d)\n",
        *rowfuncs, *colfuncs, *rowsize, *colsize); 
}


int MIC_init_devices(int nthreads_mic)
{
    // Find the number of MIC cards available
    int num_devices = 0;
#ifdef __INTEL_OFFLOAD
    num_devices = _Offload_number_of_devices();

    if(num_devices == 0)
    {
        printf("No target devices available. Exiting\n");
        exit(0);
    }
    else
    {
        int mic_id;
        for(mic_id = 0; mic_id < num_devices; mic_id++)
        {
            printf("On CPU: nthreads_mic = %d\n", nthreads_mic);
#pragma offload target(mic:mic_id) in(nthreads_mic)
            {
                //printf("On MIC: nthreads_mic = %d\n", nthreads_mic);
                omp_set_num_threads (nthreads_mic);
#pragma omp parallel num_threads(nthreads_mic)
                {
                    int tid = omp_get_thread_num();
                    //if(tid == 0)
                    //    printf("Initialized openmp threads\n");
                }
            }
        }
        printf("Number of Target devices installed: %d\n\n",num_devices);
    }
#endif
    return num_devices;
}

void MIC_copy_buffers (int num_devices, int nshells, int nnz,
                       int *shellptr, int *shellid, int *shellrid, double *shellvalue,
                       int *f_startind,
                       int *rowpos, int *colpos,
                       int *rowptr, int *colptr)
{
#ifdef __INTEL_OFFLOAD
    int mic_id;

    for(mic_id = 0; mic_id < num_devices; mic_id++)
    {
#pragma offload_transfer target(mic:mic_id) \
        in(shellptr: length(nshells + 1) ALLOC) \
        in(shellid: length(nnz) ALLOC) \
        in(shellrid: length(nnz) ALLOC) \
        in(shellvalue: length(nnz) ALLOC) \
        in(f_startind: length(nshells + 1) ALLOC) \
        in(rowpos: length(nshells) ALLOC) \
        in(colpos: length(nshells) ALLOC) \
        in(rowptr: length(nnz) ALLOC) \
        in(colptr: length(nnz) ALLOC)
    }
#endif
}

void MIC_create_matrices(int num_devices,
                         double *D1, double *D2, double *D3,
                         double *F1_hetero, double *F2_hetero,
                         double *F3_hetero,
                         int sizeD1, int sizeD2, int sizeD3,
                         int nthreads_mic)
{
#ifdef __INTEL_OFFLOAD
    int mic_id;

    printf("num_f copies = %d\n", NUM_F_COPIES_MIC(nthreads_mic));
    for(mic_id = 0; mic_id < num_devices; mic_id++)
    {
#pragma offload target(mic: mic_id) \
        in(sizeD1, sizeD2, sizeD3, nthreads_mic) \
        nocopy(D1: length(sizeD1) ALLOC) \
        nocopy(D2: length(sizeD2) ALLOC) \
        nocopy(D3: length(sizeD3) ALLOC) \
        nocopy(F1_hetero: length(sizeD1) ALLOC) \
        nocopy(F2_hetero: length(sizeD2) ALLOC) \
        nocopy(F3_hetero: length(sizeD3) ALLOC) \
        nocopy(F1: length(0) alloc_if(0) free_if(0)) \
        nocopy(F2: length(0) alloc_if(0) free_if(0)) \
        nocopy(F3: length(0) alloc_if(0) free_if(0))
        {
            F1 = (double *)_mm_malloc(sizeD1 * NUM_F_COPIES_MIC(nthreads_mic) * sizeof(double), 64);
            F2 = (double *)_mm_malloc(sizeD2 * NUM_F_COPIES_MIC(nthreads_mic) * sizeof(double), 64);
            F3 = (double *)_mm_malloc(sizeD3 * NUM_F_COPIES_MIC(nthreads_mic) * sizeof(double), 64);
        }
    }
#endif
}


void copy_double_array_CPU_to_MIC(int num_devices,
                                  double *A, int size)
{
#ifdef __INTEL_OFFLOAD
    int mic_id;

    for(mic_id = 0; mic_id < num_devices; mic_id++)
    {
#pragma offload_transfer target(mic:mic_id) \
        in(A: length(size) REUSE)
    }
#endif
}

void MIC_copy_D_matrices(int num_devices,
                         double *D1, double *D2, double *D3,
                         int sizeD1, int sizeD2, int sizeD3)
{
    copy_double_array_CPU_to_MIC(num_devices, D1, sizeD1);
    copy_double_array_CPU_to_MIC(num_devices, D2, sizeD2);
    copy_double_array_CPU_to_MIC(num_devices, D3, sizeD3);

    printf("Sent D1, D2, D3 to MIC cards\n");

}

