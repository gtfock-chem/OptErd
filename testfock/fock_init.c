#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <assert.h>

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
