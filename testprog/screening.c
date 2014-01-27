#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <assert.h>

#include "screening.h"


void schwartz_screening (BasisSet_t basis, int **shellptr,
                         int **shellid, int **shellrid,
                         double **shellvalue)
{
    ERD_t erd;
    int nshells;
    int M;
    int N;
    int dimM;
    int dimN;
    int nints;
    int iM;
    int iN;
    double mvalue;
    double allmax;
    int index;
    double *integrals;
    int *_shellptr;
    double eta;
    int nnz;
    double *_shellvalue;
    int *_shellid;
    int *_shellrid;
    double *vpairs;
    
    CInt_createERD (basis, &erd);
    nshells = CInt_getNumShells (basis);
    
    vpairs = (double *)malloc (sizeof(double) * nshells * nshells);
    assert (vpairs != NULL);

    allmax = 0.0;
    for (M = 0; M < nshells; M++)
    {
        dimM = CInt_getShellDim (basis, M);
        for (N = 0; N < nshells; N++)
        {
            dimN = CInt_getShellDim (basis, N);
            CInt_computeShellQuartet (basis, erd, M, N, M, N, &integrals, &nints);            
            mvalue = 0.0;
            if (nints != 0)
            {
                for (iM = 0; iM < dimM; iM++)
                {
                    for (iN = 0; iN < dimN; iN++)
                    {
                        index = iM * (dimN*dimM*dimN+dimN) + iN * (dimM*dimN+1);
                        if (mvalue < fabs (integrals[index]))
                        {
                            mvalue = fabs (integrals[index]);                    
                        }
                    }
                }
            }
            vpairs[M * nshells  + N] = mvalue;
            if (mvalue > allmax)
            {
                allmax = mvalue;
            }
        }
    }

    // init shellptr
    nnz = 0;
    eta = TOLSRC*TOLSRC/allmax;
    _shellptr = (int *)malloc (sizeof(int) * (nshells + 1));
    assert (_shellptr != NULL);
    memset (_shellptr, 0, sizeof(int) * (nshells + 1));
    for (M = 0; M < nshells; M++)
    {
        for (N = 0; N < nshells; N++)
        {
            mvalue = vpairs[M * nshells  + N];
            if (mvalue > eta)
            {
                nnz++;
            }
        }
        _shellptr[M + 1] = nnz;
    }
 
    _shellvalue  = (double *)malloc (sizeof(double) * nnz);
    _shellid  = (int *)malloc (sizeof(int) * nnz);
    _shellrid  = (int *)malloc (sizeof(int) * nnz);
    assert (_shellvalue != NULL &&
            _shellid != NULL &&
             _shellrid != NULL);    
    nnz = 0;   
    for (M = 0; M < nshells; M++)
    {
        for (N = 0; N < nshells; N++)
        {
            mvalue = vpairs[M * nshells  + N];
            if (mvalue > eta)
            {
                _shellvalue[nnz] = mvalue;                       
                _shellid[nnz] = N;
                _shellrid[nnz] = M;
                nnz++;
            }
        }
    }

    free (vpairs);
    CInt_destroyERD (erd);

    *shellid = _shellid;
    *shellrid = _shellrid;
    *shellptr = _shellptr;
    *shellvalue = _shellvalue;
}
