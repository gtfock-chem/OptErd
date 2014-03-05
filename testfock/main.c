#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>
#include <math.h>
#include <sys/time.h>

#include "erd_profile.h"
#include "fock_init.h"


int main (int argc, char **argv)
{
    int nnz;
    int *shellptr;
    int *shellid;
    int *shellrid;
    double *shellvalue;
    int *rowpos;
    int *colpos;
    int *rowptr;
    int *colptr;
    int rowsize;
    int colsize;
    int rowfuncs;
    int colfuncs;
    int sizeD1;
    int sizeD2;
    int sizeD3;
    int nshells;
    int *f_startind;
    int startM;
    int endM;
    int startP;
    int endP;
    int sizetask;
    int nfuncs;
    int i;
    int j;
    struct timeval tv1, tv2;
    double timepass;
    
    if (argc != 9)
    {
        printf ("Usage: %s <basis set> <xyz> <nthreads>"
                " <startM> <endN> <startP> <endP> <ntasks>\n",
                argv[0]);
        return -1;
    }

    const int nthreads = atoi (argv[3]);
    omp_set_num_threads (nthreads);

    // load basis set
    BasisSet_t basis;
    CInt_createBasisSet (&basis);
    CInt_loadBasisSet (basis, argv[1], argv[2]);
    printf ("Molecule info:\n");
    printf ("  #Atoms\t= %d\n", CInt_getNumAtoms (basis));
    printf ("  #Shells\t= %d\n", CInt_getNumShells (basis));
    printf ("  #Funcs\t= %d\n", CInt_getNumFuncs (basis));
    printf ("  #OccOrb\t= %d\n", CInt_getNumOccOrb (basis));
    printf ("  nthreads\t= %d\n", nthreads);
    nshells = CInt_getNumShells (basis);
    nfuncs = CInt_getNumFuncs (basis);
    // functions starting positions of shells
    f_startind = (int *)malloc (sizeof(int) * (nshells + 1));
    assert (NULL != f_startind);
    for (i = 0; i < nshells; i++)
    {
        f_startind[i] = CInt_getFuncStartInd (basis, i);
    }
    f_startind[nshells] = nfuncs;
    
    startM = atoi(argv[4]);
    endM = atoi(argv[5]);
    startP = atoi(argv[6]);
    endP = atoi(argv[7]);
    sizetask = atoi (argv[8]);
    assert (endM <= nshells - 1);
    assert (endP <= nshells - 1);
    assert (sizetask > 0);
    assert (startM >= 0 && startM <=endM);
    assert (startP >= 0 && startP <=endP);
    assert (startM + sizetask <= endM + 1);
    assert (startP + sizetask <= endP + 1);
    printf ("\ncreate FD for (%d:%d,0:%d|%d:%d,0:%d)\n",
        startM, endM, nshells - 1,
        startP, endP, nshells - 1);
    printf ("compute (%d:%d,0:%d|%d:%d,0:%d)\n",
        startM, startM + sizetask - 1, nshells - 1,
        startP, startP + sizetask - 1, nshells - 1);

    // initialization
    schwartz_screening (basis, &shellptr, &shellid,
                        &shellrid, &shellvalue, &nnz);  
    create_buffers (nshells, nnz, shellptr, shellid,
                    f_startind, startM, endM, startP, endP,
                    &rowpos, &colpos, &rowptr, &colptr,
                    &rowfuncs, &colfuncs, &rowsize, &colsize);
    
    // malloc D and F
    double totalFDsize = 0.0;
    sizeD1 = rowfuncs * rowsize;
    sizeD2 = colfuncs * colsize;
    sizeD3 = rowsize * colsize;   
    double *D1 = (double *)malloc (sizeof(double) * sizeD1);
    double *D2 = (double *)malloc (sizeof(double) * sizeD2); 
    double *D3 = (double *)malloc (sizeof(double) * sizeD3);
    totalFDsize += 1.0 * sizeof(double) * (sizeD1 + sizeD2 + sizeD3);
    double *F1 = (double *)malloc (sizeof(double) * sizeD1 * nthreads);
    double *F2 = (double *)malloc (sizeof(double) * sizeD2 * nthreads); 
    double *F3 = (double *)malloc (sizeof(double) * sizeD3 * nthreads);
    totalFDsize += 1.0 * sizeof(double) * (sizeD1 + sizeD2 + sizeD3) * nthreads;
    assert (D1 != NULL &&
            D2 != NULL &&
            D3 != NULL &&
            F1 != NULL &&
            F2 != NULL &&
            F3 != NULL);
    printf ("use %.3lf MB\n", totalFDsize/1024.0/1024.0);

    // init D
    #pragma omp parallel for
    for (j = 0; j < sizeD1; j++)
    {
        D1[j] = 1.0;
    }
    #pragma omp parallel for
    for (j = 0; j < sizeD2; j++)
    {
        D2[j] = 1.0;
    }
    #pragma omp parallel for 
    for (j = 0; j < sizeD3; j++)
    {
        D3[j] = 1.0;
    }
    
    ERD_t erd;
    double tolscr = TOLSRC;
    double tolscr2 = tolscr * tolscr;
    CInt_createERD (basis, &erd, nthreads);

    double* totalcalls = (double *) malloc(sizeof (double) * nthreads * 64);
    assert(totalcalls != NULL);
    double* totalnintls = (double *) malloc(sizeof (double) * nthreads * 64);
    assert(totalnintls != NULL);

    #pragma omp parallel for
    for (int i = 0; i < nthreads; i++) {
        totalcalls[i * 64] = 0.0;
        totalnintls[i * 64] = 0.0;
    }
    
    gettimeofday (&tv1, NULL);   
    /************************************************************/
    // init F
    #pragma omp parallel for
    for (j = 0; j < sizeD1; j++)
    {
        F1[j] = 0.0;
    }
    #pragma omp parallel for
    for (j = 0; j < sizeD2; j++)
    {
        F2[j] = 0.0;
    }
    #pragma omp parallel for 
    for (j = 0; j < sizeD3; j++)
    {
        F3[j] = 0.0;
    }

    // main computation
    compute_task (basis, erd, shellptr,
                  shellvalue, shellid, shellrid,
                  f_startind, rowpos, colpos, rowptr, colptr,
                  tolscr2, startM, startP,
                  startM, startM + sizetask - 1,
                  startP, startP + sizetask - 1,
                  D1, D2, D3, F1, F2, F3,
                  rowsize, colsize, colsize, sizeD1, sizeD2, sizeD3,
                  totalcalls, totalnintls);

    // reduction
    #pragma omp parallel for
    for (j = 0; j < sizeD1; j++)
    {
        for (i = 1; i < nthreads; i++)
        {
            F1[j + 0 * sizeD1] += F1[j + i * sizeD1];
        }
    }
    #pragma omp parallel for
    for (j = 0; j < sizeD2; j++)
    {
        for (i = 1; i < nthreads; i++)
        {
            F2[j + 0 * sizeD2] += F2[j + i * sizeD2];
        }
    }
    #pragma omp parallel for 
    for (j = 0; j < sizeD3; j++)
    {
        for (i = 1; i < nthreads; i++)
        {
            F3[j + 0 * sizeD3] += F3[j + i * sizeD3];
        }
    }

    /***********************************************************/
    gettimeofday (&tv2, NULL);
    timepass = (tv2.tv_sec - tv1.tv_sec) +
        (tv2.tv_usec - tv1.tv_usec)/1e6;
    printf ("takes %.3lf secs\n", timepass);

    for (int i = 1; i < nthreads; i++)
    {
        totalcalls[0 * 64] = totalcalls[0 * 64] + totalcalls[i * 64];
        totalnintls[0 * 64] = totalnintls[0 * 64] + totalnintls[i * 64];
    }
    printf ("%.4le usq, %.4le uints\n", totalcalls[0], totalnintls[0]);
    
    // use 1 if thread timing is not required
    // erd_print_profile (1);

    CInt_destroyERD (erd);
    free (rowpos);
    free (colpos);
    free (rowptr);
    free (colptr);
    free (shellptr);
    free (shellid);
    free (shellvalue);
    free (shellrid);
    free (f_startind);
    free (totalnintls);
    free (totalcalls);

    CInt_destroyBasisSet (basis);

    return 0;
}
