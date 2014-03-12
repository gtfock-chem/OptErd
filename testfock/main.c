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
    int toOffload = 0;
    int num_devices = 0;
    struct timeval tv1, tv2;
    double timepass;
    
    if (argc != 12)
    {
        printf ("Usage: %s <basis set> <xyz> <nthreads>"
                " <startM> <endN> <startP> <endP> <ntasks> <toOffload> <nthreads_mic>\n",
                argv[0]);
        return -1;
    }

    const int nthreads = atoi (argv[3]);
    omp_set_num_threads (nthreads);
    const int nthreads_mic = atoi (argv[10]);

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
    toOffload = atoi(argv[9]);
    double mic_fraction = atof(argv[11]);
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

#ifdef __INTEL_OFFLOAD
    if(toOffload)
    {
        printf("Before MIC_init_devices\n");
        num_devices = MIC_init_devices(nthreads_mic);
    }
    printf("num_devices = %d\n", num_devices);
#endif

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
    int sizeD1_aligned = ALIGNED_8(sizeD1);
    sizeD2 = colfuncs * colsize;
    int sizeD2_aligned = ALIGNED_8(sizeD2);
    sizeD3 = rowsize * colsize;   
    int sizeD3_aligned = ALIGNED_8(sizeD3);
    double *D1 = (double *)malloc (sizeof(double) * sizeD1_aligned);
    double *D2 = (double *)malloc (sizeof(double) * sizeD2_aligned); 
    double *D3 = (double *)malloc (sizeof(double) * sizeD3_aligned);
    totalFDsize += 1.0 * sizeof(double) * (sizeD1_aligned + sizeD2_aligned + sizeD3_aligned);
    double *F1 = (double *)malloc (sizeof(double) * sizeD1_aligned * nthreads);
    double *F2 = (double *)malloc (sizeof(double) * sizeD2_aligned * nthreads); 
    double *F3 = (double *)malloc (sizeof(double) * sizeD3_aligned * nthreads);
    totalFDsize += 1.0 * sizeof(double) * (sizeD1_aligned + sizeD2_aligned + sizeD3_aligned) * nthreads;
    assert (D1 != NULL &&
            D2 != NULL &&
            D3 != NULL &&
            F1 != NULL &&
            F2 != NULL &&
            F3 != NULL);
    printf ("use %.3lf MB\n", totalFDsize/1024.0/1024.0);

    double *F1_hetero = NULL;
    double *F2_hetero = NULL;
    double *F3_hetero = NULL;
#ifdef __INTEL_OFFLOAD
    if(toOffload)
    {
        F1_hetero = (double *)_mm_malloc (sizeof(double) * sizeD1_aligned * (num_devices + 1), 64);
        F2_hetero = (double *)_mm_malloc (sizeof(double) * sizeD2_aligned * (num_devices + 1), 64);
        F3_hetero = (double *)_mm_malloc (sizeof(double) * sizeD3_aligned * (num_devices + 1), 64);

        memset(F1_hetero, 0, sizeof(double) * sizeD1_aligned * (num_devices + 1));
        memset(F2_hetero, 0, sizeof(double) * sizeD2_aligned * (num_devices + 1));
        memset(F3_hetero, 0, sizeof(double) * sizeD3_aligned * (num_devices + 1));
        MIC_copy_buffers(2, nshells, nnz, shellptr, shellid, shellrid, shellvalue,
                f_startind, startM, endM, startP, endP,
                rowpos, colpos, rowptr, colptr,
                rowfuncs, colfuncs, rowsize, colsize);
        printf("Copied buffers to MIC\n");
        MIC_create_matrices(num_devices, D1, D2, D3, F1_hetero, F2_hetero, F3_hetero, sizeD1_aligned, sizeD2_aligned, sizeD3_aligned, nthreads_mic);
    }
#endif

    // init D
    #pragma omp parallel for
    for (j = 0; j < sizeD1_aligned; j++)
    {
        D1[j] = 1.0;
    }
    #pragma omp parallel for
    for (j = 0; j < sizeD2_aligned; j++)
    {
        D2[j] = 1.0;
    }
    #pragma omp parallel for 
    for (j = 0; j < sizeD3_aligned; j++)
    {
        D3[j] = 1.0;
    }
  
#ifdef __INTEL_OFFLOAD
    if(toOffload)
    {
        MIC_copy_D_matrices(num_devices, D1, D2, D3, sizeD1_aligned, sizeD2_aligned, sizeD3_aligned);
    }
#endif
    ERD_t erd;
    double tolscr = TOLSRC;
    double tolscr2 = tolscr * tolscr;
    CInt_createERD (basis, &erd, nthreads);

#ifdef __INTEL_OFFLOAD
    if(toOffload)
    {
        int mic_id;
        for (mic_id = 0; mic_id < num_devices; mic_id++)
        {
#pragma offload target(mic:mic_id) nocopy(basis_mic, erd_mic) in(nthreads_mic)
            {
                CInt_createERD (basis_mic, &erd_mic, nthreads_mic);
            }
        }
    }
#endif
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
    for (j = 0; j < sizeD1_aligned * nthreads; j++)
    {
        F1[j] = 0.0;
    }
    #pragma omp parallel for
    for (j = 0; j < sizeD2_aligned * nthreads; j++)
    {
        F2[j] = 0.0;
    }
    #pragma omp parallel for 
    for (j = 0; j < sizeD3_aligned * nthreads; j++)
    {
        F3[j] = 0.0;
    }

#ifdef __INTEL_OFFLOAD
    if(toOffload)
    {
        MIC_reset_F_matrices(num_devices, sizeD1_aligned, sizeD2_aligned, sizeD3_aligned, nthreads_mic);
    }
#endif

    long start, end;
    start = __rdtsc();
    // main computation
    compute_task_hetero (num_devices, basis, erd, shellptr,
                  shellvalue, shellid, shellrid,
                  f_startind, rowpos, colpos, rowptr, colptr,
                  tolscr2, startM, startP,
                  startM, startM + sizetask - 1,
                  startP, startP + sizetask - 1,
                  D1, D2, D3, F1, F2, F3, F1_hetero, F2_hetero, F3_hetero,
                  rowsize, colsize, colsize, sizeD1, sizeD2, sizeD3, mic_fraction,
                  totalcalls, totalnintls);
    end = __rdtsc();
    printf("Compute cycles = %lld\n", end - start);

    start = __rdtsc();
    reduce_F(num_devices, F1, F2, F3, F1_hetero, F2_hetero, F3_hetero, sizeD1, sizeD2, sizeD3,
             nthreads, nthreads_mic);
    end = __rdtsc();
    printf("Reduce cycles = %lld\n", end - start);

    start = __rdtsc();
    reduce_F_CPU_MIC(num_devices, F1, F2, F3, F1_hetero, F2_hetero, F3_hetero, sizeD1, sizeD2, sizeD3);
    end = __rdtsc();
    printf("CPU-MIC Reduce cycles = %lld\n", end - start);
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
  
#if 1
    // Print the output for comparison
    printf("F1:\n");
    for(j = 0; j < sizeD1; j++)
        printf("%E\n", F1[j]);
    printf("F2:\n");
    for(j = 0; j < sizeD2; j++)
        printf("%E\n", F2[j]);
    printf("F3:\n");
    for(j = 0; j < sizeD3; j++)
        printf("%E\n", F3[j]);
#endif
            
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
