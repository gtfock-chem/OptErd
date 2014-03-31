#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>
#include <math.h>
#include <sys/time.h>

#include "fock_init.h"
#include "fock_offload.h"


#define ALIGNED_8(size) ((((size) + 7)/8)*8)


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
    int task_startM;
    int task_endM;
    int task_startP;
    int task_endP;
    int nfuncs;
    int toOffload = 0;
    int mic_numdevs = 0;
    struct timeval tv1, tv2;
    double timepass;
    int nthreads_mic;

    if (argc != 14)
    {
        printf ("Usage: %s <basis set> <xyz> <nthreads>"
                " <startM> <endM> <startP> <endP>"
                " <task_startM> <task_endM> <task_startP> <task_endP>"
                " <toOffload> <mic_fraction>\n",
                argv[0]);
        return -1;
    }

    const int nthreads = atoi (argv[3]);
    omp_set_num_threads (nthreads);
    startM = atoi (argv[4]);
    endM = atoi (argv[5]);
    startP = atoi (argv[6]);
    endP = atoi (argv[7]);
    task_startM = atoi (argv[8]);
    task_endM = atoi (argv[9]);
    task_startP = atoi (argv[10]);
    task_endP = atoi (argv[11]);
    toOffload = atoi (argv[12]);
    double mic_fraction = atof (argv[13]);

    // load basis set and create ERD_t
    BasisSet_t basis;
    ERD_t erd;
    if (toOffload == 1)
    {
        CInt_offload_createBasisSet (&basis);
        CInt_offload_loadBasisSet (basis, argv[1], argv[2]);
        
    }
    else
    {
        CInt_createBasisSet (&basis);
        CInt_loadBasisSet (basis, argv[1], argv[2]);
    }
    printf ("Molecule info:\n");
    printf ("  #Atoms\t= %d\n", CInt_getNumAtoms (basis));
    printf ("  #Shells\t= %d\n", CInt_getNumShells (basis));
    printf ("  #Funcs\t= %d\n", CInt_getNumFuncs (basis));
    printf ("  #OccOrb\t= %d\n", CInt_getNumOccOrb (basis));
    printf ("  nthreads\t= %d\n", nthreads);

    nshells = CInt_getNumShells (basis);
    nfuncs = CInt_getNumFuncs (basis);
    // functions starting positions of shells
    f_startind = (int *) malloc (sizeof (int) * (nshells + 1));
    assert (NULL != f_startind);
    for (int i = 0; i < nshells; i++)
    {
        f_startind[i] = CInt_getFuncStartInd (basis, i);
    }
    f_startind[nshells] = nfuncs;

    assert (endM <= nshells - 1);
    assert (endP <= nshells - 1);
    assert (startM >= 0 && startM <= endM);
    assert (startP >= 0 && startP <= endP);
    assert (task_startM >= 0 && task_startM <= task_endM);
    assert (task_startP >= 0 && task_startP <= task_endP);
    assert (task_startM >= startM);
    assert (task_startP >= startP);
    assert (task_endM <= endM);
    assert (task_endP <= endP);   
    printf ("\ncreate FD for (%d:%d,0:%d|%d:%d,0:%d)\n",
            startM, endM, nshells - 1, startP, endP, nshells - 1);
    printf ("compute (%d:%d,0:%d|%d:%d,0:%d)\n",
            task_startM, task_endM, nshells - 1,
            task_startP, task_endP, nshells - 1);

    // fock initialization
    schwartz_screening (basis, &shellptr, &shellid,
                        &shellrid, &shellvalue, &nnz);
    create_buffers (nshells, nnz, shellptr, shellid,
                    f_startind, startM, endM, startP, endP,
                    &rowpos, &colpos, &rowptr, &colptr,
                    &rowfuncs, &colfuncs, &rowsize, &colsize);

    // malloc D and F
    double totalFDsize = 0.0;
    sizeD1 = rowfuncs * rowsize;
    int sizeD1_aligned = ALIGNED_8 (sizeD1);
    sizeD2 = colfuncs * colsize;
    int sizeD2_aligned = ALIGNED_8 (sizeD2);
    sizeD3 = rowsize * colsize;
    int sizeD3_aligned = ALIGNED_8 (sizeD3);
    double *D1 = (double *) _mm_malloc (sizeof (double) * sizeD1_aligned, 64);
    double *D2 = (double *) _mm_malloc (sizeof (double) * sizeD2_aligned, 64);
    double *D3 = (double *) _mm_malloc (sizeof (double) * sizeD3_aligned, 64); 
    double *VD1 = (double *) _mm_malloc (sizeof (double) * sizeD1_aligned, 64);
    double *VD2 = (double *) _mm_malloc (sizeof (double) * sizeD2_aligned, 64);
    double *VD3 = (double *) _mm_malloc (sizeof (double) * sizeD3_aligned, 64);
    double *F1 = (double *) _mm_malloc (sizeof (double) * sizeD1_aligned * nthreads, 64);
    double *F2 = (double *) _mm_malloc (sizeof (double) * sizeD2_aligned * nthreads, 64);
    double *F3 = (double *) _mm_malloc (sizeof (double) * sizeD3_aligned * nthreads, 64);   
    totalFDsize += 2.0 * sizeof (double) *
        (sizeD1_aligned + sizeD2_aligned + sizeD3_aligned);
    totalFDsize += 1.0 * sizeof (double) *
        (sizeD1_aligned + sizeD2_aligned + sizeD3_aligned) * nthreads;
    assert (D1 != NULL &&
            D2 != NULL &&
            D3 != NULL &&
            VD1 != NULL &&
            VD2 != NULL &&
            VD3 != NULL &&
            F1 != NULL &&
            F2 != NULL &&
            F3 != NULL);
    printf ("use %.3lf MB\n", totalFDsize / 1024.0 / 1024.0);

    double *F1_offload = NULL;
    double *F2_offload = NULL;
    double *F3_offload = NULL;
    if (toOffload == 1)
    {
        offload_init (nshells, nnz, shellptr, shellid, shellrid,
                      shellvalue, f_startind, rowpos, colpos, rowptr, colptr,
                      D1, D2, D3,
                      VD1, VD2, VD3,
                      &F1_offload, &F2_offload, &F3_offload,
                      sizeD1_aligned, sizeD2_aligned,
                      sizeD3_aligned,
                      &mic_numdevs, &nthreads_mic);
        CInt_offload_createERD (basis, &erd, nthreads, nthreads_mic);
        printf ("mic_numdevs %d, num_threads_mic %d\n", mic_numdevs, nthreads_mic);
    }
    else
    {
        CInt_createERD (basis, &erd, nthreads);
    }
    
    // init D
    #pragma omp parallel for
    for (int j = 0; j < sizeD1_aligned; j++)
    {
        D1[j] = 1.0;
    }
    #pragma omp parallel for
    for (int j = 0; j < sizeD2_aligned; j++)
    {
        D2[j] = 1.0;
    }
    #pragma omp parallel for
    for (int j = 0; j < sizeD3_aligned; j++)
    {
        D3[j] = 1.0;
    }

    if (toOffload == 1)
    {
        offload_copy_D (mic_numdevs, D1, sizeD1_aligned);
        offload_copy_D (mic_numdevs, D2, sizeD2_aligned);
        offload_copy_D (mic_numdevs, D3, sizeD3_aligned);
    }
    double tolscr = TOLSRC;
    double tolscr2 = tolscr * tolscr;
    double *totalcalls = (double *) malloc (sizeof (double) * nthreads * 64);
    assert (totalcalls != NULL);
    double *totalnintls = (double *) malloc (sizeof (double) * nthreads * 64);
    assert (totalnintls != NULL);

    #pragma omp parallel for
    for (int i = 0; i < nthreads; i++)
    {
        totalcalls[i * 64] = 0.0;
        totalnintls[i * 64] = 0.0;
    }
    
    /************************************************************/
    // init F
    #pragma omp parallel for
    for (int j = 0; j < sizeD1_aligned * nthreads; j++)
    {
        F1[j] = 0.0;
    }
    #pragma omp parallel for
    for (int j = 0; j < sizeD2_aligned * nthreads; j++)
    {
        F2[j] = 0.0;
    }
    #pragma omp parallel for 
    for (int j = 0; j < sizeD3_aligned * nthreads; j++)
    {
        F3[j] = 0.0;
    }
    if (toOffload)
    {
        offload_reset_F (mic_numdevs);
    }
    
    printf ("Compute tasks\n");
    gettimeofday (&tv1, NULL);
    // main computation
    if (toOffload == 1)
    {
        offload_fock_task (mic_numdevs, basis, erd, shellptr,
                           shellvalue, shellid, shellrid,
                           f_startind, rowpos, colpos, rowptr, colptr,
                           tolscr2, startM, startP,
                           task_startM, task_endM,
                           task_startP, task_endP,
                           D1, D2, D3, F1, F2, F3,
                           rowsize, colsize, colsize,
                           sizeD1_aligned, sizeD2_aligned, sizeD3_aligned,
                           mic_fraction, totalcalls, totalnintls);
    }
    else
    {
        fock_task (basis, erd, shellptr,
                   shellvalue, shellid, shellrid,
                   f_startind, rowpos, colpos, rowptr, colptr,
                   tolscr2, startM, startP,
                   task_startM, task_endM,
                   task_startP, task_endP,
                   D1, D2, D3, F1, F2, F3,
                   rowsize, colsize, colsize,
                   sizeD1_aligned, sizeD2_aligned, sizeD3_aligned,
                   totalcalls, totalnintls);
    }
    gettimeofday (&tv2, NULL);
    if (toOffload == 1)
    {
        offload_reduce_mic (mic_numdevs, F1_offload, F2_offload, F3_offload,
                            sizeD1_aligned, sizeD2_aligned, sizeD3_aligned);
    }
    // reduction on CPU
    for (int i = 1; i < nthreads; i++)
    {
        #pragma omp parallel
        {
            #pragma omp for
            for (int j = 0; j < sizeD1_aligned; j++)
            {
                F1[j + 0 * sizeD1_aligned] += F1[j + i * sizeD1_aligned];
            }
            #pragma omp for
            for (int j = 0; j < sizeD2_aligned; j++)
            {
                F2[j + 0 * sizeD2_aligned] += F2[j + i * sizeD2_aligned];
            }
            #pragma omp for 
            for (int j = 0; j < sizeD3_aligned; j++)
            {
                F3[j + 0 * sizeD3_aligned] += F3[j + i * sizeD3_aligned];
            }
        }
    }
    
    if (toOffload == 1)
    {
        offload_wait_mic (mic_numdevs);
        offload_reduce (mic_numdevs,
                        F1, F2, F3,
                        F1_offload, F2_offload, F3_offload,
                        sizeD1_aligned, sizeD2_aligned, sizeD3_aligned);
    }
    /***********************************************************/

    timepass = (tv2.tv_sec - tv1.tv_sec) + (tv2.tv_usec - tv1.tv_usec) / 1e6;
    printf ("takes %.3lf secs\n", timepass);

#if 1
    FILE *fp;
    fp = fopen ("F.dat", "w+");
    assert (fp != NULL);

    for (int j = 0; j < sizeD1; j++)
        fprintf (fp, "%le\n", F1[j]);
    for (int j = 0; j < sizeD2; j++)
        fprintf (fp, "%le\n", F2[j]);
    for (int j = 0; j < sizeD3; j++)
        fprintf (fp, "%le\n", F3[j]);

    fclose (fp);
#endif

    if (toOffload == 1)
    {
        offload_deinit (mic_numdevs, shellptr, shellid, shellrid,
                        shellvalue, f_startind, rowpos, colpos,
                        rowptr, colptr, D1, D2, D3,
                        VD1, VD2, VD3, F1_offload, F2_offload, F3_offload);
        CInt_offload_destroyERD (erd);
        CInt_offload_destroyBasisSet (basis);
    }
    else
    {
        CInt_destroyERD (erd);
        CInt_destroyBasisSet (basis);
    }
    _mm_free (F1);
    _mm_free (F2);
    _mm_free (F3);
    _mm_free (D1);
    _mm_free (D2);
    _mm_free (D3);
    _mm_free (VD1);
    _mm_free (VD2);
    _mm_free (VD3);
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

    return 0;
}
