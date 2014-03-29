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
    int toOffload = 0;
    int num_devices = 0;
    struct timeval tv1, tv2;
    double timepass;

    if (argc != 12)
    {
        printf ("Usage: %s <basis set> <xyz> <nthreads>"
                " <startM> <endN> <startP> <endP> <ntasks>"
                " <toOffload> <nthreads_mic>\n",
                argv[0]);
        return -1;
    }

    const int nthreads = atoi (argv[3]);
    omp_set_num_threads (nthreads);
    const int nthreads_mic = atoi (argv[10]);
    startM = atoi (argv[4]);
    endM = atoi (argv[5]);
    startP = atoi (argv[6]);
    endP = atoi (argv[7]);
    sizetask = atoi (argv[8]);
    toOffload = atoi (argv[9]);
    double mic_fraction = atof (argv[11]);

    // load basis set and create ERD_t
    BasisSet_t basis;
    ERD_t erd;
    if (toOffload == 1)
    {
        CInt_offload_createBasisSet (&basis);
        CInt_offload_loadBasisSet (basis, argv[1], argv[2]);
        CInt_offload_createERD (basis, &erd, nthreads, nthreads_mic);
    }
    else
    {
        CInt_createBasisSet (&basis);
        CInt_loadBasisSet (basis, argv[1], argv[2]);
        CInt_createERD (basis, &erd, nthreads);
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
    assert (sizetask > 0);
    assert (startM >= 0 && startM <= endM);
    assert (startP >= 0 && startP <= endP);
    assert (startM + sizetask <= endM + 1);
    assert (startP + sizetask <= endP + 1);
    printf ("\ncreate FD for (%d:%d,0:%d|%d:%d,0:%d)\n",
            startM, endM, nshells - 1, startP, endP, nshells - 1);
    printf ("compute (%d:%d,0:%d|%d:%d,0:%d)\n",
            startM, startM + sizetask - 1, nshells - 1,
            startP, startP + sizetask - 1, nshells - 1);

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
    double *F1 = (double *) _mm_malloc (sizeof (double) * sizeD1_aligned * nthreads, 64);
    double *F2 = (double *) _mm_malloc (sizeof (double) * sizeD2_aligned * nthreads, 64);
    double *F3 = (double *) _mm_malloc (sizeof (double) * sizeD3_aligned * nthreads, 64);   
    totalFDsize += 1.0 * sizeof (double) *
        (sizeD1_aligned + sizeD2_aligned + sizeD3_aligned);
    totalFDsize += 1.0 * sizeof (double) *
        (sizeD1_aligned + sizeD2_aligned + sizeD3_aligned) * nthreads;
    assert (D1 != NULL &&
            D2 != NULL &&
            D3 != NULL &&
            F1 != NULL &&
            F2 != NULL &&
            F3 != NULL);
    printf ("use %.3lf MB\n", totalFDsize / 1024.0 / 1024.0);

    double *F1_offload = NULL;
    double *F2_offload = NULL;
    double *F3_offload = NULL;
    if (toOffload == 1)
    {
        num_devices = _Offload_number_of_devices();
        if(num_devices == 0)
        {
            printf("No target devices available. Exiting\n");
            exit(0);
        }
        else
        {
            for(int mic_id = 0; mic_id < num_devices; mic_id++)
            {
                #pragma offload target(mic:mic_id) in(nthreads_mic)
                {
                    omp_set_num_threads (nthreads_mic);
                    #pragma omp parallel num_threads(nthreads_mic)
                    {
                        int tid = omp_get_thread_num();
                    }  
                }
            }
        }
        printf("Number of Target devices installed: %d\n\n",num_devices);
        
        if (mic_fraction * num_devices >= 1)
        {
            printf ("Invalid mic_fraction value. It should be < %.3lf\n",
                    ((double) 1) / num_devices);
            exit (0);
        }
        F1_offload =
            (double *) _mm_malloc (sizeof (double) * sizeD1_aligned * num_devices, 64);
        F2_offload =
            (double *) _mm_malloc (sizeof (double) * sizeD2_aligned * num_devices, 64);
        F3_offload =
            (double *) _mm_malloc (sizeof (double) * sizeD3_aligned * num_devices, 64);
        assert (F1_offload != NULL);
        assert (F2_offload != NULL);
        assert (F3_offload != NULL);
        offload_init (num_devices, nshells, nnz, shellptr, shellid, shellrid,
                      shellvalue, f_startind, rowpos, colpos, rowptr, colptr,
                      D1, D2, D3,
                      F1_offload, F2_offload, F3_offload,
                      sizeD1_aligned, sizeD2_aligned,
                      sizeD3_aligned, nthreads_mic);      
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
        offload_copy_D (num_devices, D1, sizeD1_aligned);
        offload_copy_D (num_devices, D2, sizeD2_aligned);
        offload_copy_D (num_devices, D3, sizeD3_aligned);
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
        offload_reset_F (num_devices);
    }

    
    printf ("Compute tasks\n");
    gettimeofday (&tv1, NULL);
    // main computation
    compute_task (num_devices, basis, erd, shellptr,
                  shellvalue, shellid, shellrid,
                  f_startind, rowpos, colpos, rowptr, colptr,
                  tolscr2, startM, startP,
                  startM, startM + sizetask - 1,
                  startP, startP + sizetask - 1,
                  D1, D2, D3, F1, F2, F3,
                  rowsize, colsize, colsize,
                  sizeD1_aligned, sizeD2_aligned, sizeD3_aligned,
                  mic_fraction, totalcalls, totalnintls, toOffload);
    gettimeofday (&tv2, NULL);
    if (toOffload == 1)
    {
        offload_reduce_mic (num_devices, F1_offload, F2_offload, F3_offload,
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
        offload_wait_mic (num_devices);
        offload_reduce (num_devices,
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
    if (toOffload == 1)
    {
        _mm_free (F1_offload);
        _mm_free (F2_offload);
        _mm_free (F3_offload);
    }
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
