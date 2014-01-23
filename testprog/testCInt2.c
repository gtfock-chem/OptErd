#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>
#include <math.h>
#include <sys/time.h>

#include <yepPredefines.h>
#include "CInt.h"
#include "screening.h"

extern uint64_t erd__set_ij_kl_pairs_ticks[256];
extern uint64_t erd__rys_roots_weights_ticks[256];
extern uint64_t erd__2d_coefficients_ticks[256];
extern uint64_t erd__2d_pq_integrals_ticks[256];
extern uint64_t erd__int2d_to_e0f0_ticks[256];
extern uint64_t erd__e0f0_pcgto_block_ticks[256];
extern uint64_t erd__xyz_to_ry_abcd_ticks[256];
extern uint64_t erd__hrr_matrix_ticks[256];
extern uint64_t erd__hrr_transform_ticks[256];
extern uint64_t erd__csgto_ticks[256];

void initProfile(int nthreads)
{
    memset(erd__set_ij_kl_pairs_ticks, 0, nthreads * sizeof(int));
    memset(erd__rys_roots_weights_ticks, 0, nthreads * sizeof(int));
    memset(erd__2d_coefficients_ticks, 0, nthreads * sizeof(int));
    memset(erd__2d_pq_integrals_ticks, 0, nthreads * sizeof(int));
    memset(erd__int2d_to_e0f0_ticks, 0, nthreads * sizeof(int));
    memset(erd__e0f0_pcgto_block_ticks, 0, nthreads * sizeof(int));
    memset(erd__xyz_to_ry_abcd_ticks, 0, nthreads * sizeof(int));
    memset(erd__hrr_matrix_ticks, 0, nthreads * sizeof(int));
    memset(erd__hrr_transform_ticks, 0, nthreads * sizeof(int));
    memset(erd__csgto_ticks, 0, nthreads * sizeof(int));
}

void printProfile(int nthreads, uint64_t freq)
{
    double per_thread_total_secs[YEP_COUNT_OF(erd__set_ij_kl_pairs_ticks)] = { 0.0 };

    printf("%22s", "erd__set_ij_kl_pairs");
    for (int i = 0; i < nthreads; i++) {
        printf("\t%.3lf", ((double) erd__set_ij_kl_pairs_ticks[i]) / freq);
        per_thread_total_secs[i] += (double) erd__set_ij_kl_pairs_ticks[i] / freq;
    }
    printf("\n");

    printf("%22s", "erd__rys_roots_weights");
    for (int i = 0; i < nthreads; i++) {
        printf("\t%.3lf", ((double) erd__rys_roots_weights_ticks[i]) / freq);
        per_thread_total_secs[i] += (double) erd__rys_roots_weights_ticks[i] / freq;
    }
    printf("\n");

    printf("%22s", "erd__2d_coefficients");
    for (int i = 0; i < nthreads; i++) {
        printf("\t%.3lf", ((double) erd__2d_coefficients_ticks[i]) / freq);
        per_thread_total_secs[i] += (double) erd__2d_coefficients_ticks[i] / freq;
    }
    printf("\n");

    printf("%22s", "erd__2d_pq_integrals");
    for (int i = 0; i < nthreads; i++) {
        printf("\t%.3lf", ((double) erd__2d_pq_integrals_ticks[i]) / freq);
        per_thread_total_secs[i] += (double) erd__2d_pq_integrals_ticks[i] / freq;
    }
    printf("\n");

    printf("%22s", "erd__int2d_to_e0f0");
    for (int i = 0; i < nthreads; i++) {
        printf("\t%.3lf", ((double) erd__int2d_to_e0f0_ticks[i]) / freq);
        per_thread_total_secs[i] += (double) erd__int2d_to_e0f0_ticks[i] / freq;
    }
    printf("\n");

    printf("%22s", "erd__e0f0_pcgto");
    for (int i = 0; i < nthreads; i++) {
        printf("\t%.3lf", ((double) erd__e0f0_pcgto_block_ticks[i]) / freq);
        per_thread_total_secs[i] += (double) erd__e0f0_pcgto_block_ticks[i] / freq;
    }
    printf("\n");

    printf("%22s", "erd__xyz_to_ry_abcd");
    for (int i = 0; i < nthreads; i++) {
        printf("\t%.3lf", ((double) erd__xyz_to_ry_abcd_ticks[i]) / freq);
        per_thread_total_secs[i] += (double) erd__xyz_to_ry_abcd_ticks[i] / freq;
    }
    printf("\n");

    printf("%22s", "erd__hrr_matrix");
    for (int i = 0; i < nthreads; i++) {
        printf("\t%.3lf", ((double) erd__hrr_matrix_ticks[i]) / freq);
        per_thread_total_secs[i] += (double) erd__hrr_matrix_ticks[i] / freq;
    }
    printf("\n");

    printf("%22s", "erd__hrr_transform");
    for (int i = 0; i < nthreads; i++) {
        printf("\t%.3lf", ((double) erd__hrr_transform_ticks[i]) / freq);
        per_thread_total_secs[i] += (double) erd__hrr_transform_ticks[i] / freq;
    }
    printf("\n");
#if 0
    printf("%22s", "erd__csgto");
    for (int i = 0; i < nthreads; i++) {
        printf("\t%.3lf", ((double) erd__csgto_ticks[i]) / freq);
        per_thread_total_secs[i] += (double) erd__csgto_ticks[i] / freq;
    }
    printf("\n");
    printf("%22s", "total");
    for (int i = 0; i < nthreads; i++) {
        printf("\t%.3lf", per_thread_total_secs[i]);
    }
#endif
    printf("\n");

}

int
main (int argc, char **argv)
{
    BasisSet_t basis;
    ERD_t *erd;
    int ns;
    int i;
    int nthreads;
    double *totalcalls = 0;
    double *totalnintls = 0;
    int *shellptr;
    int *shellid;
    int *shellrid;
    double *shellvalue;
    double fraction;
    
    struct timeval tv1;
    struct timeval tv2;
    double timepass;
    uint64_t start_clock, end_clock;
    
    if (argc != 5)
    {
        printf ("Usage: %s <basisset> <xyz> <fraction> <nthreads>\n", argv[0]);
        return -1;
    }

    start_clock = __rdtsc();
    sleep(1);
    end_clock = __rdtsc();
    uint64_t freq = end_clock - start_clock;
    printf("freq = %lld\n", freq);

    fraction = atof (argv[3]);
    assert (fraction > 0.0 && fraction <= 1.0);
    nthreads = atoi (argv[4]);
    initProfile(nthreads);
    omp_set_num_threads (nthreads);

    // load basis set
    CInt_createBasisSet (&basis);
    CInt_loadBasisSet (basis, argv[1], argv[2]);
    schwartz_screening (basis, &shellptr, &shellid, &shellrid, &shellvalue);
    
    printf ("Molecule info:\n");
    printf ("  #Atoms\t= %d\n", CInt_getNumAtoms (basis));
    printf ("  #Shells\t= %d\n", CInt_getNumShells (basis));
    printf ("  #Funcs\t= %d\n", CInt_getNumFuncs (basis));
    printf ("  #OccOrb\t= %d\n", CInt_getNumOccOrb (basis));
    printf ("  nthreads\t= %d\n", nthreads);

    erd = (ERD_t *) malloc (sizeof (ERD_t) * nthreads);
    assert (erd != NULL);
    totalcalls = (double *) malloc (sizeof (double) * nthreads * 64);
    assert (totalcalls != NULL);
    totalnintls = (double *) malloc (sizeof (double) * nthreads * 64);
    assert (totalnintls != NULL);
    
#pragma omp parallel for
    for (i = 0; i < nthreads; i++)
    {
        CInt_createERD (basis, &(erd[i]));
        totalcalls[i * 64] = 0.0;
        totalnintls[i * 64] = 0.0;
    }

    printf ("Computing integrals ...\n");

    //printf ("max memory footprint per thread = %lf KB\n",
    //    CInt_getMaxMemory (erd[0])/1024.0);
    
    ns = CInt_getNumShells (basis);
    timepass = 0.0;
    gettimeofday (&tv1, NULL);

    int slen;
    int *idx;
    int swap;
    int tmpid;
    int reallen;
    
    slen = shellptr[ns];
    idx = (int *) malloc (sizeof (int) * slen);
    assert (idx != NULL);
    srand (1234);
    for (i = 0; i < slen; i++)
    {
        idx[i] = i;
    }
    for (i = 0; i < slen - 1; i++)
    {
        swap = (int)((double)rand()/RAND_MAX * (slen - i));
        tmpid = idx[swap];
        idx[swap] = idx[slen - i - 1];
        idx[slen - i - 1] = tmpid;
    }
    reallen = (int)(slen * fraction);
    reallen = ((reallen == 0) ? 1 : reallen);    
    start_clock = __rdtsc(); 
#pragma omp parallel
    {
        int tid;
        int M;
        int N;
        int P;
        int Q;
        int k;
        int i;
        int j;
        double *integrals;
        int nints;
        double value1;
        double value2;
        int start2;
        int end2;

        tid = omp_get_thread_num ();

#pragma omp for schedule(dynamic)
        for (k = 0; k < reallen; k++)
        {
            i = idx[k];
            M = shellrid[i];
            N = shellid[i];
            value1 = shellvalue[i];
            for (P = 0; P < ns; P++)
            {                               
                start2 = shellptr[P];
                end2 = shellptr[P + 1];
                for (j = start2; j < end2; j++)
                {                
                    Q = shellid[j];
                    value2 = shellvalue[j];
                    if (M > N || P > Q || (M + N) > (P + Q))
                        continue;
                    if (fabs(value1 * value2) < TOLSRC * TOLSRC)
                        continue;
                    totalcalls[tid * 64] = totalcalls[tid * 64] + 1;

                    CInt_computeShellQuartet (basis, erd[tid], M, N, P, Q, &integrals,
                                              &nints);

                    totalnintls[tid * 64] = totalnintls[tid * 64] + nints;
                }
            }
        }
    }

    for (i = 1; i < nthreads; i++)
    {
        totalcalls[0 * 64] = totalcalls[0 * 64] + totalcalls[i * 64];
        totalnintls[0 * 64] = totalnintls[0 * 64] + totalnintls[i * 64];
    }
    end_clock = __rdtsc();
    gettimeofday (&tv2, NULL);
    uint64_t total_ticks = end_clock - start_clock;
    timepass = ((double) total_ticks) / freq;
    printf ("Done\n");
    printf ("\n");
    printf ("Number of calls: %.6le, Number of integrals: %.6le\n",
            totalcalls[0], totalnintls[0]);
    printf ("Total GigaTicks: %.3lf\n",
            (double) (total_ticks) * 1.0e-9);
    printf ("Total time: %.4lf secs\n", timepass);
    printf ("Average time per call: %.3le us\n",
            1000.0 * 1000.0 * timepass / totalcalls[0]);

    printProfile(nthreads, freq);

    for (i = 0; i < nthreads; i++)
    {
        CInt_destroyERD (erd[i]);
    }
    free (erd);
    free (totalcalls);
    free (totalnintls);
    free (shellptr);
    free (shellid);
    free (shellvalue);
    free (shellrid);
    
    CInt_destroyBasisSet (basis);

    return 0;
}
