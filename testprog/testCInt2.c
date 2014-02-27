#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <sys/time.h>

#include <screening.h>
#include <erd_profile.h>

static uint64_t get_cpu_frequency(void) {
    const uint64_t start_clock = __rdtsc();
    sleep(1);
    const uint64_t end_clock = __rdtsc();
    return end_clock - start_clock;
}

int main (int argc, char **argv)
{
    int nnz;
    int *shellptr;
    int *shellid;
    int *shellrid;
    double *shellvalue;
    if (argc != 5) {
        printf ("Usage: %s <basisset> <xyz> <fraction> <nthreads>\n", argv[0]);
        return -1;
    }

    const uint64_t freq = get_cpu_frequency();

    const double fraction = atof (argv[3]);
    assert(fraction > 0.0 && fraction <= 1.0);
    const int nthreads = atoi(argv[4]);
    #ifdef _OPENMP
    omp_set_num_threads(nthreads);
    #else
    assert(nthreads == 1);
    #endif
    
    // load basis set
    BasisSet_t basis;
    CInt_createBasisSet(&basis);
    CInt_loadBasisSet(basis, argv[1], argv[2]);
    schwartz_screening(basis, &shellptr, &shellid, &shellrid, &shellvalue, &nnz);

    printf("Molecule info:\n");
    printf("  #Atoms\t= %d\n", CInt_getNumAtoms(basis));
    printf("  #Shells\t= %d\n", CInt_getNumShells(basis));
    printf("  #Funcs\t= %d\n", CInt_getNumFuncs(basis));
    printf("  #OccOrb\t= %d\n", CInt_getNumOccOrb(basis));
    printf("  nthreads\t= %d\n", nthreads);

    ERD_t erd;
    CInt_createERD(basis, &erd, nthreads);

    double* totalcalls = (double *) malloc (sizeof (double) * nthreads * 64);
    assert(totalcalls != NULL);
    double* totalnintls = (double *) malloc (sizeof (double) * nthreads * 64);
    assert(totalnintls != NULL);
    
    #pragma omp parallel for
    for (int i = 0; i < nthreads; i++) {
        totalcalls[i * 64] = 0.0;
        totalnintls[i * 64] = 0.0;
    }
    printf("Computing integrals ...\n");

    // reset profiler
    erd_reset_profile ();
        
    //printf ("max memory footprint per thread = %lf KB\n",
    //    CInt_getMaxMemory (erd[0])/1024.0);

    const int ns = CInt_getNumShells(basis);

    srand(1234);
    /* In (fraction) cases rand() returns value not greater than computationThreshold */
    const int computationThreshold = lround(fraction * RAND_MAX);

    const uint64_t start_clock = __rdtsc(); 
    #pragma omp parallel
    {
        #ifdef _OPENMP
        const int tid = omp_get_thread_num();
        #else
        const int tid = 0;
        #endif

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < shellptr[ns]; i++) {
            const int M = shellrid[i];
            const int N = shellid[i];
            if (M > N)
                continue;

            const double value1 = shellvalue[i];
            for (int P = 0; P < ns; P++) {
                const int start2 = shellptr[P];
                const int end2 = shellptr[P + 1];
                for (int j = start2; j < end2; j++) {
                    const int Q = shellid[j];
                    if (P > Q)
                        continue;
                    if ((M + N) > (P + Q))
                        continue;

                    const double value2 = shellvalue[j];
                    if (fabs(value1 * value2) < TOLSRC * TOLSRC)
                        continue;

                    /* Sample random integer. With probability (fraction) process the shell quartet. */
                    if (rand() <= computationThreshold) {
                        double *integrals;
                        int nints;
                        CInt_computeShellQuartet(basis, erd, tid, M, N, P, Q, &integrals, &nints);

                        totalcalls[tid * 64] += 1;
                        totalnintls[tid * 64] += nints;
                    }
                }
            }
        }
    }
    const uint64_t end_clock = __rdtsc();

    for (int i = 1; i < nthreads; i++) {
        totalcalls[0 * 64] = totalcalls[0 * 64] + totalcalls[i * 64];
        totalnintls[0 * 64] = totalnintls[0 * 64] + totalnintls[i * 64];
    } 
    const uint64_t total_ticks = end_clock - start_clock;
    const double timepass = ((double) total_ticks) / freq;

    printf("Done\n");
    printf("\n");
    printf("Number of calls: %.6le, Number of integrals: %.6le\n", totalcalls[0], totalnintls[0]);
    printf("Total GigaTicks: %.3lf, freq = %.3lf GHz\n", (double) (total_ticks) * 1.0e-9, (double)freq/1.0e9);
    printf("Total time: %.4lf secs\n", timepass);
    printf("Average time per call: %.3le us\n", 1000.0 * 1000.0 * timepass / totalcalls[0]);

    // use 1 if thread timing is not required
    erd_print_profile(1);

    CInt_destroyERD(erd);
    free(totalcalls);
    free(totalnintls);
    free(shellptr);
    free(shellid);
    free(shellvalue);
    free(shellrid);

    CInt_destroyBasisSet(basis);

    return 0;
}
