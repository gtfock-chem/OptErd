#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h>

#include "erd_profile.h"


__declspec(align(256)) uint64_t erd_ticks[MAXTHREADS][NUM_TICKS + 8];

static char ticks_name[NUM_TICKS][128] = 
{
    "erd__prepare_ctr",
    "erd__set_ij_kl_pairs",
    "erd__rys_roots_weights",
    "erd__2d_coefficients",
    "erd__2d_pq_integrals",
    "erd__int2d_to_e0f0",
    "@erd__e0f0_pcgto_block",
    "erd__xyz_to_ry_abcd",
    "erd__hrr_matrix",
    "erd__hrr_transform",
    "erd__move_ry",
    "erd__spherical_transform",
    "@erd__csgto",
    
    "erd__prepare_ctr_1111",
    "erd__set_ij_kl_pairs_1111",
    "erd__ssss_pcgto_block",
    "erd__sssp_pcgto_block",
    "erd__sspp_pcgto_block",
    "erd__sppp_pcgto_block",
    "erd__pppp_pcgto_block",
    "@erd__1111_csgto"
};


void erd_reset_profile (void)
{
    int i;
    int nthreads;

    nthreads = omp_get_max_threads ();
    assert (nthreads <= MAXTHREADS);
    printf ("Reset profiler for %d threads\n", nthreads);
    #pragma omp parallel for
    for (i = 0; i < nthreads; i++)
    {
        memset (erd_ticks[i], 0, sizeof(uint64_t) * NUM_TICKS);
    }
}


void erd_print_profile (int mode)
{
    int k;
    int i;
    double total_secs = 0.0;
    int nthreads;
    uint64_t start_clock;
    uint64_t end_clock;
    uint64_t freq;

    nthreads = omp_get_max_threads ();
    start_clock = __rdtsc();
    sleep(1);
    end_clock = __rdtsc();
    freq = end_clock - start_clock;

    printf ("\n");    
    printf("freq = %.3lf GHz\n", (double)freq/1e9); 
    for (k = 0; k < NUM_TICKS; k++)
    {
        total_secs = 0.0;
        printf("%25s", ticks_name[k]);
        for (i = 0; i < nthreads; i++)
        {
            if (mode == 0)
            {
                printf("\t%.3lf", (double)erd_ticks[i][k]/freq);
            }
            total_secs += (double)erd_ticks[i][k]/freq;
        }        
        printf(":\t%.3lf", total_secs/nthreads);
        printf("\n");
    }
}


int start_timer_ (uint64_t *stime)
{
    *stime  = __rdtsc();

    return 0;
}


int end_timer_ (int *idx, uint64_t *stime)
{
    int tid = omp_get_thread_num();   
    uint64_t etime = __rdtsc();
    erd_ticks[tid][*idx] += etime - *stime;

    return 0;
}
