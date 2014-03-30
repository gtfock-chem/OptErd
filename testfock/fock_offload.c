#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <sys/time.h>
#include <offload.h>
#include <immintrin.h>
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#include "CInt.h"
#include "config.h"
#include "fock_offload.h"


#define CHUNK_SIZE 8
#define MIC_MIN_WORK_SIZE 20
#define PFD0 1
#define PFD1 12
#define ALLOC alloc_if(1) free_if(0)
#define REUSE alloc_if(0) free_if(0)
#define ONCE  alloc_if(1) free_if(1)
#define FREE  alloc_if(0) free_if(1)
#define ALIGNED_8(size) ((((size) + 7)/8)*8)
#define NUM_F_COPIES_MIC(n) (((n)+3) >> 2)
#define MY_F_COPY_MIC(tid) ((tid) >> 2)


#pragma offload_attribute(push, target(mic))

pfock_mic_t pfock_mic;

static inline void update_F (double *integrals, int dimM, int dimN, int dimP, int dimQ,
                      int flag1, int flag2, int flag3,
                      int iMN, int iPQ, int iMP, int iNP, int iMQ, int iNQ,
                      double *D1, double *D2, double *D3,
                      double *J1, double *J2, double *K3,
                      int ldX1, int ldX2, int ldX3)
{
    int iM;
    int iN;
    int iP;
    int iQ;
    double I;
    int flag4;
    int flag5;
    int flag6;
    int flag7;
    double *D_MN;
    double *D_PQ;
    double *D_MP;
    double *D_NP;
    double *D_MQ;
    double *D_NQ;
    double *J_MN;
    double *J_PQ;
    double *K_MP;
    double *K_NP;
    double *K_MQ;
    double *K_NQ;

    flag4 = (flag1 == 1 && flag2 == 1) ? 1 : 0;
    flag5 = (flag1 == 1 && flag3 == 1) ? 1 : 0;
    flag6 = (flag2 == 1 && flag3 == 1) ? 1 : 0;
    flag7 = (flag4 == 1 && flag3 == 1) ? 1 : 0;

    D_MN = D1 + iMN;
    D_PQ = D2 + iPQ;
    D_MP = D3 + iMP;
    D_MQ = D3 + iMQ;
    D_NP = D3 + iNP;
    D_NQ = D3 + iNQ;
    J_MN = J1 + iMN;
    J_PQ = J2 + iPQ;
    K_MP = K3 + iMP;
    K_MQ = K3 + iMQ;
    K_NP = K3 + iNP;
    K_NQ = K3 + iNQ;


#ifdef __MIC__
    _mm_prefetch((const char* ) D_MN, _MM_HINT_T0);
    _mm_prefetch((const char* ) D_PQ, _MM_HINT_T0);
    _mm_prefetch((const char* ) D_MP, _MM_HINT_T0);
    _mm_prefetch((const char* ) D_MQ, _MM_HINT_T0);
    _mm_prefetch((const char* ) D_NP, _MM_HINT_T0);
    _mm_prefetch((const char* ) D_NQ, _MM_HINT_T0);
    _mm_prefetch((const char* ) J_MN, _MM_HINT_T0);
    _mm_prefetch((const char* ) J_PQ, _MM_HINT_T0);
    _mm_prefetch((const char* ) K_MP, _MM_HINT_T0);
    _mm_prefetch((const char* ) K_MQ, _MM_HINT_T0);
    _mm_prefetch((const char* ) K_NP, _MM_HINT_T0);
    _mm_prefetch((const char* ) K_NQ, _MM_HINT_T0);
#endif

    for (iM = 0; iM < dimM; iM++)
    {
        for (iN = 0; iN < dimN; iN++)
        {
            int iMN = iM * ldX1 + iN;
            double j_MN = 0;
            for (iP = 0; iP < dimP; iP++)
            {
                int iMP = iM * ldX3 + iP;
                int iNP = iN * ldX3 + iP;
                double k_MP = 0;
                double k_NP = 0;

                for (iQ = 0; iQ < dimQ; iQ++)
                {
                    int iPQ = iP * ldX2 + iQ;
                    int iMQ = iM * ldX3 + iQ;
                    int iNQ = iN * ldX3 + iQ;

                    I = integrals[iM + dimM * (iN + dimN * (iP + dimP * iQ))];

                    // F(m, n) += D(p, q) * 2 * I(m, n, p, q)
                    // F(n, m) += D(p, q) * 2 * I(n, m, p, q)
                    // F(m, n) += D(q, p) * 2 * I(m, n, q, p)
                    // F(n, m) += D(q, p) * 2 * I(n, m, q, p)
                    double vMN = 1.0 * (1 + flag1 + flag2 + flag4) * D_PQ[iPQ] * I;
                    j_MN += vMN;

                    // F(p, q) += D(m, n) * 2 * I(p, q, m, n)
                    // F(p, q) += D(n, m) * 2 * I(p, q, n, m)
                    // F(q, p) += D(m, n) * 2 * I(q, p, m, n)
                    // F(q, p) += D(n, m) * 2 * I(q, p, n, m)
                    double vPQ = 1.0 * (flag3 + flag5 + flag6 +
                               flag7) * D_MN[iMN] * I;
#ifdef __MIC__
                    #pragma omp atomic update
#endif
                    J_PQ[iPQ] += vPQ;

                    // F(m, p) -= D(n, q) * I(m, n, p, q)
                    // F(p, m) -= D(q, n) * I(p, q, m, n)
                    double vMP = (1 + flag3) * 0.5 * D_NQ[iNQ] * I;
                    k_MP -= vMP;

                    // F(n, p) -= D(m, q) * I(n, m, p, q)
                    // F(p, n) -= D(q, m) * I(p, q, n, m)
                    double vNP = (flag1 + flag5) * 0.5 * D_MQ[iMQ] * I;
                    k_NP -= vNP;

                    // F(m, q) -= D(n, p) * I(m, n, q, p)
                    // F(q, m) -= D(p, n) * I(q, p, m, n)
                    double vMQ = (flag2 + flag6) * 0.5 * D_NP[iNP] * I;
#ifdef __MIC__
                    #pragma omp atomic update
#endif
                    K_MQ[iMQ] -= vMQ;

                    // F(n, q) -= D(m, p) * I(n, m, q, p)
                    // F(q, n) -= D(p, m) * I(q, p, n, m)
                    double vNQ = (flag4 + flag7) * 0.5 * D_MP[iMP] * I;
#ifdef __MIC__
                    #pragma omp atomic update
#endif
                    K_NQ[iNQ] -= vNQ;
                }
#ifdef __MIC__
                #pragma omp atomic update
#endif
                K_MP[iMP] += k_MP;
#ifdef __MIC__
                #pragma omp atomic update
#endif
                K_NP[iNP] += k_NP;
            }
#ifdef __MIC__
            #pragma omp atomic update
#endif
            J_MN[iMN] += j_MN;
        }
    }
}


static void compute_chunk (BasisSet_t basis, ERD_t erd,
                    int *shellptr, double *shellvalue,
                    int *shellid, int *shellrid, int *f_startind,
                    int *rowpos, int *colpos, int *rowptr, int *colptr,
                    double tolscr2, int startrow, int startcol,
                    int startChunkMN, int endChunkMN, int startChunkPQ, int endChunkPQ,
                    double *D1, double *D2, double *D3,
                    double *J1, double *J2, double *K3,
                    int ldX1, int ldX2, int ldX3,
                    int sizeX1, int sizeX2, int sizeX3,
                    double *totalcalls, double *totalnintls,
                    int nt)
{
    int i, j;
    double value1;
    double value2;
    int dimM;
    int dimN;
    int dimP;
    int dimQ;
    int flag1;
    int flag2;
    int flag3;
    int iX1M;
    int iX2P;
    int iX3M;
    int iX3N;
    int iX3P;
    int iX3Q;
    int iX1N;
    int iX2Q;
    int iMN;
    int iPQ;
    int iMP;
    int iNP;
    int iMQ;
    int iNQ;
    int M;
    int N;
    int P;
    int Q;
    double *integrals;
    int nints;

    for (i = startChunkMN; i < endChunkMN; i++)
    {
        M = shellrid[i];
        N = shellid[i];
        value1 = shellvalue[i];
        dimM = f_startind[M + 1] - f_startind[M];
        dimN = f_startind[N + 1] - f_startind[N];
        iX3M = rowpos[M];
        iX1M = f_startind[M] - f_startind[startrow];
        iX1N = iX3N = rowptr[i];
        iMN = iX1M * ldX1 + iX1N;
        flag1 = (value1 < 0.0) ? 1 : 0;
        for (j = startChunkPQ; j < endChunkPQ; j++)
        {
            P = shellrid[j];
            Q = shellid[j];
            if ((M > P && (M + P) % 2 == 1) ||
                    (M < P && (M + P) % 2 == 0))
                continue;
            if ((M == P) &&
                    ((N > Q && (N + Q) % 2 == 1) ||
                     (N < Q && (N + Q) % 2 == 0)))
                continue;
            value2 = shellvalue[j];
            dimP = f_startind[P + 1] - f_startind[P];
            dimQ = f_startind[Q + 1] - f_startind[Q];
            iX3P = colpos[P];
            iX2P = f_startind[P] - f_startind[startcol];
            iMP = iX3M * ldX3 + iX3P;
            iNP = iX3N * ldX3 + iX3P;
            iX3Q = iX2Q = colptr[j];
            iPQ = iX2P * ldX2 + iX2Q;
            iMQ = iX3M * ldX3 + iX3Q;
            iNQ = iX3N * ldX3 + iX3Q;
            flag3 = (M == P && Q == N) ? 0 : 1;
            flag2 = (value2 < 0.0) ? 1 : 0;
            if (fabs (value1 * value2) >= tolscr2)
            {
                CInt_computeShellQuartet (basis, erd, nt,
                        M, N, P, Q, &integrals, &nints);
                if (nints != 0)
                {
                    update_F (integrals, dimM, dimN, dimP, dimQ,
                              flag1, flag2, flag3,
                              iMN, iPQ, iMP, iNP, iMQ, iNQ,
                              D1, D2, D3, J1, J2, K3, ldX1, ldX2, ldX3);
                }
            }
        }
    }
}



static void compute_block_of_chunks (BasisSet_t basis, ERD_t erd,
                   int *shellptr, double *shellvalue,
                   int *shellid, int *shellrid, int *f_startind,
                   int *rowpos, int *colpos, int *rowptr, int *colptr,
                   double tolscr2, int startrow, int startcol,
                   int startMN, int endMN, int startPQ, int endPQ,
                   int startChunk, int endChunk, int chunksMN,
                   double *D1, double *D2, double *D3,
                   int ldX1, int ldX2, int ldX3,
                   int sizeX1, int sizeX2, int sizeX3,
                   double *totalcalls, double *totalnintls)
{
    int k;

    #pragma omp parallel
    {
        int tid;
        
        tid = omp_get_thread_num ();
        double *J1 = pfock_mic.F1[MY_F_COPY_MIC(tid)];
        double *J2 = pfock_mic.F2[MY_F_COPY_MIC(tid)];
        double *K3 = pfock_mic.F3[MY_F_COPY_MIC(tid)];

        #pragma omp for schedule(dynamic)
        for(k = startChunk; k < endChunk; k++)
        {
            int chunkIdMN = k / chunksMN;
            int chunkIdPQ = k % chunksMN;

            int startChunkMN = startMN + chunkIdMN * CHUNK_SIZE;
            int endChunkMN   = startChunkMN + CHUNK_SIZE;
            if(endChunkMN > endMN) endChunkMN = endMN;
            int startChunkPQ = startPQ + chunkIdPQ * CHUNK_SIZE;
            int endChunkPQ   = startChunkPQ + CHUNK_SIZE;
            if(endChunkPQ > endPQ) endChunkPQ = endPQ;

            compute_chunk (basis, erd,
                    shellptr, shellvalue,
                    shellid, shellrid, f_startind,
                    rowpos, colpos, rowptr, colptr,
                    tolscr2, startrow, startcol,
                    startChunkMN, endChunkMN, startChunkPQ, endChunkPQ,
                    D1, D2, D3,
                    J1, J2, K3,
                    ldX1, ldX2, ldX3,
                    sizeX1, sizeX2, sizeX3,
                    totalcalls, totalnintls, tid);
        }
    }
}

#pragma offload_attribute(pop)


void offload_reset_F (int num_devices)
{
    for(int mic_id = 0; mic_id < num_devices; mic_id++)
    {
        #pragma offload target(mic:mic_id)\
                nocopy(pfock_mic)
        {
            int num_F = pfock_mic.num_F;
            double **F1 = pfock_mic.F1;
            double **F2 = pfock_mic.F2;
            double **F3 = pfock_mic.F3;
            int sizeD1 = pfock_mic.sizeD1;
            int sizeD2 = pfock_mic.sizeD2;
            int sizeD3 = pfock_mic.sizeD3;
            #pragma omp parallel
            {
                for (int i = 1; i < num_F; i++)
                {
                    #pragma omp for
                    for (int j = 0; j < sizeD1; j++)
                    {
                        F1[i][j] = 0.0;
                    }
                    #pragma omp for
                    for (int j = 0; j < sizeD2; j++)
                    {
                        F2[i][j] = 0.0;
                    }
                    #pragma omp for
                    for (int j = 0; j < sizeD3; j++)
                    {
                        F3[i][j] = 0.0;
                    }
                }
            }
        }
    }
}


void offload_fock_task (int num_devices,
                        BasisSet_t basis, ERD_t erd,
                        int *shellptr, double *shellvalue,
                        int *shellid, int *shellrid, int *f_startind,
                        int *rowpos, int *colpos, int *rowptr, int *colptr,
                        double tolscr2, int startrow, int startcol,
                        int startM, int endM, int startP, int endP,
                        double *D1, double *D2, double *D3,
                        double *F1, double *F2, double *F3,
                        int ldX1, int ldX2, int ldX3,
                        int sizeX1, int sizeX2, int sizeX3, double mic_fraction,
                        double *totalcalls, double *totalnintls)
{
    int startMN;
    int endMN;
    int startPQ;
    int endPQ;
    startMN = shellptr[startM];
    endMN = shellptr[endM + 1];
    startPQ = shellptr[startP];
    endPQ = shellptr[endP + 1];

    int chunksMN = (endMN - startMN) / CHUNK_SIZE;
    int chunksPQ = (endPQ - startPQ) / CHUNK_SIZE;
    int totalChunks = chunksMN * chunksPQ;
    int head = 0;
    mic_fraction = 0;
    int initialChunksMIC = totalChunks * mic_fraction;
    
    head = num_devices * initialChunksMIC;    
    #pragma omp parallel
    {
        int my_chunk;
        int tid;        
        tid = omp_get_thread_num ();
        double *J1 = &(F1[tid * sizeX1]);
        double *J2 = &(F2[tid * sizeX2]);
        double *K3 = &(F3[tid * sizeX3]);

        if(0)//tid == 0)
        {
            int signalled[num_devices];
            int mic_id;
            int startChunk = 0;
            int endChunk;
            int dummy_tag = 0;
            int *finish_tag = &dummy_tag;
            long start, end;
            for(mic_id = 0; mic_id < num_devices; mic_id++)
            {
                int endChunk = startChunk + initialChunksMIC;
                start = __rdtsc();
                #pragma offload target(mic:mic_id) \
                    nocopy(basis_mic, erd_mic) \
                    in(shellptr, shellvalue, shellid, shellrid: length(0) REUSE) \
                    in(f_startind, rowpos, colpos, rowptr, colptr: length(0) REUSE) \
                    in(tolscr2, startrow, startcol) \
                    in(startMN, endMN, startPQ, endPQ) \
                    in(startChunk, endChunk, chunksMN) \
                    in(D1, D2, D3: length(0) REUSE) \
                    in(ldX1, ldX2, ldX3, sizeX1, sizeX2, sizeX3) \
                    nocopy(totalcalls, totalnintls) \
                    signal(finish_tag)
                compute_block_of_chunks (basis_mic, erd_mic,
                        shellptr, shellvalue,
                        shellid, shellrid, f_startind,
                        rowpos, colpos, rowptr, colptr,
                        tolscr2, startrow, startcol,
                        startMN, endMN, startPQ, endPQ,
                        startChunk, endChunk, chunksMN, 
                        D1, D2, D3,
                        ldX1, ldX2, ldX3,
                        sizeX1, sizeX2, sizeX3,
                        totalcalls, totalnintls);
                end = __rdtsc();
                signalled[mic_id] = 1;
                printf("Launched %d tasks on mic:%d in %lld cycles\n",
                    initialChunksMIC, mic_id, end - start);
                startChunk += initialChunksMIC;
            }

            int num_devices_active = num_devices;

            while(num_devices_active > 0)
            {
                for(mic_id = 0; mic_id < num_devices; mic_id++)
                {
                    if(signalled[mic_id] == 1)
                    {
                        int sig = _Offload_signaled (mic_id, finish_tag);
                        int chunksMIC = 0;

                        if (sig != 0)
                        {
                            signalled[mic_id] = 0;
                            #pragma omp critical
                            {
                                int remTasks = totalChunks - head;
                                chunksMIC = remTasks * mic_fraction;
                                if(chunksMIC < MIC_MIN_WORK_SIZE)
                                    chunksMIC = 0;
                                my_chunk = head;
                                head += chunksMIC;
                            }

                            if(chunksMIC > 0)
                            {
                                startChunk = my_chunk;
                                endChunk = startChunk + chunksMIC;
                                #pragma offload target(mic:mic_id) \
                                    nocopy(basis_mic, erd_mic) \
                                    in(shellptr, shellvalue, shellid, shellrid: length(0) REUSE) \
                                    in(f_startind, rowpos, colpos, rowptr, colptr: length(0) REUSE) \
                                    in(tolscr2, startrow, startcol) \
                                    in(startMN, endMN, startPQ, endPQ) \
                                    in(startChunk, endChunk, chunksMN) \
                                    in(D1, D2, D3: length(0) REUSE) \
                                    in(ldX1, ldX2, ldX3, sizeX1, sizeX2, sizeX3) \
                                    nocopy(totalcalls, totalnintls) \
                                    signal(finish_tag)
                                compute_block_of_chunks (basis_mic, erd_mic,
                                        shellptr, shellvalue,
                                        shellid, shellrid, f_startind,
                                        rowpos, colpos, rowptr, colptr,
                                        tolscr2, startrow, startcol,
                                        startMN, endMN, startPQ, endPQ,
                                        startChunk, endChunk, chunksMN, 
                                        D1, D2, D3,
                                        ldX1, ldX2, ldX3,
                                        sizeX1, sizeX2, sizeX3,
                                        totalcalls, totalnintls);
                                signalled[mic_id] = 1;
                                printf("Launched %d tasks on mic:%d\n", chunksMIC, mic_id);
                            }
                            else
                            {
                                num_devices_active--;
                            }

                        }
                    }
                }
            }
        }
        else
        {
            while(1)
            {
                #pragma omp critical
                {
                    my_chunk = head;
                    head++;
                }
                if(my_chunk >= totalChunks) break;

                int chunkIdMN = my_chunk / chunksMN;
                int chunkIdPQ = my_chunk % chunksMN;

                int startChunkMN =  startMN + chunkIdMN * CHUNK_SIZE;
                int endChunkMN   = startChunkMN + CHUNK_SIZE;
                if(endChunkMN > endMN) endChunkMN = endMN;
                int startChunkPQ =  startPQ + chunkIdPQ * CHUNK_SIZE;
                int endChunkPQ   = startChunkPQ + CHUNK_SIZE;
                if(endChunkPQ > endPQ) endChunkPQ = endPQ;
                compute_chunk (basis, erd,
                        shellptr, shellvalue,
                        shellid, shellrid, f_startind,
                        rowpos, colpos, rowptr, colptr,
                        tolscr2, startrow, startcol,
                        startChunkMN, endChunkMN, startChunkPQ, endChunkPQ,
                        D1, D2, D3,
                        J1, J2, K3,
                        ldX1, ldX2, ldX3,
                        sizeX1, sizeX2, sizeX3,
                        totalcalls, totalnintls, tid);
            }
        }
    } /* #pragma omp parallel */

}


void offload_reduce_mic (int num_devices,
                         double *F1_offload,
                         double *F2_offload,
                         double *F3_offload,
                         int sizeD1, int sizeD2, int sizeD3)
{
    for(int mic_id = 0; mic_id < num_devices; mic_id++)
    {
        #pragma offload target(mic:mic_id)\
            nocopy(pfock_mic)\
            out(F1_offload[0:sizeD1]: into(F1_offload[mic_id*sizeD1:sizeD1]))\
            out(F2_offload[0:sizeD2]: into(F2_offload[mic_id*sizeD2:sizeD2]))\
            out(F3_offload[0:sizeD3]: into(F3_offload[mic_id*sizeD3:sizeD3]))\
            signal(mic_id)
        {
            int num_F = pfock_mic.num_F;
            double **F1 = pfock_mic.F1;
            double **F2 = pfock_mic.F2;
            double **F3 = pfock_mic.F3;
            int sizeD1 = pfock_mic.sizeD1;
            int sizeD2 = pfock_mic.sizeD2;
            int sizeD3 = pfock_mic.sizeD3;
            #pragma omp parallel
            {
                // reduction on MIC
                for (int i = 1; i < num_F; i++)
                {
                    #pragma omp for
                    for (int j = 0; j < sizeD1; j++)
                    {
                        F1[0][j] += F1[i][j];
                    }
                    #pragma omp for
                    for (int j = 0; j < sizeD2; j++)
                    {
                        F2[0][j] += F2[i][j];
                    }
                    #pragma omp for 
                    for (int j = 0; j < sizeD3; j++)
                    {
                        F3[0][j] += F3[i][j];
                    }
                }
            }
        }
    }
}


void offload_reduce (int num_devices,
                     double *F1, double *F2, double *F3,
                     double *F1_offload, double *F2_offload, double *F3_offload,
                     int sizeD1, int sizeD2, int sizeD3)
{
    int i, j;

    for(i = 0; i < num_devices; i++)
    {
        #pragma omp parallel
        {
            #pragma omp for
            for (j = 0; j < sizeD1; j++)
            {
                F1[j + 0 * sizeD1] += F1_offload[j + i * sizeD1];
            }
            #pragma omp for
            for (j = 0; j < sizeD2; j++)
            {
                F2[j + 0 * sizeD2] += F2_offload[j + i * sizeD2];
            }
            #pragma omp for 
            for (j = 0; j < sizeD3; j++)
            {
                F3[j + 0 * sizeD3] += F3_offload[j + i * sizeD3];
            }
        }
    }
}


void offload_init (int nshells, int nnz,
                   int *shellptr, int *shellid, 
                   int *shellrid, double *shellvalue,
                   int *f_startind,
                   int *rowpos, int *colpos,
                   int *rowptr, int *colptr,
                   double *D1, double *D2, double *D3,
                   double *VD1, double *VD2, double *VD3,
                   double **_F1_offload,
                   double **_F2_offload,
                   double **_F3_offload,
                   int sizeD1, int sizeD2, int sizeD3,
                   int *_mic_numdevs,
                   int *_nthreads_mic)
{
    int mic_numdevs;
    int nthreads_mic;
    double *F1_offload;
    double *F2_offload;
    double *F3_offload;

    mic_numdevs = _Offload_number_of_devices();
    if (mic_numdevs <= 0)
    {
        *_mic_numdevs = 0;
        return;
    }

    #pragma offload target(mic:0) out(nthreads_mic)
    {
        nthreads_mic = omp_get_max_threads ();
    }
    for(int mic_id = 0; mic_id < mic_numdevs; mic_id++)
    {
        #pragma offload target(mic:mic_id) in(nthreads_mic)
        {
            omp_set_num_threads (nthreads_mic);
        }
    }
        
    F1_offload = 
        (double *) PFOCK_MALLOC (sizeof (double) * sizeD1 * mic_numdevs, 64);
    F2_offload =
        (double *) PFOCK_MALLOC (sizeof (double) * sizeD2 * mic_numdevs, 64);
    F3_offload =
        (double *) PFOCK_MALLOC (sizeof (double) * sizeD3 * mic_numdevs, 64);
    if (F1_offload == NULL ||
        F2_offload == NULL ||
        F3_offload == NULL)
    {
        *_mic_numdevs = 0;
        return;        
    }

    *_mic_numdevs = mic_numdevs;
    *_nthreads_mic = nthreads_mic;
    *_F1_offload = F1_offload;
    *_F2_offload = F2_offload;
    *_F3_offload = F3_offload;

    for(int mic_id = 0; mic_id < mic_numdevs; mic_id++)
    {
        #pragma offload target(mic: mic_id) \
            nocopy(pfock_mic)\
            in(sizeD1, sizeD2, sizeD3, nthreads_mic) \
            nocopy(D1: length(sizeD1) ALLOC) \
            nocopy(D2: length(sizeD2) ALLOC) \
            nocopy(D3: length(sizeD3) ALLOC) \
            nocopy(F1_offload: length(sizeD1) ALLOC) \
            nocopy(F2_offload: length(sizeD2) ALLOC) \
            nocopy(F3_offload: length(sizeD3) ALLOC) \
            in(shellptr: length(nshells + 1) ALLOC) \
            in(shellid: length(nnz) ALLOC)  \
            in(shellrid: length(nnz) ALLOC) \
            in(shellvalue: length(nnz) ALLOC) \
            in(f_startind: length(nshells + 1) ALLOC) \
            in(rowpos: length(nshells) ALLOC) \
            in(colpos: length(nshells) ALLOC) \
            in(rowptr: length(nnz) ALLOC) \
            in(colptr: length(nnz) ALLOC)
        {
            int num_F;
            pfock_mic.sizeD1 = sizeD1;
            pfock_mic.sizeD2 = sizeD2;
            pfock_mic.sizeD3 = sizeD3;
            pfock_mic.nthreads = nthreads_mic;
            num_F = pfock_mic.num_F = NUM_F_COPIES_MIC (nthreads_mic);
            pfock_mic.F1 = (double **)malloc (sizeof(double *) * num_F);
            pfock_mic.F2 = (double **)malloc (sizeof(double *) * num_F);
            pfock_mic.F3 = (double **)malloc (sizeof(double *) * num_F);
            assert (pfock_mic.F1 != NULL &&
                    pfock_mic.F2 != NULL &&
                    pfock_mic.F3 != NULL);
            pfock_mic.F1[0] = F1_offload;
            pfock_mic.F2[0] = F2_offload;
            pfock_mic.F3[0] = F3_offload;
            for (int i = 1; i < num_F; i++) 
            {
                pfock_mic.F1[i] = (double *)PFOCK_MALLOC (sizeof(double) * sizeD1, 64);
                pfock_mic.F2[i] = (double *)PFOCK_MALLOC (sizeof(double) * sizeD2, 64);
                pfock_mic.F3[i] = (double *)PFOCK_MALLOC (sizeof(double) * sizeD3, 64);
                assert (pfock_mic.F1[i] != NULL &&
                        pfock_mic.F2[i] != NULL &&
                        pfock_mic.F3[i] != NULL);
            }           
        }
    }
}


void offload_deinit (int mic_numdevs,
                     int *shellptr, int *shellid, 
                     int *shellrid, double *shellvalue,
                     int *f_startind,
                     int *rowpos, int *colpos,
                     int *rowptr, int *colptr,
                     double *D1, double *D2, double *D3,
                     double *VD1, double *VD2, double *VD3,
                     double *F1_offload,
                     double *F2_offload,
                     double *F3_offload)                     
{
    for(int mic_id = 0; mic_id < mic_numdevs; mic_id++)
    {
        #pragma offload target(mic: mic_id) \
            nocopy(pfock_mic)\
            nocopy(D1: length(0) FREE) \
            nocopy(D2: length(0) FREE) \
            nocopy(D3: length(0) FREE) \
            nocopy(F1_offload: length(0) FREE) \
            nocopy(F2_offload: length(0) FREE) \
            nocopy(F3_offload: length(0) FREE) \
            nocopy(shellptr: length(0) FREE) \
            nocopy(shellid: length(0) FREE)  \
            nocopy(shellrid: length(0) FREE) \
            nocopy(shellvalue: length(0) FREE) \
            nocopy(f_startind: length(0) FREE) \
            nocopy(rowpos: length(0) FREE) \
            nocopy(colpos: length(0) FREE) \
            nocopy(rowptr: length(0) FREE) \
            nocopy(colptr: length(0) FREE)
        {
            int num_F;
            num_F = pfock_mic.num_F;
            for (int i = 1; i < num_F; i++) 
            {
                PFOCK_FREE (pfock_mic.F1[i]);
                PFOCK_FREE (pfock_mic.F2[i]);
                PFOCK_FREE (pfock_mic.F3[i]);
            }
            free (pfock_mic.F1);
            free (pfock_mic.F2);
            free (pfock_mic.F3);
        }
    }
    
    PFOCK_FREE (F1_offload);
    PFOCK_FREE (F2_offload);
    PFOCK_FREE (F3_offload);
}


void offload_copy_D (int num_devices, double *D, int sizeD)
{
    for (int mic_id = 0; mic_id < num_devices; mic_id++)
    {
        #pragma offload_transfer target(mic:mic_id) \
                in(D: length(sizeD) REUSE)
    }
}


void offload_wait_mic (int num_devices)
{
    for(int mic_id = 0; mic_id < num_devices; mic_id++)
    {
        #pragma offload_wait target(mic:mic_id) wait(mic_id)
    }
}
