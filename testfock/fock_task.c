#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <sys/time.h>

#include "CInt.h"
#include "fock_init.h"

#include <offload.h>
#define CHUNK_SIZE 8
#define PFD0 1
#define PFD1 12


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <immintrin.h>

static inline void update_F (double *integrals, int dimM, int dimN, int dimP, int dimQ,
                      int flag1, int flag2, int flag3,
                      int iMN, int iPQ, int iMP, int iNP, int iMQ, int iNQ,
                      double *D1, double *D2, double *D3,
                      double *J1, double *J2, double *K3, int ldX1, int ldX2, int ldX3)
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
                    //atomic_dadd(&J_PQ[iPQ], vPQ)

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
                    //atomic_dadd(&K_MQ[iMQ], -vMQ)

                    // F(n, q) -= D(m, p) * I(n, m, q, p)
                    // F(q, n) -= D(p, m) * I(q, p, n, m)
                    double vNQ = (flag4 + flag7) * 0.5 * D_MP[iMP] * I;
#ifdef __MIC__
                    #pragma omp atomic update
#endif
                    K_NQ[iNQ] -= vNQ;
                    //atomic_dadd(&K_NQ[iNQ], -vNQ)
                }
#ifdef __MIC__
                #pragma omp atomic update
#endif
                K_MP[iMP] += k_MP;
#ifdef __MIC__
                #pragma omp atomic update
#endif
                K_NP[iNP] += k_NP;
                //atomic_dadd(&K_MP[iMP], k_MP)
                //atomic_dadd(&K_NP[iNP], k_NP)
            }
#ifdef __MIC__
            #pragma omp atomic update
#endif
            J_MN[iMN] += j_MN;
            //atomic_dadd(&J_MN[iMN],  j_MN)
        }
    }
}


void compute_chunk (BasisSet_t basis, ERD_t erd,
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
            //totalcalls[nt * 64] += 1;
            if (fabs (value1 * value2) >= tolscr2)
            {
                CInt_computeShellQuartet (basis, erd, nt,
                        M, N, P, Q, &integrals, &nints);
                //totalnintls[nt * 64] += nints;
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



void compute_block_of_chunks (BasisSet_t basis, ERD_t erd,
                   int *shellptr, double *shellvalue,
                   int *shellid, int *shellrid, int *f_startind,
                   int *rowpos, int *colpos, int *rowptr, int *colptr,
                   double tolscr2, int startrow, int startcol,
                   int startMN, int endMN, int startPQ, int endPQ,
                   int startChunk, int endChunk, int chunksMN,
                   double *D1, double *D2, double *D3,
                   double *F1, double *F2, double *F3,
                   int ldX1, int ldX2, int ldX3,
                   int sizeX1, int sizeX2, int sizeX3,
                   double *totalcalls, double *totalnintls,
                   int nthreads_mic)
{
    int k;

#pragma omp parallel
    {
        int tid;
        
        tid = omp_get_thread_num ();
        double *J1 = &(F1[(tid >> 2) * sizeX1]);
        double *J2 = &(F2[(tid >> 2) * sizeX2]);
        double *K3 = &(F3[(tid >> 2) * sizeX3]);

        #pragma omp for schedule(dynamic)
        for(k = startChunk; k < endChunk; k++)
        {
            int chunkIdMN = k / chunksMN;
            int chunkIdPQ = k % chunksMN;

            int startChunkMN = startMN + chunkIdMN * CHUNK_SIZE;
            int endChunkMN   = startChunkMN + CHUNK_SIZE;
            int startChunkPQ = startPQ + chunkIdPQ * CHUNK_SIZE;
            int endChunkPQ   = startChunkPQ + CHUNK_SIZE;

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

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

void compute_task_cpu (int num_devices,
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
                   double *totalcalls, double *totalnintls,
                   int nthreads, int nthreads_mic)
{
    int startMN;
    int endMN;
    int startPQ;
    int endPQ;
    startMN = shellptr[startM];
    endMN = shellptr[endM + 1];
    startPQ = shellptr[startP];
    endPQ = shellptr[endP + 1];

    //printf("startMN = %d, endMN = %d\n", startMN, endMN);
    //printf("startPQ = %d, endPQ = %d\n", startPQ, endPQ);

    int chunksMN = (endMN - startMN) / CHUNK_SIZE;
    int chunksPQ = (endPQ - startPQ) / CHUNK_SIZE;
    int totalChunks = chunksMN * chunksPQ;

    int initialChunksMIC = totalChunks * mic_fraction;
/*
    int sx;
    for(sx = startMN; sx < endMN; sx++)
    {
        printf("M = %d, N = %d\n", shellrid[sx], shellid[sx]);
    }
    for(sx = startPQ; sx < endPQ; sx++)
    {
        printf("P = %d, Q = %d\n", shellrid[sx], shellid[sx]);
    }
    */
#if 0
                compute_block_of_chunks (basis, erd,
                        shellptr, shellvalue,
                        shellid, shellrid, f_startind,
                        rowpos, colpos, rowptr, colptr,
                        tolscr2, startrow, startcol,
                        startMN, endMN, startPQ, endPQ,
                        0, totalChunks, chunksMN, 
                        D1, D2, D3,
                        F1, F2, F3,
                        ldX1, ldX2, ldX3,
                        sizeX1, sizeX2, sizeX3,
                        totalcalls, totalnintls);
#endif
#if 1
    int head = num_devices * initialChunksMIC;
    #pragma omp parallel
    {
        int my_chunk;
        int tid;
        
        tid = omp_get_thread_num ();
        double *J1 = &(F1[tid * sizeX1]);
        double *J2 = &(F2[tid * sizeX2]);
        double *K3 = &(F3[tid * sizeX3]);


        if(tid == 0)
        {
#ifdef __INTEL_OFFLOAD
#if 1
            int signalled[num_devices];
            int mic_id;
            int startChunk = 0;
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
                nocopy(F1_mic, F2_mic, F3_mic) \
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
                        F1_mic, F2_mic, F3_mic,
                        ldX1, ldX2, ldX3,
                        sizeX1, sizeX2, sizeX3,
                        totalcalls, totalnintls, nthreads_mic);
                end = __rdtsc();
                signalled[mic_id] = 0;
                printf("Processed %d tasks on mic:%d in %lld cycles\n", initialChunksMIC, mic_id, end - start);
                startChunk += initialChunksMIC;
            }
#if 1
            int num_devices_done = 0;

            while(num_devices_done < num_devices)
            {
                for(mic_id = 0; mic_id < num_devices; mic_id++)
                {
                    int sig = _Offload_signaled (mic_id, finish_tag);
                    //printf("%d) signaled = %d\n", mic_id, sig);
                    if (sig != 0 && signalled[mic_id] == 0)
                    {
                        signalled[mic_id] = 1;
                        num_devices_done++;
                    }
                }
            }
#endif
#if 0
            while(1)
            {
#pragma omp critical
                {
                    my_chunk = head;
                    head+=MIC_WORK_SIZE;
                }
                if(my_chunk >= totalChunks) break;

                int chunkIdMN = my_chunk / chunksMN;
                int chunkIdPQ = my_chunk % chunksMN;

                int startChunkMN = startMN + chunkIdMN * CHUNK_SIZE;
                int endChunkMN   = startChunkMN + CHUNK_SIZE;
                int startChunkPQ = startPQ + chunkIdPQ * CHUNK_SIZE;
                int endChunkPQ   = startChunkPQ + CHUNK_SIZE;

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
#endif
#endif
#endif
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

                int startChunkMN = chunkIdMN * CHUNK_SIZE;
                int endChunkMN   = startChunkMN + CHUNK_SIZE;
                int startChunkPQ = chunkIdPQ * CHUNK_SIZE;
                int endChunkPQ   = startChunkPQ + CHUNK_SIZE;

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
        //printf("tid = %d\n", nt);
    } /* #pragma omp parallel */
#endif

}


void compute_task_hetero (int num_devices,
                   BasisSet_t basis, ERD_t erd,
                   int *shellptr, double *shellvalue,
                   int *shellid, int *shellrid, int *f_startind,
                   int *rowpos, int *colpos, int *rowptr, int *colptr,
                   double tolscr2, int startrow, int startcol,
                   int startM, int endM, int startP, int endP,
                   double *D1, double *D2, double *D3,
                   double *F1, double *F2, double *F3,
                   double *F1_hetero, double *F2_hetero, double *F3_hetero,
                   int ldX1, int ldX2, int ldX3,
                   int sizeX1, int sizeX2, int sizeX3, double mic_fraction,
                   double *totalcalls, double *totalnintls)
{
    //int mic_id;
    long start, end;

    //int sizetask = endM - startM + 1;
    //int sizetask_MIC = sizetask * mic_fraction;
    //sizetask_MIC = ((sizetask_MIC + 63)/64) * 64;
    //int curStartP = startP;

    start = __rdtsc();
    //int sizetaskP = endP - curStartP + 1;
    printf("Launching tasks on CPU and MIC\n");
    compute_task_cpu (num_devices, basis, erd, shellptr,
            shellvalue, shellid, shellrid,
            f_startind, rowpos, colpos, rowptr, colptr,
            tolscr2, startrow, startcol,
            startM, endM,
            startP, endP,
            D1, D2, D3, F1, F2, F3,
            ldX1, ldX2, ldX3, sizeX1, sizeX2, sizeX3,
            mic_fraction,
            totalcalls, totalnintls,
            24, 240);
    end = __rdtsc();
    printf("Finished tasks on CPU+%dMIC in %lld cycles\n", num_devices, end - start);

}

void reduce_F(int num_devices,
              double *F1, double *F2, double *F3,
              double *F1_hetero, double *F2_hetero, double *F3_hetero,
              int sizeD1, int sizeD2, int sizeD3,
              int nthreads, int nthreads_mic)
{
    int dummy_tag = 0;
    int *finish_tag = &dummy_tag;

    long start, end;
#ifdef __INTEL_OFFLOAD
    int mic_id;
    for(mic_id = 0; mic_id < num_devices; mic_id++)
    {
        printf("Reducing on mic:%d\n", mic_id);
        start = __rdtsc();
#pragma offload target(mic:mic_id) \
        in(sizeD1, sizeD2, sizeD3, nthreads_mic) \
        nocopy(F1_mic, F2_mic, F3_mic) \
        signal(finish_tag)
        {
            int j;
#pragma omp parallel
            {
                int i;
                // reduction
                for (i = 1; i < nthreads_mic; i++)
                {
#pragma omp for
                    for (j = 0; j < sizeD1; j++)
                    {
                        //_mm_prefetch((const char* ) &F1_mic[j + (i + PFD0)* sizeD1], _MM_HINT_T0);
                        //_mm_prefetch((const char* ) &F1_mic[j + (i + PFD1)* sizeD1], _MM_HINT_T1);
                        F1_mic[j + 0 * sizeD1] += F1_mic[j + i * sizeD1];
                    }
#pragma omp for
                    for (j = 0; j < sizeD2; j++)
                    {
                        //_mm_prefetch((const char* ) &F2_mic[j + (i + PFD0)* sizeD2], _MM_HINT_T0);
                        //_mm_prefetch((const char* ) &F2_mic[j + (i + PFD1)* sizeD2], _MM_HINT_T1);
                        F2_mic[j + 0 * sizeD2] += F2_mic[j + i * sizeD2];
                    }
#pragma omp for 
                    for (j = 0; j < sizeD3; j++)
                    {
                        //_mm_prefetch((const char* ) &F3_mic[j + (i + PFD0)* sizeD3], _MM_HINT_T0);
                        //_mm_prefetch((const char* ) &F3_mic[j + (i + PFD1)* sizeD3], _MM_HINT_T1);
                        F3_mic[j + 0 * sizeD3] += F3_mic[j + i * sizeD3];
                    }
                }
            }
        }
        end = __rdtsc();
        printf("mic:%d reduce cycles = %lld\n", mic_id, end - start);
    }
#endif

    int i, j;
    printf("sizeD1 = %d, sizeD2 = %d, sizeD3 = %d\n", sizeD1, sizeD2, sizeD3);
    start = __rdtsc();
    // reduction
    for (i = 1; i < nthreads; i++)
    {
#pragma omp parallel
        {
#pragma omp for
            for (j = 0; j < sizeD1; j++)
            {
                F1[j + 0 * sizeD1] += F1[j + i * sizeD1];
            }
#pragma omp for
            for (j = 0; j < sizeD2; j++)
            {
                F2[j + 0 * sizeD2] += F2[j + i * sizeD2];
            }
#pragma omp for 
            for (j = 0; j < sizeD3; j++)
            {
                F3[j + 0 * sizeD3] += F3[j + i * sizeD3];
            }
        }
    }
    end = __rdtsc();
    printf("CPU reduce cycles = %lld\n", end - start);

#ifdef __INTEL_OFFLOAD
    for(mic_id = 0; mic_id < num_devices; mic_id++)
    {
#pragma offload_wait target(mic:mic_id) wait(finish_tag)
    }
#endif

}


void reduce_F_CPU_MIC(int num_devices,
                      double *F1, double *F2, double *F3,
                      double *F1_hetero, double *F2_hetero, double *F3_hetero,
                      int sizeD1, int sizeD2, int sizeD3)
{
    int dummy_tag = 0;
    int *finish_tag = &dummy_tag;
    int mic_id;
    int i, j;

    long start, end;
#ifdef __INTEL_OFFLOAD
    for(mic_id = 0; mic_id < num_devices; mic_id++)
    {
        start = __rdtsc();
#pragma offload target(mic:mic_id) \
        in(sizeD1, sizeD2, sizeD3) \
        nocopy(F1_mic, F2_mic, F3_mic) \
        out(F1_hetero[0:sizeD1]: into(F1_hetero[mic_id * sizeD1:sizeD1])) \
        out(F2_hetero[0:sizeD2]: into(F2_hetero[mic_id * sizeD2:sizeD2])) \
        out(F3_hetero[0:sizeD3]: into(F3_hetero[mic_id * sizeD3:sizeD3])) \
        signal(finish_tag)
        {
            memcpy(F1_hetero, F1_mic, sizeD1 * sizeof(double));
            memcpy(F2_hetero, F2_mic, sizeD2 * sizeof(double));
            memcpy(F3_hetero, F3_mic, sizeD3 * sizeof(double));
        }
        end = __rdtsc();
        printf("Copy F to CPU: mic:%d, %lld cycles\n", mic_id, end - start);
    }
    for(mic_id = 0; mic_id < num_devices; mic_id++)
    {
#pragma offload_wait target(mic:mic_id) wait(finish_tag)
    }
#endif
    printf("F1[0]: cpu: %E, mic0:%E, mic1:%E\n", F1[0], F1_hetero[0], F1_hetero[sizeD1]);
    start = __rdtsc();
    for(i = 0; i < num_devices; i++)
    {
#pragma omp parallel
        {
#pragma omp for
            for (j = 0; j < sizeD1; j++)
            {
                F1[j + 0 * sizeD1] += F1_hetero[j + i * sizeD1];
            }
#pragma omp for
            for (j = 0; j < sizeD2; j++)
            {
                F2[j + 0 * sizeD2] += F2_hetero[j + i * sizeD2];
            }
#pragma omp for 
            for (j = 0; j < sizeD3; j++)
            {
                F3[j + 0 * sizeD3] += F3_hetero[j + i * sizeD3];
            }
        }
    }
    end = __rdtsc();
    printf("Reduced F on CPU in %lld cycles\n", end - start);
}
