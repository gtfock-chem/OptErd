#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <macdecls.h>
#include <sys/time.h>

#include "CInt.h"


static void update_F (double *integrals, int dimM, int dimN, int dimP, int dimQ,
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

    for (iM = 0; iM < dimM; iM++)
    {
        for (iN = 0; iN < dimN; iN++)
        {
            for (iP = 0; iP < dimP; iP++)
            {
                for (iQ = 0; iQ < dimQ; iQ++)
                {
                   I = integrals[iM + dimM * (iN + dimN * (iP + dimP * iQ))];                    
                    
                    // F(m, n) += D(p, q) * 2 * I(m, n, p, q)
                    // F(n, m) += D(p, q) * 2 * I(n, m, p, q)
                    // F(m, n) += D(q, p) * 2 * I(m, n, q, p)
                    // F(n, m) += D(q, p) * 2 * I(n, m, q, p)
                    J_MN[iM * ldX1 + iN] += 1.0 * (1 + flag1 + flag2 + flag4) *
                            D_PQ[iP * ldX2 + iQ] * I;
                    
                    // F(p, q) += D(m, n) * 2 * I(p, q, m, n)
                    // F(p, q) += D(n, m) * 2 * I(p, q, n, m)
                    // F(q, p) += D(m, n) * 2 * I(q, p, m, n)
                    // F(q, p) += D(n, m) * 2 * I(q, p, n, m)
                    J_PQ[iP * ldX2 + iQ] += 1.0 * (flag3 + flag5 + flag6 + flag7) *
                            D_MN[iM * ldX1 + iN] * I;
                    // F(m, p) -= D(n, q) * I(m, n, p, q)
                    // F(p, m) -= D(q, n) * I(p, q, m, n)
                    K_MP[iM * ldX3 + iP] -= (1 + flag3) *
                            0.5 * D_NQ[iN * ldX3 + iQ] * I;
                    // F(n, p) -= D(m, q) * I(n, m, p, q)
                    // F(p, n) -= D(q, m) * I(p, q, n, m)
                    K_NP[iN * ldX3 + iP] -= (flag1 + flag5) *
                            0.5 * D_MQ[iM * ldX3 + iQ] * I;
                    // F(m, q) -= D(n, p) * I(m, n, q, p)
                    // F(q, m) -= D(p, n) * I(q, p, m, n)
                    K_MQ[iM * ldX3 + iQ] -= (flag2 + flag6) *
                            0.5 * D_NP[iN * ldX3 + iP] * I;

                    // F(n, q) -= D(m, p) * I(n, m, q, p)
                    // F(q, n) -= D(p, m) * I(q, p, n, m)
                    K_NQ[iN * ldX3 + iQ] -= (flag4 + flag7) *
                           0.5 * D_MP[iM * ldX3 + iP] * I;
                }
            }
        }
    }
}


void compute_task (BasisSet_t basis, ERD_t erd,
                   int *shellptr, double *shellvalue,
                   int *shellid, int *shellrid, int *f_startind,
                   int *rowpos, int *colpos, int *rowptr, int *colptr,
                   double tolscr2, int startrow, int startcol,
                   int startM, int endM, int startP, int endP,
                   double *D1, double *D2, double *D3,
                   double *F1, double *F2, double *F3,
                   int ldX1, int ldX2, int ldX3,
                   int sizeX1, int sizeX2, int sizeX3)
{
    int startMN;
    int endMN;
    int startPQ;
    int endPQ;
    int sizeX1;
    int sizeX2;
    int sizeX3;
    int nsq;
    int nitl;
    
    startMN = shellptr[startM];
    endMN = shellptr[endM + 1];
    startPQ = shellptr[startP];
    endPQ = shellptr[endP + 1];

    #pragma omp parallel
    {
        int i;
        int j;
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
        int nt;
   
        nt = omp_get_thread_num ();
        double *J1 = &(F1[nt * sizeX1]);
        double *J2 = &(F2[nt * sizeX2]);
        double *K3 = &(F3[nt * sizeX3]);
        mynsq = 0;
        mynitl = 0;
        
        #pragma omp for schedule(dynamic)
        for (i = startMN; i < endMN; i++)
        {
            M = shellrid[i];
            N = shellid[i];
            value1 = shellvalue[i];            
            dimM = f_startind[M + 1] - f_startind[M];
            dimN = f_startind[N + 1] - f_startind[N];
            iX3M = rowpos[M]; 
            iX1M = f_startind[M] - f_startind[startrow];           
            iX1N = iX3N = rowptr[i];
            iMN = iX1M * ldX1+ iX1N;
            flag1 = (value1 < 0.0) ? 1 : 0;
            for (j = startPQ; j < endPQ; j++)
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
                dimQ =  f_startind[Q + 1] - f_startind[Q];
                iX3P = colpos[P];
                iX2P = f_startind[P] - f_startind[startcol];
                iMP = iX3M * ldX3 + iX3P;
                iNP = iX3N * ldX3 + iX3P;              
                iX3Q = iX2Q = colptr[j];
                iPQ = iX2P * ldX2+ iX2Q;
                iMQ = iX3M * ldX3 + iX3Q;
                iNQ = iX3N * ldX3 + iX3Q;                            
                flag3 = (M == P && Q == N) ? 0 : 1;                    
                flag2 = (value2 < 0.0) ? 1 : 0;
                if (fabs(value1 * value2) >= tolscr2)
                {                     
                    CInt_computeShellQuartet (basis, erd, nt,
                                              M, N, P, Q, &integrals, &nints);
                    if (nints != 0)
                    {
                        update_F (integrals, dimM, dimN, dimP, dimQ,
                                  flag1, flag2, flag3,
                                  iMN, iPQ, iMP, iNP, iMQ, iNQ,
                                  D1, D2, D3, J1, J2, K3,
                                  ldX1, ldX2, ldX3);
                    }
                }
            }
        }
    } /* #pragma omp parallel */
}
