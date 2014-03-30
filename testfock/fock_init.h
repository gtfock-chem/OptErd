#ifndef __FOCK_INIT_H__
#define __FOCK_INIT_H__


#include "CInt.h"


#define TOLSRC 1e-10


void schwartz_screening (BasisSet_t basis, int **shellptrOut,
                         int **shellidOut, int **shellridOut,
                         double **shellvalueOut, int *nnzOut);

void create_buffers (int nshells, int nnz,
                     int *shellptr, int *shellid, int *f_startind,
                     int startM, int endM, int startP, int endP,
                     int **rowposOut, int **colposOut,
                     int **rowptrOut, int **colptrOut,
                     int *rowfuncs, int *colfuncs,
                     int *rowsize, int *colsize);

void fock_task (BasisSet_t basis, ERD_t erd,
                int *shellptr, double *shellvalue,
                int *shellid, int *shellrid, int *f_startind,
                int *rowpos, int *colpos, int *rowptr, int *colptr,
                double tolscr2,
                int startrow, int startcol,
                int startM, int endM, int startP, int endP,
                double *D1, double *D2, double *D3,
                double *F1, double *F2, double *F3,
                int ldX1, int ldX2, int ldX3,
                int sizeX1, int sizeX2, int sizeX3,
                double *nsq, double *nitl);


#endif /* __FOCK_INIT_H__ */
