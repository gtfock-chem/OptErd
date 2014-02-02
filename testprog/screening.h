#ifndef __SCREENING_H__
#define __SCREENING_H__


#if defined(OPTERD_TEST_REFERENCE)
#include "../legacy/include/CInt.h"
#else
#include "../include/CInt.h"
#endif


#define TOLSRC 1e-10


void schwartz_screening (BasisSet_t basis, int **shellptr,
                         int **shellid, int **shellrid,
                         double **shellvalue, int *nnz);


#endif /* __SCREENING_H__ */
