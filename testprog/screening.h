#ifndef __SCREENING_H__
#define __SCREENING_H__


#include "CInt.h"


#define TOLSRC 1e-10


void schwartz_screening (BasisSet_t basis, int **shellptr,
                         int **shellid, int **shellrid,
                         double **shellvalue);


#endif /* __SCREENING_H__ */
