#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__HRR_TRANSFORM */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation performs a HRR transformation on a */
/*                batch of contracted cartesian gaussian integrals: */
/*                                      nxyzet */
/*                     y (m,nxyzab)  =   sum   x (m,i) * rot (i,nxyzab) */
/*                                       i=1 */
/*                where rot is the HRR transformation matrix. Due to */
/*                the very sparse nature of this matrix, only those */
/*                i indices in the summation are addressed which */
/*                correspond to nonzero HRR transformation matrix */
/*                elements. */
/*                  Input: */
/*                     M          =  # of elements not involved in the */
/*                                   transformation (invariant indices) */
/*                     NROW       =  maximum # of nonzero row elements */
/*                                   per column in transformation matrix */
/*                     NXYZET     =  dimension of the cartesian e0-part */
/*                     NXYZAB     =  dimension of the resulting cartesian */
/*                                   ab-part */
/*                     NXYZA      =  dimension of the cartesian part due */
/*                                   to the a-shell */
/*                     NXYZB      =  dimension of the cartesian part due */
/*                                   to the b-shell */
/*                     LROW (N)   =  # of nonzero entries in column N of */
/*                                   ROT and ROW matrix. N ranges from */
/*                                   1 to NXYZB */
/*                     ROW (I,N)  =  I-th nonzero row label of column N */
/*                                   in ROT matrix. N ranges from 1 to */
/*                                   NXYZAB */
/*                     ROT (I,N)  =  I-th nonzero HRR transformation */
/*                                   matrix element of column N. N ranges */
/*                                   from 1 to NXYZB */
/*                     X          =  batch of untransformed integrals */
/*                                   (m,e0) */
/*                  Output: */
/*                     Y          =  batch of HRR transformed integrals */
/*                                   (m,ab) */
/* ------------------------------------------------------------------------ */
int erd__hrr_transform (int m, int nrow,
                        int nxyza, int nxyzb,
                        int *lrow, int *row,
                        double *rot, double *x, double *y)
{
    int a;
    int b;
    int i;
    int j;
    int n;
    double rot1;
    int mrow;
    int xcol1;

/*             ...perform the HRR transformation. One of the main */
/*                properties of this transformation is that the */
/*                last nonzero element of the HRR transformation */
/*                matrix is always equal to 1. Hence we can skip */
/*                the multiplication with that element. */
/*                Use basic row grouping of the transformation */
/*                to improve cache line reusing. */
    n = 0;
    for (b = 0; b < nxyzb; ++b)
    {
        mrow = lrow[b];
        for (a = 0; a < nxyza; ++a)
        {
            for (j = 0; j < m; ++j)
            {
                y[j + n * m] = 0.0;
            }           
            for (i = 0; i < mrow; i++)
            {
                xcol1 = row[i + n * nrow] - 1;
                rot1 = rot[i + b * nrow];
                for (j = 0; j < m; ++j)
                {
                    y[j + n * m] += rot1 * x[j + xcol1 * m];
                }
            }
            ++n;
        }      
    }

    return 0;
}
