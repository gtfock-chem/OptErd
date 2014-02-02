#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


#pragma offload_attribute(push, target(mic))

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__SPHERICAL_TRANSFORM */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation performs a simple cartesian to spherical */
/*                transformation step on a batch of contracted cartesian */
/*                gaussian integrals: */
/*                         y (m,r) = sum  mat (i,r) * x (m,i) */
/*                                    i */
/*                where i is an element of the cartesian basis, r an */
/*                element of the spherical basis, mat is the cartesian */
/*                -> spherical transformation matrix and m are the bra */
/*                elements not involved. Ordering of the cartesian */
/*                monomial basis is assumed to be: */
/*                         e f g     a b c */
/*                        X Y Z  >  X Y Z   , if e > a */
/*                                            for e=a if f > b */
/*                                            for e=a and f=b if g > c */
/*                  Input: */
/*                    M           =  # of elements not involved in the */
/*                                   transformation (invariant indices) */
/*                    NROW        =  maximum # of nonzero rows in the */
/*                                   transformation matrix */
/*                    NXYZ        =  # of cartesian monomials xyz for */
/*                                   the shell to be transformed */
/*                    NRY         =  # of spherical functions ry for */
/*                                   the shell to be transformed */
/*                    LROW (R)    =  # of xyz-monomials contributing */
/*                                   to the R-th ry-component */
/*                    ROW (I,R)   =  I-th xyz-monomial row index */
/*                                   containing nonzero contribution to */
/*                                   the R-th ry-component */
/*                    ROT (I,R)   =  I-th nonzero xyz-monomial to R-th */
/*                                   ry-component transformation matrix */
/*                                   element */
/*                    X           =  input batch of cartesian integrals */
/*                  Output: */
/*                    Y           =  output batch of spherical integrals */
/* ------------------------------------------------------------------------ */
int erd__spherical_transform (
			      int m, int nrow, int nry,
                              int *lrow, int *row, double *rot,
                              double *x, double *y)
{   
    int i;
    int n;
    int r;
    double rot1;
    int mrow;
    int xcol1;

/*             ...perform the cartesian -> spherical transformation. */
/*                Use basic row grouping of the transformation */
/*                to improve cache line reusing. */
    for (r = 0; r < nry; ++r)
    {
        mrow = lrow[r];

        for (n = 0; n < m; ++n)
        {
            y[r * m + n] = 0.0;
        }
        
        for (i = 0; i < mrow; ++i)
        {
            xcol1 = row[r * nrow + i] - 1;
            rot1 = rot[r * nrow + i];
            for (n = 0; n < m; ++n)
            {
                y[r * m + n] += rot1 * x[xcol1 * m + n];
            }
        }
    }

    return 0;
}

#pragma offload_attribute(pop)
