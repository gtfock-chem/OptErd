#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__NORMALIZE_CARTESIAN */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation normalizes a cartesian monomial part */
/*                of a batch of contracted cartesian gaussian integrals. */
/*                The normalization factors are xyz-exponent dependent */
/*                and are given for a specific monomial as: */
/*                                    ______________________ */
/*                     l m n         /       2^(2L) */
/*                    x y z -->     / ----------------------- */
/*                                \/ (2l-1)!!(2m-1)!!(2n-1)!! */
/*                where L = l+m+n. The best way to deal with these */
/*                factors for a complete set of monomials for fixed L */
/*                is to split up each factor into its l-,m- and n- */
/*                component: */
/*                       _______        _______        _______ */
/*                      / 2^(2l)       / 2^(2m)       / 2^(2n) */
/*                     / -------  *   / -------  *   / ------- */
/*                   \/ (2l-1)!!    \/ (2m-1)!!    \/ (2n-1)!! */
/*                These factors are passed in argument as NORM. */
/*                  Input: */
/*                    M           =  # of elements not involved in the */
/*                                   normalization (invariant indices) */
/*                    NXYZ        =  # of monomials for shell L */
/*                    L           =  the shell type for which the */
/*                                   normalization will be done */
/*                    NORM        =  individual normalization factors */
/*                                   for each monomial exponent */
/*                    BATCH       =  batch of unnormalized integrals */
/*                  Output: */
/*                    BATCH       =  batch of normalized integrals */
/* ------------------------------------------------------------------------ */
int erd__normalize_cartesian (int m, int l,
                              double *norm, double *batch)
{
    int batch_offset;

    int i, n, x, y, z, ybeg;
    double xnorm, scalar;

    batch_offset = 1 + m * 1;
    batch -= batch_offset;

    n = 0;
    for (x = l; x >= 0; x--)
    {
        xnorm = norm[x];
        ybeg = l - x;
        for (y = ybeg; y >= 0; --y)
        {
            z = ybeg - y;
            scalar = xnorm * norm[y] * norm[z];
            ++n;
            for (i = 1; i <= m; ++i)
            {
                batch[i + n * m] =
                    scalar * batch[i + n * m];
            }
        }
    }

    return 0;
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
