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
                        int nxyzet, int nxyzab,
                        int nxyza, int nxyzb,
                        int *lrow, int *row,
                        double *rot, double *x, double *y)
{
    int row_offset, x_offset, y_offset, rot_offset;

    int a, b, i, j, n;
    double rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8;
    int mrow, xcol1, xcol2, xcol3, xcol4, xcol5, xcol6, xcol7,
        xcol8, nleft, nstep;

/*             ...perform the HRR transformation. One of the main */
/*                properties of this transformation is that the */
/*                last nonzero element of the HRR transformation */
/*                matrix is always equal to 1. Hence we can skip */
/*                the multiplication with that element. */
/*                Use basic row grouping of the transformation */
/*                to improve cache line reusing. */
    x_offset = 1 + m * 1;
    x -= x_offset;
    y_offset = 1 + m * 1;
    y -= y_offset;
    row_offset = 1 + nrow * 1;
    row -= row_offset;
    rot_offset = 1 + nrow * 1;
    rot -= rot_offset;
    --lrow;

    n = 1;
    for (b = 1; b <= nxyzb; ++b)
    {
        mrow = lrow[b];
        switch (MIN(mrow, 9))
        {
        case 1:
            goto L1;
        case 2:
            goto L2;
        case 3:
            goto L3;
        case 4:
            goto L4;
        case 5:
            goto L5;
        case 6:
            goto L6;
        case 7:
            goto L7;
        case 8:
            goto L8;
        case 9:
            goto L9;
        }

      L1:
        for (a = 1; a <= nxyza; ++a)
        {
            xcol1 = row[n * nrow + 1];
            for (j = 1; j <= m; ++j)
            {
                y[j + n * m] = x[j + xcol1 * m];
            }
            ++n;
        }
        goto L100;


      L2:
        for (a = 1; a <= nxyza; ++a)
        {
            xcol1 = row[n * nrow + 1];
            xcol2 = row[n * nrow + 2];
            rot1 = rot[b * nrow + 1];
            for (j = 1; j <= m; ++j)
            {
                y[j + n * m] = rot1 * x[j + xcol1 * m] +
                    x[j + xcol2 * m];
            }
            ++n;
        }
        goto L100;

      L3:
        for (a = 1; a <= nxyza; ++a)
        {
            xcol1 = row[n * nrow + 1];
            xcol2 = row[n * nrow + 2];
            xcol3 = row[n * nrow + 3];
            rot1 = rot[b * nrow + 1];
            rot2 = rot[b * nrow + 2];
            for (j = 1; j <= m; ++j)
            {
                y[j + n * m] = rot1 * x[j + xcol1 * m] +
                    rot2 * x[j + xcol2 * m] +
                    x[j + xcol3 * m];
            }
            ++n;
        }
        goto L100;

      L4:
        for (a = 1; a <= nxyza; ++a)
        {
            xcol1 = row[n * nrow + 1];
            xcol2 = row[n * nrow + 2];
            xcol3 = row[n * nrow + 3];
            xcol4 = row[n * nrow + 4];
            rot1 = rot[b * nrow + 1];
            rot2 = rot[b * nrow + 2];
            rot3 = rot[b * nrow + 3];
            for (j = 1; j <= m; ++j)
            {
                y[j + n * m] = rot1 * x[j + xcol1 * m] +
                    rot2 * x[j + xcol2 * m] +
                    rot3 * x[j + xcol3 * m] +
                    x[j + xcol4 * m];
            }
            ++n;
        }
        goto L100;

      L5:
        for (a = 1; a <= nxyza; ++a)
        {
            xcol1 = row[n * nrow + 1];
            xcol2 = row[n * nrow + 2];
            xcol3 = row[n * nrow + 3];
            xcol4 = row[n * nrow + 4];
            xcol5 = row[n * nrow + 5];
            rot1 = rot[b * nrow + 1];
            rot2 = rot[b * nrow + 2];
            rot3 = rot[b * nrow + 3];
            rot4 = rot[b * nrow + 4];
            for (j = 1; j <= m; ++j)
            {
                y[j + n * m] = rot1 * x[j + xcol1 * m] +
                    rot2 * x[j + xcol2 * m] +
                    rot3 * x[j + xcol3 * m] +
                    rot4 * x[j + xcol4 * m] +
                    x[j + xcol5 * m];
            }
            ++n;
        }
        goto L100;

      L6:
        for (a = 1; a <= nxyza; ++a)
        {
            xcol1 = row[n * nrow + 1];
            xcol2 = row[n * nrow + 2];
            xcol3 = row[n * nrow + 3];
            xcol4 = row[n * nrow + 4];
            xcol5 = row[n * nrow + 5];
            xcol6 = row[n * nrow + 6];
            rot1 = rot[b * nrow + 1];
            rot2 = rot[b * nrow + 2];
            rot3 = rot[b * nrow + 3];
            rot4 = rot[b * nrow + 4];
            rot5 = rot[b * nrow + 5];
            for (j = 1; j <= m; ++j)
            {
                y[j + n * m] = rot1 * x[j + xcol1 * m] +
                    rot2 * x[j + xcol2 * m] +
                    rot3 * x[j + xcol3 * m] +
                    rot4 * x[j + xcol4 * m] +
                    rot5 * x[j + xcol5 * m] +
                    x[j + xcol6 * m];
            }
            ++n;
        }
        goto L100;

      L7:
        for (a = 1; a <= nxyza; ++a)
        {
            xcol1 = row[n * nrow + 1];
            xcol2 = row[n * nrow + 2];
            xcol3 = row[n * nrow + 3];
            xcol4 = row[n * nrow + 4];
            xcol5 = row[n * nrow + 5];
            xcol6 = row[n * nrow + 6];
            xcol7 = row[n * nrow + 7];
            rot1 = rot[b * nrow + 1];
            rot2 = rot[b * nrow + 2];
            rot3 = rot[b * nrow + 3];
            rot4 = rot[b * nrow + 4];
            rot5 = rot[b * nrow + 5];
            rot6 = rot[b * nrow + 6];
            for (j = 1; j <= m; ++j)
            {
                y[j + n * m] = rot1 * x[j + xcol1 * m] +
                    rot2 * x[j + xcol2 * m] +
                    rot3 * x[j + xcol3 * m] +
                    rot4 * x[j + xcol4 * m] +
                    rot5 * x[j + xcol5 * m] +
                    rot6 * x[j + xcol6 * m] +
                    x[j + xcol7 * m];
            }
            ++n;
        }
        goto L100;

      L8:
        for (a = 1; a <= nxyza; ++a)
        {
            xcol1 = row[n * nrow + 1];
            xcol2 = row[n * nrow + 2];
            xcol3 = row[n * nrow + 3];
            xcol4 = row[n * nrow + 4];
            xcol5 = row[n * nrow + 5];
            xcol6 = row[n * nrow + 6];
            xcol7 = row[n * nrow + 7];
            xcol8 = row[n * nrow + 8];
            rot1 = rot[b * nrow + 1];
            rot2 = rot[b * nrow + 2];
            rot3 = rot[b * nrow + 3];
            rot4 = rot[b * nrow + 4];
            rot5 = rot[b * nrow + 5];
            rot6 = rot[b * nrow + 6];
            rot7 = rot[b * nrow + 7];
            for (j = 1; j <= m; ++j)
            {
                y[j + n * m] = rot1 * x[j + xcol1 * m] + rot2 * x[j
                                                                            +
                                                                            xcol2
                                                                            *
                                                                            m]
                    + rot3 * x[j + xcol3 * m] + rot4 * x[j +
                                                              xcol4 *
                                                              m] +
                    rot5 * x[j + xcol5 * m] + rot6 * x[j +
                                                            xcol6 * m] +
                    rot7 * x[j + xcol7 * m] + x[j + xcol8 * m];
            }
            ++n;
        }
        goto L100;


/*             ...# of rows > 8. Perform transformations in bundles */
/*                of 8 rows. */


      L9:
        for (a = 1; a <= nxyza; ++a)
        {
            xcol1 = row[mrow - 7 + n * nrow];
            xcol2 = row[mrow - 6 + n * nrow];
            xcol3 = row[mrow - 5 + n * nrow];
            xcol4 = row[mrow - 4 + n * nrow];
            xcol5 = row[mrow - 3 + n * nrow];
            xcol6 = row[mrow - 2 + n * nrow];
            xcol7 = row[mrow - 1 + n * nrow];
            xcol8 = row[mrow + n * nrow];
            rot1 = rot[mrow - 7 + b * nrow];
            rot2 = rot[mrow - 6 + b * nrow];
            rot3 = rot[mrow - 5 + b * nrow];
            rot4 = rot[mrow - 4 + b * nrow];
            rot5 = rot[mrow - 3 + b * nrow];
            rot6 = rot[mrow - 2 + b * nrow];
            rot7 = rot[mrow - 1 + b * nrow];
            for (j = 1; j <= m; ++j)
            {
                y[j + n * m] = rot1 * x[j + xcol1 * m] + rot2 * x[j
                                                                            +
                                                                            xcol2
                                                                            *
                                                                            m]
                    + rot3 * x[j + xcol3 * m] + rot4 * x[j +
                                                              xcol4 *
                                                              m] +
                    rot5 * x[j + xcol5 * m] + rot6 * x[j +
                                                            xcol6 * m] +
                    rot7 * x[j + xcol7 * m] + x[j + xcol8 * m];
            }
            nleft = mrow - 8;
            nstep = nleft / 8;
            for (i= 1; i <= nstep; ++i)
            {
                xcol1 = row[nleft - 7 + n * nrow];
                xcol2 = row[nleft - 6 + n * nrow];
                xcol3 = row[nleft - 5 + n * nrow];
                xcol4 = row[nleft - 4 + n * nrow];
                xcol5 = row[nleft - 3 + n * nrow];
                xcol6 = row[nleft - 2 + n * nrow];
                xcol7 = row[nleft - 1 + n * nrow];
                xcol8 = row[nleft + n * nrow];
                rot1 = rot[nleft - 7 + b * nrow];
                rot2 = rot[nleft - 6 + b * nrow];
                rot3 = rot[nleft - 5 + b * nrow];
                rot4 = rot[nleft - 4 + b * nrow];
                rot5 = rot[nleft - 3 + b * nrow];
                rot6 = rot[nleft - 2 + b * nrow];
                rot7 = rot[nleft - 1 + b * nrow];
                rot8 = rot[nleft + b * nrow];
                for (j = 1; j <= m; ++j)
                {
                    y[j + n * m] = y[j + n * m] + rot1 * x[j +
                                                                     xcol1 *
                                                                     m] +
                        rot2 * x[j + xcol2 * m] + rot3 * x[j +
                                                                xcol3 *
                                                                m] +
                        rot4 * x[j + xcol4 * m] + rot5 * x[j +
                                                                xcol5 *
                                                                m] +
                        rot6 * x[j + xcol6 * m] + rot7 * x[j +
                                                                xcol7 *
                                                                m] +
                        rot8 * x[j + xcol8 * m];
                }
                nleft += -8;
            }
            for (i = nleft; i >= 1; --i)
            {
                xcol1 = row[i + n * nrow];
                rot1 = rot[i + b * nrow];
                for (j = 1; j <= m; ++j)
                {
                    y[j + n * m] += rot1 * x[j + xcol1 * m];
                }
            }
            ++n;
        }
      L100:
        ;
    }

    return 0;
}
