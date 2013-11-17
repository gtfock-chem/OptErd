#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


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
int erd__spherical_transform (int m, int nrow, int nxyz, int nry,
                              int *lrow, int *row, double *rot,
                              double *x, double *y)
{
    int row_offset, x_offset, y_offset, rot_offset;
    
    /* Local variables */
    int i, n, r;
    double rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8;
    int mrow, xcol1, xcol2, xcol3, xcol4, xcol5, xcol6, xcol7,
        xcol8, nbase, nleft, nstep, nrest;

/*             ...perform the cartesian -> spherical transformation. */
/*                Use basic row grouping of the transformation */
/*                to improve cache line reusing. */
    x_offset = 1 + m * 1;
    x -= x_offset;
    y_offset = 1 + m * 1;
    y -= y_offset;
    rot_offset = 1 + nrow * 1;
    rot -= rot_offset;
    row_offset = 1 + nrow * 1;
    row -= row_offset;
    --lrow;

    for (r = 1; r <= nry; ++r)
    {
        mrow = lrow[r];
        switch (MIN (mrow, 9))
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
        xcol1 = row[r * nrow + 1];
        rot1 = rot[r * nrow + 1];
        for (n = 1; n <= m; ++n)
        {
            y[n + r * m] = rot1 * x[n + xcol1 * m];
        }
        goto L100;

      L2:
        xcol1 = row[r * nrow + 1];
        xcol2 = row[r * nrow + 2];
        rot1 = rot[r * nrow + 1];
        rot2 = rot[r * nrow + 2];
        for (n = 1; n <= m; ++n)
        {
            y[n + r * m] = rot1 * x[n + xcol1 * m] + rot2 * x[n +
                                                                          xcol2
                                                                          *
                                                                          m];
        }
        goto L100;

      L3:
        xcol1 = row[r * nrow + 1];
        xcol2 = row[r * nrow + 2];
        xcol3 = row[r * nrow + 3];
        rot1 = rot[r * nrow + 1];
        rot2 = rot[r * nrow + 2];
        rot3 = rot[r * nrow + 3];
        for (n = 1; n <= m; ++n)
        {
            y[n + r * m] = rot1 * x[n + xcol1 * m] + rot2 * x[n +
                                                                          xcol2
                                                                          *
                                                                          m]
                + rot3 * x[n + xcol3 * m];
        }
        goto L100;

      L4:
        xcol1 = row[r * nrow + 1];
        xcol2 = row[r * nrow + 2];
        xcol3 = row[r * nrow + 3];
        xcol4 = row[r * nrow + 4];
        rot1 = rot[r * nrow + 1];
        rot2 = rot[r * nrow + 2];
        rot3 = rot[r * nrow + 3];
        rot4 = rot[r * nrow + 4];
        for (n = 1; n <= m; ++n)
        {
            y[n + r * m] = rot1 * x[n + xcol1 * m] + rot2 * x[n +
                                                                          xcol2
                                                                          *
                                                                          m]
                + rot3 * x[n + xcol3 * m] + rot4 * x[n + xcol4 * m];
        }
        goto L100;

      L5:
        xcol1 = row[r * nrow + 1];
        xcol2 = row[r * nrow + 2];
        xcol3 = row[r * nrow + 3];
        xcol4 = row[r * nrow + 4];
        xcol5 = row[r * nrow + 5];
        rot1 = rot[r * nrow + 1];
        rot2 = rot[r * nrow + 2];
        rot3 = rot[r * nrow + 3];
        rot4 = rot[r * nrow + 4];
        rot5 = rot[r * nrow + 5];
        for (n = 1; n <= m; ++n)
        {
            y[n + r * m] = rot1 * x[n + xcol1 * m] + rot2 * x[n +
                                                                          xcol2
                                                                          *
                                                                          m]
                + rot3 * x[n + xcol3 * m] + rot4 * x[n +
                                                          xcol4 * m] +
                rot5 * x[n + xcol5 * m];
        }
        goto L100;

      L6:
        xcol1 = row[r * nrow + 1];
        xcol2 = row[r * nrow + 2];
        xcol3 = row[r * nrow + 3];
        xcol4 = row[r * nrow + 4];
        xcol5 = row[r * nrow + 5];
        xcol6 = row[r * nrow + 6];
        rot1 = rot[r * nrow + 1];
        rot2 = rot[r * nrow + 2];
        rot3 = rot[r * nrow + 3];
        rot4 = rot[r * nrow + 4];
        rot5 = rot[r * nrow + 5];
        rot6 = rot[r * nrow + 6];
        for (n = 1; n <= m; ++n)
        {
            y[n + r * m] = rot1 * x[n + xcol1 * m] + rot2 * x[n +
                                                                          xcol2
                                                                          *
                                                                          m]
                + rot3 * x[n + xcol3 * m] + rot4 * x[n +
                                                          xcol4 * m] +
                rot5 * x[n + xcol5 * m] + rot6 * x[n + xcol6 * m];
        }
        goto L100;

      L7:
        xcol1 = row[r * nrow + 1];
        xcol2 = row[r * nrow + 2];
        xcol3 = row[r * nrow + 3];
        xcol4 = row[r * nrow + 4];
        xcol5 = row[r * nrow + 5];
        xcol6 = row[r * nrow + 6];
        xcol7 = row[r * nrow + 7];
        rot1 = rot[r * nrow + 1];
        rot2 = rot[r * nrow + 2];
        rot3 = rot[r * nrow + 3];
        rot4 = rot[r * nrow + 4];
        rot5 = rot[r * nrow + 5];
        rot6 = rot[r * nrow + 6];
        rot7 = rot[r * nrow + 7];
        for (n = 1; n <= m; ++n)
        {
            y[n + r * m] = rot1 * x[n + xcol1 * m] + rot2 * x[n +
                                                                          xcol2
                                                                          *
                                                                          m]
                + rot3 * x[n + xcol3 * m] + rot4 * x[n +
                                                          xcol4 * m] +
                rot5 * x[n + xcol5 * m] + rot6 * x[n + xcol6 * m] +
                rot7 * x[n + xcol7 * m];
        }
        goto L100;

      L8:
        xcol1 = row[r * nrow + 1];
        xcol2 = row[r * nrow + 2];
        xcol3 = row[r * nrow + 3];
        xcol4 = row[r * nrow + 4];
        xcol5 = row[r * nrow + 5];
        xcol6 = row[r * nrow + 6];
        xcol7 = row[r * nrow + 7];
        xcol8 = row[r * nrow + 8];
        rot1 = rot[r * nrow + 1];
        rot2 = rot[r * nrow + 2];
        rot3 = rot[r * nrow + 3];
        rot4 = rot[r * nrow + 4];
        rot5 = rot[r * nrow + 5];
        rot6 = rot[r * nrow + 6];
        rot7 = rot[r * nrow + 7];
        rot8 = rot[r * nrow + 8];
        for (n = 1; n <= m; ++n)
        {
            y[n + r * m] = rot1 * x[n + xcol1 * m] + rot2 * x[n +
                                                                          xcol2
                                                                          *
                                                                          m]
                + rot3 * x[n + xcol3 * m] + rot4 * x[n +
                                                          xcol4 * m] +
                rot5 * x[n + xcol5 * m] + rot6 * x[n + xcol6 * m] +
                rot7 * x[n + xcol7 * m] + rot8 * x[n + xcol8 * m];
        }
        goto L100;


/*             ...# of rows > 8. Perform transformations in bundles */
/*                of 8 rows. */
      L9:
        xcol1 = row[r * nrow + 1];
        xcol2 = row[r * nrow + 2];
        xcol3 = row[r * nrow + 3];
        xcol4 = row[r * nrow + 4];
        xcol5 = row[r * nrow + 5];
        xcol6 = row[r * nrow + 6];
        xcol7 = row[r * nrow + 7];
        xcol8 = row[r * nrow + 8];
        rot1 = rot[r * nrow + 1];
        rot2 = rot[r * nrow + 2];
        rot3 = rot[r * nrow + 3];
        rot4 = rot[r * nrow + 4];
        rot5 = rot[r * nrow + 5];
        rot6 = rot[r * nrow + 6];
        rot7 = rot[r * nrow + 7];
        rot8 = rot[r * nrow + 8];
        for (n = 1; n <= m; ++n)
        {
            y[n + r * m] = rot1 * x[n + xcol1 * m] + rot2 * x[n +
                                                                          xcol2
                                                                          *
                                                                          m]
                + rot3 * x[n + xcol3 * m] + rot4 * x[n +
                                                          xcol4 * m] +
                rot5 * x[n + xcol5 * m] + rot6 * x[n + xcol6 * m] +
                rot7 * x[n + xcol7 * m] + rot8 * x[n + xcol8 * m];
        }
        nbase = 8;
        nleft = mrow - 8;
        nstep = nleft / 8;
        nrest = nleft % 8;
        for (i = 1; i <= nstep; ++i)
        {
            xcol1 = row[nbase + 1 + r * nrow];
            xcol2 = row[nbase + 2 + r * nrow];
            xcol3 = row[nbase + 3 + r * nrow];
            xcol4 = row[nbase + 4 + r * nrow];
            xcol5 = row[nbase + 5 + r * nrow];
            xcol6 = row[nbase + 6 + r * nrow];
            xcol7 = row[nbase + 7 + r * nrow];
            xcol8 = row[nbase + 8 + r * nrow];
            rot1 = rot[nbase + 1 + r * nrow];
            rot2 = rot[nbase + 2 + r * nrow];
            rot3 = rot[nbase + 3 + r * nrow];
            rot4 = rot[nbase + 4 + r * nrow];
            rot5 = rot[nbase + 5 + r * nrow];
            rot6 = rot[nbase + 6 + r * nrow];
            rot7 = rot[nbase + 7 + r * nrow];
            rot8 = rot[nbase + 8 + r * nrow];
            for (n = 1; n <= m; ++n)
            {
                y[n + r * m] = y[n + r * m] + rot1 * x[n +
                                                                     xcol1 *
                                                                     m] +
                    rot2 * x[n + xcol2 * m] + rot3 * x[n +
                                                            xcol3 * m] +
                    rot4 * x[n + xcol4 * m] + rot5 * x[n +
                                                            xcol5 * m] +
                    rot6 * x[n + xcol6 * m] + rot7 * x[n +
                                                            xcol7 * m] +
                    rot8 * x[n + xcol8 * m];
            }
            nbase += 8;
        }
        for (i = 1; i <= nrest; ++i)
        {
            xcol1 = row[nbase + i + r * nrow];
            rot1 = rot[nbase + i + r * nrow];
            for (n = 1; n <= m; ++n)
            {
                y[n + r * m] += rot1 * x[n + xcol1 * m];
            }
        }
      L100:
        ;
    }

    return 0;
}
