#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "erd.h"


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__HRR_MATRIX */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__HRR_STEP */
/*  DESCRIPTION : This operation constructs the whole HRR transformation */
/*                matrix T, which will contain the information for */
/*                the complete transformation in xyz-basis: */
/*                             (e0| --> (ab|  ;  e = a+b ; a>=b */
/*                The matrix is constructed stepwise, starting from */
/*                a unit matrix and operating on its columns with */
/*                elementary HRR steps. */
/*                  Input: */
/*                    NROTHRR   = maximum # of elements of T and ROW */
/*                                matrix expected during construction */
/*                                of final T and ROW */
/*                    NCOLHRR   = maximum # of columns of T and ROW */
/*                                matrix expected during construction */
/*                                of final T and ROW */
/*                    NXYZET    = monomial dimension of the (e0| part. */
/*                    NXYZA     = monomial dimension for shell a */
/*                    NXYZP     = monomial dimension for shell p=a+b */
/*                    SHELLx    = shell types for x=a,b,p */
/*                    NABCOOR   = # of nonzero coordinate differences */
/*                                between nuclear centers A and B */
/*                    ABx       = the coordinate differences between */
/*                                nuclear centers A and B for x=X,Y,Z */
/*                    WORK      = scratch space used for assembly of */
/*                                final NROW vector and T and ROW */
/*                                matrices */
/*                  Output: */
/*                    IN1,IN2   = starting index position of final */
/*                                NROW vector (IN1) and final T and */
/*                                ROW matrices (IN2) in the big NROW */
/*                                T and ROW arrays (which are 2x the */
/*                                maximum size) */
/*                    NROWOUT   = maximum # of nonzero row labels */
/*                                of matrix T and ROW */
/*                    NROW      = vector containing # of nonzero */
/*                                entries in columns of T and ROW matrix */
/*                    ROW       = the nonzero row labels of matrix T */
/*                    T         = the HRR transformation matrix */
/*                The xyz-basis for the a- and b-parts in columns */
/*                of matrix T will be ordered such that a preceeds b. */
/* ------------------------------------------------------------------------ */
__attribute__((target(mic))) int erd__hrr_matrix (
		     int nrothrr, int ncolhrr,
                     int nxyzet, int nxyza, int nxyzp,
                     int shella, int shellb, int shellp,
                     int nabcoor, double abx, double aby, double abz,
                     int *work, int *in1, int *in2,
                     int *nrowout, int *nrow,
                     int *row, double *t)
{
    int add[3] = {0, 0, 1};

    int i, j, l, m, n, ngh, out1, out2, ngho, base1, base2,
        tleap;
    int nxyzg, nxyzh, nxyzi, shellg, shellh, nrowin, nxyzgo,
        nxyzho;

/* ------------------------------------------------------------------------ */
/*             ...accumulate T. */
    *in1 = 0;
    *in2 = 0;
    out1 = *in1 + ncolhrr;
    out2 = *in2 + nrothrr;

/*             ...form initial 'unit' T. */
    for (i = 0; i < nxyzet; ++i)
    {
        nrow[i] = 1;
        row[i] = i + 1;
        t[i] = 1.0;
    }
/*             ...build up the HRR transformation matrix + data. */
    ngh = nxyzet;
    nxyzg = nxyzet;
    nxyzh = 1;
    nxyzgo = nxyzet;
    nxyzho = 1;
    nxyzi = nxyzp;
    shellg = shellp;
    nrowin = 1;
    *nrowout = 1;
    for (shellh = 1; shellh <= shellb; ++shellh)
    {
        nxyzgo -= nxyzi;
        nxyzho = nxyzho + shellh + 1;
        ngho = nxyzgo * nxyzho;
        if (nabcoor == 3)
        {
            m = shellh / 3 + 1;
            *nrowout += m * (m + add[shellh % 3]);
        }
        else if (nabcoor == 2)
        {
            *nrowout = *nrowout + shellh / 2 + 1;
        }
        else if (nabcoor == 1)
        {
            ++(*nrowout);
        }
        erd__hrr_step (ngh, ngho, nrowin, *nrowout,
                       nxyza, nxyzi, nxyzg, nxyzh, nxyzgo,
                       shella, shellg, shellh - 1,
                       abx, aby, abz,
                       work, &nrow[*in1], &row[*in2], &t[*in2],
                       &nrow[out1], &row[out2], &t[out2]);
        nxyzh = nxyzho;
        if (shellh != shellb)
        {
            ngh = ngho;
            nxyzg = nxyzgo;
            nxyzi = nxyzi - shellg - 1;
            --shellg;
            nrowin = *nrowout;
        }
        i = *in1;
        *in1 = out1;
        out1 = i;
        i = *in2;
        *in2 = out2;
        out2 = i;
    }

/*             ...the resulting T matrix of dimension NROWOUT x (NXYZA */
/*                x NXYZH) has the following structure: the columns */
/*                are such that they come in NXYZA identical copies */
/*                NXYZH times. Hence we can condense the T matrix into */
/*                only NROWOUT x NXYZH distinct elements. The same */
/*                applies to the array NROW of size NXYZA x NXYZH, which */
/*                also allows for condensation into only NXYZH elements. */
    l = 1;
    m = 0;
    n = 0;
    base1 = *in1;
    base2 = *in2;
    tleap = *nrowout * nxyza;
    for (j = 0; j < nxyzh; ++j)
    {
        nrow[base1 + j] = nrow[base1 + l];
        for (i = 0; i < *nrowout; ++i)
        {
            t[base2 + m + i] = t[base2 + n + i];
        }
        l += nxyza;
        m += *nrowout;
        n += tleap;
    }

    *in1 = *in1 + 1;
    *in2 = *in2 + 1;
    return 0;
}
