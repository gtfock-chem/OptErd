#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__XYZ_TO_RY_ABCD */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__XYZ_TO_RY_MATRIX */
/*  DESCRIPTION : This operation generates all the info needed for */
/*                cartesian -> spherical transformation of a batch */
/*                of cartesian integrals (AB|CD). No duplicate info */
/*                is generated, i.e. the above transmitted pointers */
/*                of the type Z00x, I0x1 and I0x2 with x=A,B,C,D may */
/*                coincide if the shell types among the A,B,C,D are */
/*                equal (see below). The transformation data is placed */
/*                in two arrays (int and flp) at the appropriate */
/*                locations. */
/*                Input: */
/*                    NXYZx = dimension of xyz-monomial basis */
/*                            corresponding to the shell SHELLx */
/*                            for x=A,B,C,D. */
/*                     NRYx = dimension of ry-spherical basis */
/*                            corresponding to the shell SHELLx */
/*                            for x=A,B,C,D. */
/*                   SHELLx = the shell types for A,B,C,D. */
/*                I(Z)START = Starting location for the int (flp) */
/*                            data. */
/*                Output: */
/*                    NROWx = maximum # of xyz-monomials contributing */
/*                            to the ry-components for x=A,B,C,D. */
/*                    NROTx = maximum # of elements in transformation */
/*                            matrix for x=A,B,C,D. This is equal to */
/*                            NROWx times NRYx. */
/*                     Z00x = pointer for the transformation matrix */
/*                            elements for x=A,B,C,D. */
/*                     I0x1 = pointer for the # of row indices leading */
/*                            to non-zero elements in the transformation */
/*                            matrix for x=A,B,C,D. */
/*                     I0x2 = pointer for the row index labels of the */
/*                            non-zero elements in the transformation */
/*                            matrix for x=A,B,C,D. */
/*                 I(Z)USED = # of int (flp) words used. */
/*                 I(Z)CORE = The int (flp) arrays holding the */
/*                            transformation data. */
/*                If any of the A,B,C,D shell labels are equal, their */
/*                offsets will be set equal in the following sequence, */
/*                governed by their order of usage: */
/*                    1) If shell C = D, then: */
/*                              Z00C = Z00D */
/*                              I0C1 = I0D1 */
/*                              I0C2 = I0D2 */
/*                    2) If shell B = C or D, then: */
/*                              Z00B = Z00C or Z00D */
/*                              I0B1 = I0C1 or I0D1 */
/*                              I0B2 = I0C2 or I0D2 */
/*                    3) If shell A = B or C or D, then: */
/*                              Z00A = Z00B or Z00C or Z00D */
/*                              I0A1 = I0B1 or I0C1 or I0D1 */
/*                              I0A2 = I0B2 or I0C2 or I0D2 */
/*                Only mutually different transformation matrices + */
/*                associated data are generated. */
/* ------------------------------------------------------------------------ */
int erd__xyz_to_ry_abcd (int nxyza, int nxyzb, int nxyzc, int nxyzd,
                         int nrya, int nryb, int nryc, int nryd,
                         int shella, int shellb,
                         int shellc, int shelld,
                         int istart, int zstart,
                         int *nrowa, int *nrowb,
                         int *nrowc, int *nrowd,
                         int *nrota, int *nrotb,
                         int *nrotc, int *nrotd,
                         int *z00a, int *z00b, int *z00c, int *z00d,
                         int *i0a1, int *i0b1, int *i0c1, int *i0d1,
                         int *i0a2, int *i0b2, int *i0c2, int *i0d2,
                         int *iused, int *zused,
                         int *icore, double *zcore)
{
    int z0at, z0bt, z0ct, z0dt;
    int bdata, cdata;

    --zcore;
    --icore;

/*             ...shell D data. */
    *iused = 0;
    *zused = 0;
    if (shelld > 1)
    {
        *nrowd = (shelld / 2 + 1) * (shelld / 2 + 2) / 2;
        *nrotd = *nrowd * nryd;
        *z00d = zstart;
        z0dt = *z00d + *nrotd;
        *i0d1 = istart;
        *i0d2 = *i0d1 + nryd;
        erd__xyz_to_ry_matrix (nxyzd, *nrowd, shelld,
                               &zcore[z0dt], &icore[*i0d1],
                               &icore[*i0d2], &zcore[*z00d]);
        *iused = nryd + *nrotd;
        *zused = *nrotd;
    }
    else
    {
        *nrowd = 0;
        *nrotd = 0;
        *z00d = zstart;
        *i0d2 = istart;
    }

/*             ...shell C data. */
    cdata = 0;
    if (shellc > 1)
    {
        if (shellc == shelld)
        {
            *z00c = *z00d;
            *i0c1 = *i0d1;
            *i0c2 = *i0d2;
            *nrowc = *nrowd;
            *nrotc = *nrotd;
        }
        else
        {
            *nrowc = (shellc / 2 + 1) * (shellc / 2 + 2) / 2;
            *nrotc = *nrowc * nryc;
            *z00c = *z00d + *nrotd;
            z0ct = *z00c + *nrotc;
            *i0c1 = *i0d2 + *nrotd;
            *i0c2 = *i0c1 + nryc;
            erd__xyz_to_ry_matrix (nxyzc, *nrowc, shellc,
                                   &zcore[z0ct], &icore[*i0c1],
                                   &icore[*i0c2], &zcore[*z00c]);
            *iused = *iused + nryc + *nrotc;
            *zused += *nrotc;
            cdata = 1;
        }
    }
    else
    {
        *nrowc = 0;
        *nrotc = 0;
        *z00c = *z00d;
        *i0c2 = *i0d2;
    }

/*             ...shell B data (being careful, using SHELLD data if */
/*                SHELLC data is not present!). */
    bdata = 0;
    if (shellb > 1)
    {
        if (shellb == shellc)
        {
            *z00b = *z00c;
            *i0b1 = *i0c1;
            *i0b2 = *i0c2;
            *nrowb = *nrowc;
            *nrotb = *nrotc;
        }
        else if (shellb == shelld)
        {
            *z00b = *z00d;
            *i0b1 = *i0d1;
            *i0b2 = *i0d2;
            *nrowb = *nrowd;
            *nrotb = *nrotd;
        }
        else
        {
            *nrowb = (shellb / 2 + 1) * (shellb / 2 + 2) / 2;
            *nrotb = *nrowb * nryb;
            if (cdata)
            {
                *z00b = *z00c + *nrotc;
                z0bt = *z00b + *nrotb;
                *i0b1 = *i0c2 + *nrotc;
                *i0b2 = *i0b1 + nryb;
            }
            else
            {
                *z00b = *z00d + *nrotd;
                z0bt = *z00b + *nrotb;
                *i0b1 = *i0d2 + *nrotd;
                *i0b2 = *i0b1 + nryb;
            }
            erd__xyz_to_ry_matrix (nxyzb, *nrowb, shellb,
                                   &zcore[z0bt], &icore[*i0b1],
                                   &icore[*i0b2], &zcore[*z00b]);
            *iused = *iused + nryb + *nrotb;
            *zused += *nrotb;
            bdata = 1;
        }
    }
    else
    {
        *nrowb = 0;
        *nrotb = 0;
        *z00b = *z00c;
        *i0b2 = *i0c2;
    }


/*             ...shell A data (being careful, using SHELLC data if */
/*                SHELLB data is not present or using SHELLD data if */
/*                also SHELLC data is not present!). */
    if (shella > 1)
    {
        if (shella == shellb)
        {
            *z00a = *z00b;
            *i0a1 = *i0b1;
            *i0a2 = *i0b2;
            *nrowa = *nrowb;
            *nrota = *nrotb;
        }
        else if (shella == shellc)
        {
            *z00a = *z00c;
            *i0a1 = *i0c1;
            *i0a2 = *i0c2;
            *nrowa = *nrowc;
            *nrota = *nrotc;
        }
        else if (shella == shelld)
        {
            *z00a = *z00d;
            *i0a1 = *i0d1;
            *i0a2 = *i0d2;
            *nrowa = *nrowd;
            *nrota = *nrotd;
        }
        else
        {
            *nrowa = (shella / 2 + 1) * (shella / 2 + 2) / 2;
            *nrota = *nrowa * nrya;
            if (bdata)
            {
                *z00a = *z00b + *nrotb;
                z0at = *z00a + *nrota;
                *i0a1 = *i0b2 + *nrotb;
                *i0a2 = *i0a1 + nrya;
            }
            else if (cdata)
            {
                *z00a = *z00c + *nrotc;
                z0at = *z00a + *nrota;
                *i0a1 = *i0c2 + *nrotc;
                *i0a2 = *i0a1 + nrya;
            }
            else
            {
                *z00a = *z00d + *nrotd;
                z0at = *z00a + *nrota;
                *i0a1 = *i0d2 + *nrotd;
                *i0a2 = *i0a1 + nrya;
            }
            erd__xyz_to_ry_matrix (nxyza, *nrowa, shella,
                                   &zcore[z0at], &icore[*i0a1],
                                   &icore[*i0a2], &zcore[*z00a]);
            *iused = *iused + nrya + *nrota;
            *zused += *nrota;
        }
    }
     
    return 0;
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
