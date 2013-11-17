#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__PREPARE_CTR */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine prepares the for contracting the */
/*                primitive batches. Everything that needs still to be */
/*                done at the stage when calling this routine should be */
/*                placed here. At the moment, we have the following: */
/*                       i) copy the exponential prefactors */
/*                      ii) generate all the A,B,C,D norms */
/*                     iii) add the contribution due to the */
/*                          overall integral prefactor and the */
/*                          s-/p-shell type norm into one of */
/*                          the A,B,C,D norms, depending on */
/*                          which has the least elements */
/*                      iv) initialize the contraction batch */
/*                  Input: */
/*                    NCSIZE       =  size of the contraction batch */
/*                    NIJ(KL)      =  total # of ij(kl) primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair A,B(C,D) */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for contraction shells x = A,B,C,D */
/*                    SHELLx       =  the shell types for contraction */
/*                                    shells x = A,B,C,D */
/*                    ALPHAx       =  the primitive exponents for */
/*                                    contraction shells x = A,B,C,D */
/*                    PREFACT      =  overall prefactor for all integrals */
/*                    SPNORM       =  normalization factor due to */
/*                                    presence of s- and p-type shells. */
/*                                    For each s-type shell there is a */
/*                                    factor of 1 and for each p-type */
/*                                    shell a factor of 2 */
/*                    EQUALxy      =  indicates, if csh x and csh y are */
/*                                    considered to be equal for the */
/*                                    pairs xy = AB and CD */
/*                    BLOCKED      =  if false, there will be no need */
/*                                    to block the contraction step over */
/*                                    the set of primitives and thus as */
/*                                    a consequence there is no need to */
/*                                    initialize the contraction batch */
/*                    RHO          =  NIJ exponential prefactors rho(a,b) */
/*                                    + NKL exponential prefactors */
/*                                    rho(c,d), in that order */
/*                  Output: */
/*                    NORMx        =  the normalization factors due to */
/*                                    the primitive exponents for the */
/*                                    contraction shells x = A,B,C,D */
/*                    RHOAB(CD)    =  the complete set of NIJ (NKL) */
/*                                    exponential prefactors between */
/*                                    contraction shells A and B */
/*                                    (C and D) */
/*                    CBATCH       =  contraction batch initialized */
/*                                    to zero (if needed) */
/* ------------------------------------------------------------------------ */
int erd__prepare_ctr (int ncsize, int nij, int nkl,
                  int npgtoa, int npgtob,
                  int npgtoc, int npgtod,
                  int shella, int shellb,
                  int shellc, int shelld,
                  double *alphaa, double *alphab,
                  double *alphac, double *alphad,
                  double prefact, double spnorm,
                  int equalab, int equalcd,
                  int blocked, double *rho,
                  double *norma, double *normb,
                  double *normc, double *normd,
                  double *rhoab, double *rhocd, double *cbatch)
{
    int n;
    int npmin;
    double factor;
    double power;

/*             ...copy the exponential prefactors. This has to be */
/*                done with care, as the memory location for the */
/*                final RHOAB and RHOCD arrays might overlap with */
/*                the input RHO array. The following diagram shows */
/*                how such overlap might happen: */
/*                       RHO  ->  |  NIJ  |  NKL  | */
/*                                 \       \       \ */
/*                                  \       \       \ */
/*                                   | RHOAB | RHOCD | */
/*                We are always safe, if we start copying from the */
/*                last element of RHO downwards. */
    for (n = nkl - 1; n >= 0; n--)
    {
        rhocd[n] = rho[nij + n];
    }
    for (n = nij - 1; n >= 0; n--)
    {
        rhoab[n] = rho[n];
    }

    // normalize a and b
    power = (double)shella * 0.5 + 0.75;
    for (n = 0; n < npgtoa; n++)
    {
        if (alphaa[n] == 0.0)
        {
            norma[n] = 1.0;
        }
        else
        {
            norma[n] = pow (alphaa[n], power);
        }
    }
    if (equalab == 1)
    {
        memcpy (normb, norma, sizeof(double) * npgtoa);
    }
    else
    {
        power = (double)shellb * 0.5 + 0.75;
        for (n = 0; n < npgtob; n++)
        {
            if (alphab[n] == 0.0)
            {
                normb[n] = 1.0;
            }
            else
            {
                normb[n] = pow (alphab[n], power);
            }
        }
    }

    // c and d
    power = (double)shellc * 0.5 + 0.75;
    for (n = 0; n < npgtoc; n++)
    {
        if (alphac[n] == 0.0)
        {
            normc[n] = 1.0;
        }
        else
        {
            normc[n] = pow (alphac[n], power);
        }
    }
    if (equalcd == 1)
    {
        memcpy (normd, normc, sizeof(double) * npgtoc);
    }
    else
    {
        power = (double)shelld * 0.5 + 0.75;
        for (n = 0; n < npgtod; n++)
        {
            if (alphad[n] == 0.0)
            {
                normd[n] = 1.0;
            }
            else
            {
                normd[n] = pow (alphad[n], power);
            }
        }
    }
    


/*             ...rescale one of the A,B,C,D norms, which has the */
/*                least number of elements. */
    factor = prefact * spnorm;
    npmin = npgtoa < npgtob ? npgtoa : npgtob;
    npmin = npmin < npgtoc ? npmin : npgtoc;
    npmin = npmin < npgtod ? npmin : npgtod;
    if (npgtoa == npmin)
    {
        for (n = 0; n < npgtoa; n++)
        {
            norma[n] = factor * norma[n];
        }
    }
    else if (npgtob == npmin)
    {
        for (n = 0; n < npgtob; n++)
        {
            normb[n] = factor * normb[n];
        }
    }
    else if (npgtoc == npmin)
    {
        for (n = 0; n < npgtoc; n++)
        {
            normc[n] = factor * normc[n];
        }
    }
    else
    {
        for (n = 0; n < npgtod; n++)
        {
            normd[n] = factor * normd[n];
        }
    }


/*             ...initialize contraction batch (if necessary). */
    if (blocked)
    {
        memset (cbatch, 0, sizeof(double) * ncsize);
    }

    return 0;
}
