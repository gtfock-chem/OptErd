#include <stdio.h>
#include <stdlib.h>


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__TRANSPOSE_BATCH */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine transposes a batch of NROW x NCOL */
/*                integrals to an output batch of NCOL x NROW integrals. */
/*                The routine uses tiling for better use of cache. */
/*                  Input: */
/*                    NROW,NCOL    =  # of rows and columns in input */
/*                                    integral batch */
/*                    TILE         =  tile size */
/*                    BATCH        =  input integral batch */
/*                  Output: */
/*                    OBATCH       =  output integral batch */
/* ------------------------------------------------------------------------ */
int erd__transpose_batch (int nrow, int ncol,
                          double *batch, double *obatch)
{
    int i;
    int j;
    
    for (j = 0; j < nrow; j++)
    {
        for (i = 0; i < ncol; i++)
        {
            obatch[i + j * ncol] = batch[j + i * nrow];
        }
    }

    return 0;
}