#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__MOVE_RY */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__TRANSPOSE_BATCH */
/*                ERD__MAP_IJKL_TO_IKJL */
/*  DESCRIPTION : This routine moves all ry-components located on the */
/*                far right in array X to a specific position to the left */
/*                in array Y: */
/*                     X (1,2,3,...,RY) ---> Y ( 1, 2, 3,...) */
/*                                   |          |  |  |   | */
/*                                    --->------^--^--^...^ */
/*                The part of X that is not moved (i.e. the # of */
/*                invariant indices to the left in X) has been calculated */
/*                beforehand and is transmitted through argument. */
/*                  Input: */
/*                       NINTGRL    =  total # of integrals */
/*                       NINDEX     =  total # of possible target places */
/*                       NOTMOVE    =  inactive # of indices */
/*                       MOVE       =  # of indices that will be moved */
/*                       NRY        =  # of ry-components to be moved */
/*                       INDEX      =  place 1,2,3,... to which the */
/*                                     ry-components will be placed in */
/*                                     array Y (must be within the range */
/*                                     1 =< INDEX =< NINDEX) */
/*                       TILE       =  the level 1 cache tile */
/*                       X          =  initial set of integrals */
/*                       IXOFF (I)  =  NINDEX-fold array indicating total */
/*                                     # of elements preceeding place */
/*                                     I = 1,2,3,... before the move */
/*                                     (with IXOFF(1)=1) */
/*                  Output: */
/*                       IXOFF (I)  =  updated NINDEX-fold array */
/*                                     indicating total # of elements */
/*                                     preceeding place I = 1,2,3,... */
/*                                     after the move */
/*                       Y          =  final set of integrals */
/* ------------------------------------------------------------------------ */
int erd__move_ry (int nindex, int notmove, int move, int nry,
                  int index, double *x,
                  int *ixoff, double *y)
{
    int i;
/*             ...check, if the move is simply a transposition */
/*                or a more general ijkl -> ikjl move. */
    if (notmove == 1)
    {
        erd__transpose_batch (move, nry, x, y);
    }
    else
    {
        erd__map_ijkl_to_ikjl (notmove, move, nry, 1, x, y);
    }

    for (i = index; i < nindex; i++)
    {
        ixoff[i] *= nry;
    }

    return 0;
}