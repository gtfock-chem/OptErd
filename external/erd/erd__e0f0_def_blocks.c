#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "erd.h"
#include "erdutil.h"

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__E0F0_DEF_BLOCKS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine handles the memory partitions for the */
/*                generation of the entire contracted (e0|f0) batch. */
/*                It determines the block size partitions (if any) for */
/*                the ij and kl exponent paris of the whole primitive */
/*                [e0|f0] batch and returns pointers to the flp data */
/*                sections needed by the (e0|f0) generation. */
/*                The routine is also equipped with the possibility */
/*                of returning just the minimum / optimum flp memory */
/*                needed, without evaluating the (e0|f0) generation flp */
/*                pointers. This is useful for establishing just the */
/*                overall memory size needed for the (e0|f0) generation. */
/*                The keyword MEMORY has to be set true for this case */
/*                and obviously the value of ZMAX passed is irrelevant. */
/*                  Input: */
/*                    ZMAX         =  maximum flp memory */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for csh x=A,B,C,D */
/*                    SHELLP(Q)    =  the shell sums for csh: A+B (C+D) */
/*                    NIJ(KL)      =  total # of ij(kl) primitive */
/*                                    index pairs corresponding to */
/*                                    the csh pairs A,B(C,D) */
/*                    NRS(TU)      =  total # of rs(tu) contraction */
/*                                    index pairs corresponding to */
/*                                    the csh pairs A,B(C,D) */
/*                    NRSTU        =  total # of rstu contraction */
/*                                    index quadruplets (= NRS*NTU) */
/*                    NGQP         =  # of gaussian quadrature points */
/*                                    (roots) */
/*                    NGQSCR       =  size of gaussian quadrature */
/*                                    scratch space needed to calculate */
/*                                    all the quadrature roots */
/*                    NXYZT        =  total # of cartesian monomial */
/*                                    quadruplets */
/*                    L1CACHE      =  Size of level 1 cache in units of */
/*                                    8 Byte */
/*                    NCTROW       =  minimum # of rows that are */
/*                                    accepted for blocked contractions */
/*                    MEMORY       =  if this keyword is true, the */
/*                                    routine will only determine the */
/*                                    minimum / optimum flp memory and */
/*                                    store these values into NIJBLK and */
/*                                    NKLBLK, respectively (see below) */
/*                  Output: */
/*                    NIJ(KL)BLK   =  if MEMORY is false, they contain */
/*                                    the block sizes for the ij(kl) */
/*                                    primitive indices in order to */
/*                                    perform efficient contractions. */
/*                                    if MEMORY is true, they will have */
/*                                    the values for the minimum / */
/*                                    optimum flp memory */
/*                    NxSIZE       =  size of the primitive integral */
/*                                    block (x=P), contracted integrals */
/*                                    (x=C) and working (x=W) arrays */
/*                    NINT2D       =  space needed for the 2D X,Y,Z */
/*                                    integral arrays */
/*                    MXPRIM       =  the maximum # of primitives */
/*                                    between all i,j,k,l primitives, */
/*                                    i.e. = max (i,j,k,l) */
/*                    MNPRIM       =  the minimum # of primitives */
/*                                    between i and j primitives and */
/*                                    k and l primitives and form the */
/*                                    maximum between these two values, */
/*                                    i.e. = max (min(i,j),min(k,l)) */
/*                    Z.....       =  pointers for space partition of */
/*                                    the flp array (see below) */
/*                The space partitioning of the flp array will be */
/*                as follows: */
/*                 |  Zone 1  |  Zone 2  |  Zone 3  |  Zone 4  | */
/*                   Zone 1 : final (e0|f0) batch */
/*                   Zone 2 : block [e0|f0] and half transformed (e0|f0) */
/*                   Zone 3 : complete set of A,B,C,D norms and */
/*                            complete set of AB and CD exponential */
/*                            prefactors */
/*                   Zone 4 : i) scratch for block [e0|f0] generation */
/*                           ii) scratch for partial (e0|f0) generation */
/*                Memory allocation offsets for the primitive [e0|f0] */
/*                batches generation + contraction to (e0|f0): */
/*                  --- Zone 1 --- */
/*                   ZCBATCH = offset for contracted (e0|f0) batch */
ERD_OFFLOAD void erd__e0f0_def_blocks(int zmax, int npgtoa, int npgtob,
                          int npgtoc, int npgtod,
                          int shellp, int shellq,
                          int nij, int nkl, int ngqp, int ngqscr,
                          int nxyzt, int memory, int *nint2d, int *zcbatch)
{
/*             ...if the MEMORY keyword is activated, the routine */
/*                will only determine the optimum and minimum flp */
/*                memory (in that order), place them respectively */
/*                into the NKLBLK and NIJBLK variables and exit. */
    
    const int zone12 = PAD_LEN(nxyzt);

    const int mijkl = nij * nkl;
    const int mgqijkl = ngqp * mijkl;
    const int pad_mgqijkl = PAD_LEN(mgqijkl);

    if (memory) {
        *nint2d = zone12;
    } else {
        *nint2d = pad_mgqijkl * (shellp + 1) * (shellq + 1);
        assert (zone12 <= zmax);
        
        *zcbatch = 1;
    }
}
