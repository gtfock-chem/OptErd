#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__1111_DEF_BLOCKS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine handles the memory partitions for the */
/*                generation of the entire contracted (12|34) batches */
/*                for integrals over s- and/or p-functions only. */
/*                It determines the block size partitions (if any) for */
/*                the ij and kl exponent pairs of the whole primitive */
/*                [12|34] batch and returns pointers to the flp data */
/*                sections needed by the (12|34) generation. */
/*                The routine is also equipped with the possibility */
/*                of returning just the minimum / optimum flp memory */
/*                needed, without evaluating the (12|34) generation flp */
/*                pointers. This is useful for establishing just the */
/*                overall memory size needed for the (12|34) generation. */
/*                The keyword MEMORY has to be set true for this case */
/*                and obviously the value of ZMAX passed is irrelevant. */
/*                  Input: */
/*                    ZMAX         =  maximum flp memory */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for csh x=1,2,3,4 */
/*                    NIJ(KL)      =  total # of ij(kl) primitive */
/*                                    index pairs corresponding to */
/*                                    the csh pairs 1,2(3,4) */
/*                    NRS(TU)      =  total # of rs(tu) contraction */
/*                                    index pairs corresponding to */
/*                                    the csh pairs 1,2(3,4) */
/*                    NRSTU        =  total # of rstu contraction */
/*                                    index quadruplets (= NRS*NTU) */
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
/*                   Zone 1: final (12|34) batch */
/*                   Zone 2: block [12|34] and half transformed (12|34) */
/*                   Zone 3 : complete set of 1,2,3,4 norms and */
/*                            complete set (after screening!) of */
/*                            exponential prefactors */
/*                   Zone 4 :  i) scratch for block [12|34] generation */
/*                            ii) scratch for partial (12|34) generation */
/*                Memory allocation offsets for the primitive [12|34] */
/*                batches generation + contraction to (12|34): */
/*                  --- Zone 1 --- */
/*                   ZCBATCH = offset for contracted (12|34) batch */
/*                  --- Zone 2 --- */
/*                   ZPBATCH = offset for (blocked) [12|34] batch (K4) */
/*                  --- Zone 3 --- */
/*                   ZNORM1 = offset for all 1 norms */
/*                   ZNORM2 = offset for all 2 norms */
/*                   ZNORM3 = offset for all 3 norms */
/*                   ZNORM4 = offset for all 4 norms */
/*                   ZRHO12 = offset for all 12 exponential prefactors */
/*                   ZRHO34 = offset for all 34 exponential prefactors */
/*                  --- Zone 4: for block [12|34] only --- */
/*                   ZP = offset for (blocked) P values (K2) */
/*                   ZPX = offset for (blocked) PX values (K2) */
/*                   ZPY = offset for (blocked) PY values (K2) */
/*                   ZPZ = offset for (blocked) PZ values (K2) */
/*                   ZSCPK2 = offset for (blocked) P scale values (K2) */
/*                   ZQ = offset for (blocked) Q values (K2) */
/*                   ZQX = offset for (blocked) QX values (K2) */
/*                   ZQY = offset for (blocked) QY values (K2) */
/*                   ZQZ = offset for (blocked) QZ values (K2) */
/*                   ZSCQK2 = offset for (blocked) Q scale values (K2) */
/*                  --- Zone 4: for contraction only --- */
/*                   ZWORK = offset for contraction working array */
/* ------------------------------------------------------------------------ */
int erd__1111_def_blocks (int zmax, int npgto1, int npgto2,
                          int npgto3, int npgto4,
                          int nij, int nkl, int nxyzt,
                          int memory,
                          int *zcbatch, 
                          int *znorm1, int *znorm2,
                          int *znorm3, int *znorm4,
                          int *zrho12, int *zrho34,
                          int *zp, int *zpx, int *zpy, int *zpz,
                          int *zscpk2, int *zq, int *zqx,
                          int *zqy, int *zqz, int *zscqk2)
{
    int mij;
    int mkl;
    int zone3;
    int zone4;
    int zone12;
    int ncsize;
   
    zone3 = npgto1 + npgto2 + npgto3 + npgto4 + nij + nkl;
    ncsize = nxyzt;
    
/*             ...if the MEMORY keyword is activated, the routine */
/*                will only determine the optimum and minimum flp */
/*                memory (in that order), place them respectively */
/*                into the NKLBLK and NIJBLK variables and exit. */
    mij = nij;
    mkl = nkl;
    if (memory)
    {
        zone12 = ncsize;
        zone4 = (mij + mkl) * 5;
        *zcbatch = zone12 + zone3 + zone4;
    }
    else
    {
/*             ...the actual fitting into the maximum memory given. */ 
        zone12 = ncsize;
        zone4 = (mij + mkl) * 5;
        assert (zone12 + zone3 + zone4 <= zmax);
                     
/*             ...generate the memory allocation pointers. */
        *zcbatch = 0;
        *znorm1 = *zcbatch + ncsize;
        *znorm2 = *znorm1 + npgto1;
        *znorm3 = *znorm2 + npgto2;
        *znorm4 = *znorm3 + npgto3;
        *zrho12 = *znorm4 + npgto4;
        *zrho34 = *zrho12 + nij;
        *zp = *zrho34 + nkl;
        *zpx = *zp + mij;
        *zpy = *zpx + mij;
        *zpz = *zpy + mij;
        *zscpk2 = *zpz + mij;
        *zq = *zscpk2 + mij;
        *zqx = *zq + mkl;
        *zqy = *zqx + mkl;
        *zqz = *zqy + mkl;
        *zscqk2 = *zqz + mkl;
    }

    return 0;
}
