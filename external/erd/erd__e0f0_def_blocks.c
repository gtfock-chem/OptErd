#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "erd.h"


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
/*                  --- Zone 2 --- */
/*                   ZPBATCH = offset for (blocked) [e0|f0] batch (K4) */
/*                  --- Zone 3 --- */
/*                   ZNORMA = offset for all A norms */
/*                   ZNORMB = offset for all B norms */
/*                   ZNORMC = offset for all C norms */
/*                   ZNORMD = offset for all D norms */
/*                   ZRHOAB = offset for all AB exponential prefactors */
/*                   ZRHOCD = offset for all CD exponential prefactors */
/*                  --- Zone 4: for block [e0|f0] generation only --- */
/*                   ZP = offset for (blocked) P values (K2) */
/*                   ZPX = offset for (blocked) PX values (K2) */
/*                   ZPY = offset for (blocked) PY values (K2) */
/*                   ZPZ = offset for (blocked) PZ values (K2) */
/*                   ZPAX = offset for (blocked) PAX values (K2) */
/*                   ZPAY = offset for (blocked) PAY values (K2 */
/*                   ZPAZ = offset for (blocked) PAZ values (K2) */
/*                   ZPINVHF = offset for (blocked) 1/2P values (K2) */
/*                   ZSCPK2 = offset for (blocked) P scale values (K2) */
/*                   ZQ = offset for (blocked) Q values (K2) */
/*                   ZQX = offset for (blocked) QX values (K2) */
/*                   ZQY = offset for (blocked) QY values (K2) */
/*                   ZQZ = offset for (blocked) QZ values (K2) */
/*                   ZQCX = offset for (blocked) QCX values (K2) */
/*                   ZQCY = offset for (blocked) QCY values (K2 */
/*                   ZQCZ = offset for (blocked) QCZ values (K2) */
/*                   ZQINVHF = offset for (blocked) 1/2Q values (K2) */
/*                   ZSCQK2 = offset for (blocked) Q scale values (K2) */
/*                   ZRTS = offset for (blocked) quad roots (K4) */
/*                   ZWTS = offset for (blocked) quad weights (K4) */
/*                   ZGQSCR = offset for (blocked) quad scratch (K4) */
/*                   ZTVAL = offset for (blocked) T exponents (K4) */
/*                   ZPQPINV = offset for (blocked) 1/(P+Q) values (K4) */
/*                   ZSCPQK4 = offset for (blocked) PQ scale values (K4) */
/*                   ZB00 = offset for (blocked) VRR coeff B00 (K4) */
/*                   ZB01 = offset for (blocked) VRR coeff B01 (K4) */
/*                   ZB10 = offset for (blocked) VRR coeff B10 (K4) */
/*                   ZC00X = offset for (blocked) VRR coeff C00X (K4) */
/*                   ZC00Y = offset for (blocked) VRR coeff C00Y (K4) */
/*                   ZC00Z = offset for (blocked) VRR coeff C00Z (K4) */
/*                   ZD00X = offset for (blocked) VRR coeff D00X (K4) */
/*                   ZD00Y = offset for (blocked) VRR coeff D00Y (K4) */
/*                   ZD00Z = offset for (blocked) VRR coeff D00Z (K4) */
/*                   ZINT2DX = offset for (blocked) 2DX integrals (K4) */
/*                   ZINT2DY = offset for (blocked) 2DY integrals (K4) */
/*                   ZINT2DZ = offset for (blocked) 2DZ integrals (K4) */
/*                  --- Zone 4: for contraction only --- */
/*                   ZWORK = offset for contraction working array */
/* ------------------------------------------------------------------------ */
int erd__e0f0_def_blocks (int zmax, int npgtoa, int npgtob,
                          int npgtoc, int npgtod,
                          int shellp, int shellq,
                          int nij, int nkl, int nrs, int ntu,
                          int nrstu, int ngqp, int ngqscr,
                          int nxyzt, int l1cache, int nctrow,
                          int memory,
                          int *nijblk, int *nklblk,
                          int *npsize, int *ncsize, int *nwsize,
                          int *nint2d, int *mxprim, int *mnprim,
                          int *zcbatch, int *zpbatch, int *zwork,
                          int *znorma, int *znormb, int *znormc,
                          int *znormd, int *zrhoab, int *zrhocd,
                          int *zp, int *zpx, int *zpy,
                          int *zpz, int *zpax, int *zpay,
                          int *zpaz, int *zpinvhf, int *zscpk2,
                          int *zq, int *zqx, int *zqy,
                          int *zqz, int *zqcx, int *zqcy,
                          int *zqcz, int *zqinvhf, int *zscqk2,
                          int *zrts, int *zwts, int *zgqscr,
                          int *ztval, int *zpqpinv, int *zscpqk4,
                          int *zb00, int *zb01, int *zb10,
                          int *zc00x, int *zc00y, int *zc00z,
                          int *zd00x, int *zd00y, int *zd00z,
                          int *zint2dx, int *zint2dy,
                          int *zint2dz)
{
    int mij;
    int mkl;
    int wcol;
    int zone3;
    int zone4;
    int mijkl;
    int zone12;
    int mrskl;
    int zone4b;
    int zone4c;
    int mgqijkl;
    int npminrs;
    int npmintu;

    mij = nij;
    mkl = nkl;
    npminrs = MIN (npgtoa, npgtob);
    npmintu = MIN (npgtoc, npgtod);
    *ncsize = nxyzt * nrstu;
    *mnprim = MAX (npminrs, npmintu);    
    *mxprim = MAX (npgtoa, npgtob);
    *mxprim = MAX (*mxprim, npgtoc);
    *mxprim = MAX (*mxprim, npgtod);    
    zone3 = npgtoa + npgtob + npgtoc + npgtod + nij + nkl;


/*             ...if the MEMORY keyword is activated, the routine */
/*                will only determine the optimum and minimum flp */
/*                memory (in that order), place them respectively */
/*                into the NKLBLK and NIJBLK variables and exit. */
    if (memory)
    {
        mijkl = mij * mkl;
        mgqijkl = ngqp * mijkl;
        *nint2d = mgqijkl * (shellp + 1) * (shellq + 1);
        mrskl = nrs * mkl;
        *npsize = nxyzt * MAX (mijkl, mrskl);
        *nwsize = *npsize;
        zone12 = *npsize + *ncsize;
        zone4b = ngqscr + mijkl * 2 +
            (mij + mkl) * 9 + mgqijkl * 12 +
            *nint2d * 3;
        wcol = *mnprim;
        zone4c = *nwsize + nctrow * wcol;
        zone4 = MAX (zone4b, zone4c);
        *nklblk = zone12 + zone3 + zone4;
        
        mij = 1;
        mkl = 1;
        mijkl = mij * mkl;
        mgqijkl = ngqp * mijkl;
        *nint2d = mgqijkl * (shellp + 1) * (shellq + 1);
        mrskl = nrs * mkl;
        *npsize = nxyzt * MAX (mijkl, mrskl);
        *nwsize = *npsize;
        zone12 = *npsize + *ncsize;
        zone4b = ngqscr + mijkl * 2 + (mij + mkl) * 9 +
            mgqijkl * 12 + *nint2d * 3;
        wcol = *mnprim;
        zone4c = *nwsize + nctrow * wcol;
        zone4 = MAX (zone4b, zone4c);
        *nijblk = zone12 + zone3 + zone4;
    }
    else
    {
        mijkl = mij * mkl;
        mgqijkl = ngqp * mijkl;
        *nint2d = mgqijkl * (shellp + 1) * (shellq + 1);
        mrskl = nrs * mkl;
        *npsize = nxyzt * MAX (mijkl, mrskl);
        *nwsize = *npsize;
        zone12 = *npsize + *ncsize;
        zone4b = ngqscr + mijkl * 2 +
                 (mij + mkl) * 9 + mgqijkl * 12 +
                 *nint2d * 3;
        wcol = *mnprim;
        zone4c = *nwsize + nctrow * wcol;
        zone4 = MAX (zone4b, zone4c);
        assert (zone12 + zone3 + zone4 <= zmax);

        *nijblk = mij;
        *nklblk = mkl;
        *nwsize = zmax - zone12 - zone3;
/*             ...generate the memory allocation pointers. */
        *zcbatch = 1;
        *zpbatch = *ncsize + 1;
        *znorma = *zpbatch + *npsize;
        *znormb = *znorma + npgtoa;
        *znormc = *znormb + npgtob;
        *znormd = *znormc + npgtoc;
        *zrhoab = *znormd + npgtod;
        *zrhocd = *zrhoab + nij;
        *zp = *zrhocd + nkl;
        *zpx = *zp + mij;
        *zpy = *zpx + mij;
        *zpz = *zpy + mij;
        *zpax = *zpz + mij;
        *zpay = *zpax + mij;
        *zpaz = *zpay + mij;
        *zpinvhf = *zpaz + mij;
        *zscpk2 = *zpinvhf + mij;
        *zq = *zscpk2 + mij;
        *zqx = *zq + mkl;
        *zqy = *zqx + mkl;
        *zqz = *zqy + mkl;
        *zqcx = *zqz + mkl;
        *zqcy = *zqcx + mkl;
        *zqcz = *zqcy + mkl;
        *zqinvhf = *zqcz + mkl;
        *zscqk2 = *zqinvhf + mkl;
        *zrts = *zscqk2 + mkl;
        *zwts = *zrts + mgqijkl;
        *zgqscr = *zwts + mgqijkl;
        *ztval = *zgqscr + ngqscr;
        *zpqpinv = *ztval + mijkl;
        *zscpqk4 = *zpqpinv + mijkl;
        *zb00 = *zscpqk4 + mgqijkl;
        *zb01 = *zb00 + mgqijkl;
        *zb10 = *zb01 + mgqijkl;
        *zc00x = *zb10 + mgqijkl;
        *zc00y = *zc00x + mgqijkl;
        *zc00z = *zc00y + mgqijkl;
        *zd00x = *zc00z + mgqijkl;
        *zd00y = *zd00x + mgqijkl;
        *zd00z = *zd00y + mgqijkl;
        *zint2dx = *zd00z + mgqijkl;
        *zint2dy = *zint2dx + *nint2d;
        *zint2dz = *zint2dy + *nint2d;
        *zwork = *zp;
    }

    return 0;
}


int erd__e0f0_def_blocks_ (int * zmax, int * npgtoa,
                        int * npgtob, int * npgtoc, int * npgtod,
                        int * shellp, int * shellq, int * nij,
                        int * nkl, int * nrs, int * ntu,
                        int * nrstu, int * ngqp, int * ngqscr,
                        int * nxyzt, int * l1cache, int * nctrow,
                        int * memory, int * nijblk, int * nklblk,
                        int * npsize, int * ncsize, int * nwsize,
                        int * nint2d, int * mxprim, int * mnprim,
                        int * zcbatch, int * zpbatch, int * zwork,
                        int * znorma, int * znormb, int * znormc,
                        int * znormd, int * zrhoab, int * zrhocd,
                        int * zp, int * zpx, int * zpy,
                        int * zpz, int * zpax, int * zpay,
                        int * zpaz, int * zpinvhf, int * zscpk2,
                        int * zq, int * zqx, int * zqy,
                        int * zqz, int * zqcx, int * zqcy,
                        int * zqcz, int * zqinvhf, int * zscqk2,
                        int * zrts, int * zwts, int * zgqscr,
                        int * ztval, int * zpqpinv, int * zscpqk4,
                        int * zb00, int * zb01, int * zb10,
                        int * zc00x, int * zc00y, int * zc00z,
                        int * zd00x, int * zd00y, int * zd00z,
                        int * zint2dx, int * zint2dy,
                        int * zint2dz)
{
    erd__e0f0_def_blocks (*zmax, *npgtoa, *npgtob,
                          *npgtoc, *npgtod,
                          *shellp, *shellq,
                          *nij, *nkl, *nrs, *ntu,
                          *nrstu, *ngqp, *ngqscr,
                          *nxyzt, *l1cache, *nctrow,
                          *memory,
                          nijblk, nklblk,
                          npsize, ncsize, nwsize,
                          nint2d, mxprim, mnprim,
                          zcbatch, zpbatch, zwork,
                          znorma, znormb, znormc,
                          znormd, zrhoab, zrhocd,
                          zp, zpx, zpy,
                          zpz, zpax, zpay,
                          zpaz, zpinvhf, zscpk2,
                          zq, zqx, zqy,
                          zqz, zqcx, zqcy,
                          zqcz, zqinvhf, zscqk2,
                          zrts, zwts, zgqscr,
                          ztval, zpqpinv, zscpqk4,
                          zb00, zb01, zb10,
                          zc00x, zc00y, zc00z,
                          zd00x, zd00y, zd00z,
                          zint2dx, zint2dy,
                          zint2dz);

    return 0;
}