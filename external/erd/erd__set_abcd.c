#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "erd.h"
#include "erdutil.h"

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__SET_ABCD */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine handles the logistics on how to evaluate */
/*                the (12|34) integral batch in the most efficient way. */
/*                It performs the label map: */
/*                             (12|34) --> (AB|CD) */
/*                according to certain criteria that have to be met for */
/*                efficiency. The freedom we have in making the internal */
/*                association 1,2,3,4 -> A,B,C,D follows from the 8-fold */
/*                permutational symmetry of the integrals in (12|34): */
/*                     (12|34) = (21|34) = (12|43) = (21|43) = */
/*                     (34|12) = (43|12) = (34|21) = (43|21) */
/*                where the first line has only 12 -> 21 and 34 -> 43 */
/*                switches and the second line is obtained from the */
/*                first by bra <-> ket transpositions. */
/*                The type of switch to be applied is simply governed */
/*                by the demand that the final A,B,C,D shell labels */
/*                obey the relations A>=B and C>=D, since this means */
/*                the least amount of work (# of steps) for the HRR */
/*                procedure. The decision to perform a bra <-> ket */
/*                transposition comes from handling memory allocations */
/*                during both HRR on the bra and ket sides. */
/*                  Input (x = 1,2,3 and 4): */
/*                   NCGTOx        =  # of contractions for csh x */
/*                   NPGTOx        =  # of primitives per contraction */
/*                                    for csh x */
/*                   SHELLx        =  the shell type for csh x */
/*                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers */
/*                                    y = 1,2,3 and 4 */
/*                   EXPx          =  primitive exponents for csh x */
/*                   CCx           =  contraction coeffs for csh x */
/*                   SPHERIC       =  is true, if spherical integrals */
/*                                    are wanted, false if cartesian */
/*                                    ones are wanted */
/*                  Output (x = A,B,C and D): */
/*                   NCGTOx        =  # of contractions for csh x */
/*                   NPGTOx        =  # of primitives per contraction */
/*                                    for csh x */
/*                   SHELLx        =  the shell type for csh x */
/*                   SHELLy        =  the shell sums: y = P,Q,T = */
/*                                    A+B,C+D,P+Q */
/*                   MXSHELL       =  the largest (maximum) shell type */
/*                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers */
/*                                    y = A,B,C and D */
/*                   ATOMIC        =  indicates, if purely atomic */
/*                                    integrals will be evaluated */
/*                   ATOMxy        =  indicates, if atoms x and y are */
/*                                    considered to be equal for the */
/*                                    pairs xy = AB and CD */
/*                   EQUALxy       =  indicates, if csh x and csh y are */
/*                                    considered to be equal for the */
/*                                    pairs xy = AB and CD */
/*                   xxX,xxY,xxZ   =  the x,y,z-coordinate differences */
/*                                    for centers xx = AB and CD */
/*                   NxxCOOR       =  # of non-zero x,y,z-coordinate */
/*                                    differences for centers xx = AB */
/*                                    and CD */
/*                   RNxxSQ        =  square of the magnitude of the */
/*                                    distance between centers xx = AB */
/*                                    and CD */
/*                   SPNORM        =  normalization factor due to */
/*                                    presence of s- and p-type shells. */
/*                                    For each s-type shell there is a */
/*                                    factor of 1 and for each p-type */
/*                                    shell a factor of 2 */
/*                   NXYZx         =  # of cartesian monomials for csh x */
/*                   NXYZE(F)T     =  sum of # of cartesian monomials */
/*                                    for all shells in the range */
/*                                    E = A,...,A+B and in the range */
/*                                    F = C,...,C+D */
/*                   NXYZy         =  # of cartesian monomials for */
/*                                    y = P,Q shells */
/*                   NRYx          =  # of spherical functions for csh x */
/*                   INDEXx        =  index A,B,C,D -> 1,2,3,4 map */
/*                   SWAPxy        =  is .true. for xy = 12 and 34, if */
/*                                    a swap 1 <-> 2 and 3 <-> 4 has */
/*                                    been performed */
/*                   SWAPRS(TU)    =  is set .true. if the contraction */
/*                                    order of the primitives pair AB(CD) */
/*                                    will be performed in reverse order */
/*                                    BA(DC) for efficiency reasons */
/*                   TR1234        =  is .true., if a bra <-> ket */
/*                                    transposition has been applied */
/*                   LEXPx         =  pointers to locate appropriate */
/*                                    section of the exponent array */
/*                                    corresponding to csh x */
/*                   LCCx          =  pointers to locate appropriate */
/*                                    section of the contraction coeff */
/*                                    array corresponding to csh x */
/*                   LCCSEGx       =  pointers to locate appropriate */
/*                                    section of the lowest and highest */
/*                                    primitive index array defining */
/*                                    segmented contraction boundaries */
/*                                    for csh x */
/*                   NXYZHRR       =  maximum dimension of cartesian and */
/*                                    spherical components during the */
/*                                    entire HRR contraction procedure */
/*                   NCOLHRR       =  maximum # of HRR rotation matrix */
/*                                    columns needed to generate the */
/*                                    final HRR rotation matrix */
/*                   NROTHRR       =  maximum # of HRR rotation matrix */
/*                                    elements needed to generate the */
/*                                    final HRR rotation matrix */
/*                   EMPTY         =  int flag, indicating if an */
/*                                    empty batch of integrals is */
/*                                    expected. */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__set_abcd(
    uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
    uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
    bool atomic,
    double x1, double y1, double z1,
    double x2, double y2, double z2, 
    double x3, double y3, double z3,
    double x4, double y4, double z4, bool spheric,
    uint32_t *restrict npgtoa_ptr, uint32_t *restrict npgtob_ptr, uint32_t *restrict npgtoc_ptr, uint32_t *restrict npgtod_ptr,
    uint32_t *restrict shella_ptr, uint32_t *restrict shellb_ptr, uint32_t *restrict shellc_ptr, uint32_t *restrict shelld_ptr,
    double *restrict xa_ptr, double *restrict ya_ptr, double *restrict za_ptr,
    double *restrict xb_ptr, double *restrict yb_ptr, double *restrict zb_ptr,
    double *restrict xc_ptr, double *restrict yc_ptr, double *restrict zc_ptr,
    double *restrict xd_ptr, double *restrict yd_ptr, double *restrict zd_ptr,
    uint32_t *restrict nxyza_ptr, uint32_t *restrict nxyzb_ptr, uint32_t *restrict nxyzc_ptr, uint32_t *restrict nxyzd_ptr,
    uint32_t *restrict nxyzet_ptr, uint32_t *restrict nxyzft_ptr,
    uint32_t *restrict nrya_ptr, uint32_t *restrict nryb_ptr, uint32_t *restrict nryc_ptr, uint32_t *restrict nryd_ptr,
    uint32_t *restrict nabcoor_ptr, uint32_t *restrict ncdcoor_ptr,
    uint32_t *restrict ncolhrr, uint32_t *restrict nrothrr,
    uint32_t *restrict nxyzhrr, bool *restrict empty, bool *restrict tr1234)
{
/*             ...generate all 1,2,3,4 data. Decide as early as */
/*                possible, if a zero batch of integrals is expected. */
    *empty = false;

    const uint32_t preshellp = shell1 + shell2;
    const uint32_t preshellq = shell3 + shell4;
    const uint32_t shellt = preshellp + preshellq;

    const uint32_t mxshell = max4x32u(shell1, shell2, shell3, shell4);
    const bool case1 = (shellt % 2 == 1);
    const bool case2 = spheric && (2 * mxshell > shellt);
    if (atomic && (case1 || case2)) {
        *empty = true;
        return;
    }

/*             ...determine csh equality between center pairs 1,2 */
/*                and 3,4 in increasing order of complexity: */
/*                 centers -> shells -> exponents -> ctr coefficients */
/*             ...set the cartesian and spherical dimensions. In case */
/*                no spherical transformations are wanted, set the */
/*                corresponding dimensions equal to the cartesian ones. */
    uint32_t nxyz1 = ((shell1 + 1) * (shell1 + 2)) / 2;
    uint32_t nxyz2 = ((shell2 + 1) * (shell2 + 2)) / 2;
    uint32_t nxyz3 = ((shell3 + 1) * (shell3 + 2)) / 2;
    uint32_t nxyz4 = ((shell4 + 1) * (shell4 + 2)) / 2;
    uint32_t nry1 = 2 * shell1 + 1;
    uint32_t nry2 = 2 * shell2 + 1;
    uint32_t nry3 = 2 * shell3 + 1;
    uint32_t nry4 = 2 * shell4 + 1;
    if (!spheric) {
        nry1 = nxyz1;
        nry2 = nxyz2;
        nry3 = nxyz3;
        nry4 = nxyz4;
    }

/*             ...decide on the 1 <-> 2 and/or 3 <-> 4 swapping. */

/*             ...calculate NXYZHRR for the two possible HRR */
/*                and (if any) cartesian -> spherical transformation */
/*                sequences: */
/*                    i) initial dimension:  NXYZE * NXYZF */
/*                   ii) perform HRR on 34:  NXYZE * NXYZ3 * NXYZ4 */
/*                  iii) cart -> sph on 34:  NXYZE * NRY3 * NRY4 */
/*                   iv) perform HRR on 12:  NXYZ1 * NXYZ2 * NRY3 * NRY4 */
/*                    v) cart -> sph on 12:  NRY1 * NRY2 * NRY3 * NRY4 */
/*                    i) initial dimension:  NXYZE * NXYZF */
/*                   ii) perform HRR on 12:  NXYZ1 * NXYZ2 * NXYZF */
/*                  iii) cart -> sph on 12:  NRY1 * NRY2 * NXYZF */
/*                   iv) perform HRR on 34:  NRY1 * NRY2 * NXYZ3 * NXYZ4 */
/*                    v) cart -> sph on 34:  NRY1 * NRY2 * NRY3 * NRY4 */
/*                The only dimension increasing steps are ii) and iv) */
/*                in both cases. Hence we first find the maximum */
/*                between ii) and iv) for both sequences and then we */
/*                take the overall minimum of these two maxima. */
/*                Since the order of sequence of the HRR on the A,B,C,D */
/*                labels is CD followed by AB, the overall minimum */
/*                will decide if to perform a bra <-> ket transposition */
/*                on the 12/34 labels. */
    const uint32_t preshella = max32u(shell1, shell2);
    const uint32_t preshellb = min32u(shell1, shell2);
    const uint32_t preshellc = max32u(shell3, shell4);
    const uint32_t preshelld = min32u(shell3, shell4);
    const uint32_t nxyze = (preshellp + 1) * (preshellp + 2) * (preshellp + 3) / 6 - preshella * (preshella + 1) * (preshella + 2) / 6;
    const uint32_t nxyzf = (preshellq + 1) * (preshellq + 2) * (preshellq + 3) / 6 - preshellc * (preshellc + 1) * (preshellc + 2) / 6;

    if ((preshellb | preshelld) == 0) {
        *nxyzhrr = nxyze * nxyzf;
        *tr1234 = false;
    } else {
        const uint32_t nhrr1st = max32u(nxyze * nxyz3 * nxyz4, nxyz1 * nxyz2 * nry3 * nry4);
        const uint32_t nhrr2nd = max32u(nxyzf * nxyz1 * nxyz2, nxyz3 * nxyz4 * nry1 * nry2);
        *nxyzhrr = min32u(nhrr1st, nhrr2nd);
        *tr1234 = nhrr1st > nhrr2nd;
    }

/*             ...according to the previously gathered info, set the */
/*                new A,B,C,D shells, # of primitives + contraction */
/*                coeffs as well as pointers to the alpha exponents */
/*                and contraction coefficients. Also set the info for */
/*                evaluation of the [e0|f0] batches and for the HRR */
/*                steps later on. */
    uint32_t nxyzet = nxyze, nxyzft = nxyzf;
    if (*tr1234) {
        ERD_SWAP(nxyzet, nxyzft);

        ERD_SWAP(x1, x3);
        ERD_SWAP(y1, y3);
        ERD_SWAP(z1, z3);
        ERD_SWAP(shell1, shell3);
        ERD_SWAP(npgto1, npgto3);
        ERD_SWAP(nxyz1, nxyz3);
        ERD_SWAP(nry1, nry3);

        ERD_SWAP(x2, x4);
        ERD_SWAP(y2, y4);
        ERD_SWAP(z2, z4);
        ERD_SWAP(shell2, shell4);
        ERD_SWAP(npgto2, npgto4);
        ERD_SWAP(nxyz2, nxyz4);
        ERD_SWAP(nry2, nry4);
    }
    *nxyzet_ptr = nxyzet;
    *nxyzft_ptr = nxyzft;
    
    double xa = x1, xb = x2, xc = x3, xd = x4;
    double ya = y1, yb = y2, yc = y3, yd = y4;
    double za = z1, zb = z2, zc = z3, zd = z4;
    uint32_t shella = shell1, shellb = shell2, shellc = shell3, shelld = shell4;
    uint32_t npgtoa = npgto1, npgtob = npgto2, npgtoc = npgto3, npgtod = npgto4;
    uint32_t nxyza = nxyz1, nxyzb = nxyz2, nxyzc = nxyz3, nxyzd = nxyz4;
    uint32_t nrya = nry1, nryb = nry2, nryc = nry3, nryd = nry4;

    if (shell1 < shell2) {
        ERD_SWAP(xa, xb);
        ERD_SWAP(ya, yb);
        ERD_SWAP(za, zb);
        ERD_SWAP(shella, shellb);
        ERD_SWAP(npgtoa, npgtob);
        ERD_SWAP(nxyza, nxyzb);
        ERD_SWAP(nrya, nryb);
    }
    if (shell3 < shell4) {
        ERD_SWAP(xc, xd);
        ERD_SWAP(yc, yd);
        ERD_SWAP(zc, zd);
        ERD_SWAP(shellc, shelld);
        ERD_SWAP(npgtoc, npgtod);
        ERD_SWAP(nxyzc, nxyzd);
        ERD_SWAP(nryc, nryd);
    }
    *shella_ptr = shella;
    *shellb_ptr = shellb;
    *shellc_ptr = shellc;
    *shelld_ptr = shelld;
    *npgtoa_ptr = npgtoa;
    *npgtob_ptr = npgtob;
    *npgtoc_ptr = npgtoc;
    *npgtod_ptr = npgtod;
    *nxyza_ptr = nxyza;
    *nxyzb_ptr = nxyzb;
    *nxyzc_ptr = nxyzc;
    *nxyzd_ptr = nxyzd;
    *nrya_ptr = nrya;
    *nryb_ptr = nryb;
    *nryc_ptr = nryc;
    *nryd_ptr = nryd;
    *xa_ptr = xa;
    *xb_ptr = xb;
    *xc_ptr = xc;
    *xd_ptr = xd;
    *ya_ptr = ya;
    *yb_ptr = yb;
    *yc_ptr = yc;
    *yd_ptr = yd;
    *za_ptr = za;
    *zb_ptr = zb;
    *zc_ptr = zc;
    *zd_ptr = zd;
    uint32_t nabcoor = 3;
    if (xa == xb) {
        nabcoor--;
    }
    if (ya == yb) {
        nabcoor--;
    }
    if (za == zb) {
        nabcoor--;
    }
    uint32_t ncdcoor = 3;
    if (xc == xd) {
        ncdcoor--;
    }
    if (yc == yd) {
        ncdcoor--;
    }
    if (zc == zd) {
        ncdcoor--;
    }
    const uint32_t shellp = shella + shellb;
    const uint32_t shellq = shellc + shelld;
    const uint32_t nxyzp = (shellp + 1) * (shellp + 2) / 2;
    const uint32_t nxyzq = (shellq + 1) * (shellq + 2) / 2;
    
    *ncolhrr = 0;
    *nrothrr = 0;
    if (shellb != 0) {
        const uint32_t ngh = nxyzet;
        uint32_t nxyzgo = nxyzet;
        uint32_t nxyzho = 1;
        uint32_t nxyzi = nxyzp;
        uint32_t shellg = shellp;
        uint32_t nrow = 1;
        uint32_t ncol = ngh;
        uint32_t nrot = ngh;
        for (uint32_t shellh = 1; shellh <= shellb; shellh++) {
            nxyzgo -= nxyzi;
            nxyzho += shellh + 1;
            const uint32_t ngho = nxyzgo * nxyzho;
            switch (nabcoor) {
                case 3:
                {
                    const uint32_t m = shellh / 3 + 1;
                    nrow += m * (m + ((shellh % 3) == 2));
                    break;
                }
                case 2:
                    nrow += shellh / 2 + 1;
                    break;
                case 1:
                    ++nrow;
                    break;
            }
            ncol = max32u(ngho, ncol);
            nrot = max32u(nrow * ngho, nrot);
            nxyzi = nxyzi - shellg - 1;
            --shellg;
        }
        *ncolhrr = ncol;
        *nrothrr = nrot;
    }
/*             ...next find maximum values for the HRR on the CD-part */
/*                and set overall maximum values. */
    if (shelld != 0) {
        const uint32_t ngh = nxyzft;
        uint32_t nxyzgo = nxyzft;
        uint32_t nxyzho = 1;
        uint32_t nxyzi = nxyzq;
        uint32_t shellg = shellq;
        uint32_t nrow = 1;
        uint32_t ncol = ngh;
        uint32_t nrot = ngh;
        for (uint32_t shellh = 1; shellh <= shelld; shellh++) {
            nxyzgo -= nxyzi;
            nxyzho += shellh + 1;
            const uint32_t ngho = nxyzgo * nxyzho;
            switch (ncdcoor) {
                case 3:
                {
                    const uint32_t m = shellh / 3 + 1;
                    nrow += m * (m + ((shellh % 3) == 2));
                    break;
                }
                case 2:
                    nrow += shellh / 2 + 1;
                    break;
                case 1:
                    nrow += 1;
                    break;
            }
            ncol = max32u(ngho, ncol);
            nrot = max32u(nrow * ngho, nrot);
            nxyzi -= shellg + 1;
            --shellg;
        }
        *ncolhrr = max32u(ncol, *ncolhrr);
        *nrothrr = max32u(nrot, *nrothrr);
    }
    *nabcoor_ptr = nabcoor;
    *ncdcoor_ptr = ncdcoor;
}
