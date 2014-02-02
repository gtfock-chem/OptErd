#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "erd.h"

#pragma offload_attribute(push, target(mic))

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__INT2D_TO_E0F0 */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine assembles the set of batches of cartesian */
/*                eris [E0|F0] , E = A to P, F = C to Q, adding up all */
/*                the contributions from all the 2D PQ integrals. */
/*                The routine uses the reduced Rys multiplication scheme */
/*                as suggested in R.Lindh, U.Ryu and B.Liu, JCP 95, 5889. */
/*                This scheme reuses intermediate products between */
/*                2DX and 2DY integrals, which can be achieved by having */
/*                the outer loops run over all possible x and y monomial */
/*                parts and the inner loops over all allowed E and F */
/*                shell combinations. The price to pay for such loop */
/*                ordering is the scattered addressing of locations */
/*                within the batch array, which has its rows and columns */
/*                ordered such that the E and F shells are increasing */
/*                and within each E and F shell the monomials are */
/*                ordered such that x>y>z in exponents. */
/*                An example follows: */
/*                ------------------- */
/*                Let E = 0,2 and F = 0,1. Then we have the left and */
/*                right hand of the batch array ordered as follows: */
/*                          left xyz        right xyz */
/*                             000             000 */
/*                             ---             --- */
/*                             100             100 */
/*                             010             010 */
/*                             001             001 */
/*                             --- */
/*                             200 -> -5 */
/*                             110 */
/*                             101 -> -3 */
/*                             020 */
/*                             011 */
/*                             002 -> 0 */
/*                The batch would thus have dimensions 10 x 4. For */
/*                the left side (and analogous for the right side) the */
/*                reduced multiplication scheme would have its outer */
/*                most loop run over x=0,2, followed by the next loop */
/*                y=0,2-x. The innermost loop would then run over the */
/*                allowed shells E=E(max),max(E(min),x+y). In this */
/*                case all x,y-pairs can be reused for all appropriate */
/*                shell combinations. */
/*                To find the address of a specific x,y,z,E combination */
/*                inside the batch array, we first note that the z-part */
/*                is dependent on the x,y-parts and is hence not needed. */
/*                Lets look at the E-part first. The E-part is evaluated */
/*                from its dimension formula (E+1)*(E+2)/2. Organizing */
/*                the inner E-loop to run from E(max) always, the */
/*                dimension for E(max) is passed as the argument NXYZP */
/*                and all lower E dimensions are calculated by the */
/*                formula relating dimensions between E and E+1: */
/*                           dim(E+1) = dim(E) + E + 2 */
/*                In this way multiplications in the E-part can be */
/*                entirely avoided. The x,y-part is defined as the */
/*                part which has to be subtracted from dim(E) to */
/*                reach the xyz monomial position inside the E shell. */
/*                It can be divided into an x-part and a y-part. The */
/*                x-part is given by the fomula: */
/*                          x-part = - x*E + x(x-3)/2 */
/*                and for the example above has been given for E=2 */
/*                and marked with arrows ->. The last term of the x-part */
/*                involves 1 multiplication and division, however it */
/*                can be changed to: */
/*                                               x-1 */
/*                          x-part = - x*E - x + sum i */
/*                                               i=0 */
/*                and clever additions inside the x-loop avoid the use */
/*                of multiplications and divisions. The y-part is trivial */
/*                and is simply equal to -y. The overall conclusion is */
/*                thus that the location of a specific x,y,z,E quadruple */
/*                inside the batch comes at the cost of one x*E(max) */
/*                multiplication in the outermost x-loops, since the */
/*                other x*E ones can again be reached via stepwise */
/*                subtraction of x from x*E(max). */
/*                Due to the very computational intensive steps inside */
/*                the x,y,z,E loops, special sections of identical */
/*                x,y,z,E loop structers have been given for each */
/*                # of roots =< 9, thus saving considerable computing */
/*                time over the general case. */
/*                For comments on how the x,y,z,E loop structures are */
/*                coded please refer to the general root case. */
/*                  Input: */
/*                    SHELLx      =  shell types for individual csh */
/*                                   x=A,C and csh sums P=A+B,Q=C+D */
/*                    NGQP        =  # of gaussian quadrature points */
/*                                   (roots) */
/*                    NEXQ        =  current # of exponent quadruplets */
/*                    NGQEXQ      =  product of # of gaussian quadrature */
/*                                   points times exponent quadruplets */
/*                    NXYZE(F)T   =  sum of # of cartesian monomials */
/*                                   for all shells in the range */
/*                                   E = A,...,P=A+B and in the range */
/*                                   F = C,...,Q=C+D */
/*                    NXYZy       =  # of cartesian monomials for */
/*                                   y = P,Q shells */
/*                    INT2Dx      =  all current 2D PQ integrals for */
/*                                   each cartesian component */
/*                                   (x = X,Y,Z) */
/*                    TEMP1(2)    =  scratch arrays holding intermediate */
/*                                   2D PQ integral products */
/*                    SCALE       =  the NGQEXQ scaling factors */
/*                  Output: */
/*                    BATCH       =  batch of primitive cartesian */
/*                                   [E0|F0] integrals corresponding */
/*                                   to all current exponent quadruplets */
/* ------------------------------------------------------------------------ */
int erd__int2d_to_e0f0 (int shella, int shellp, int shellc, int shellq,
                        int ngqexq,
                        int nxyzet, int nxyzft, int nxyzp, int nxyzq,
                        double *int2dx, double *int2dy, double *int2dz,
                        double *temp1, double *temp2, double *scale,
                        double *batch)
{
    int i, j, k, m, se, sf, xe, ye, ze, xf, yf, zf, xep, xfp;
    int xye, xyf, xyep, xyfp, seend, sfend, yeend, yfend, xemax,
        xfmax, nxyze, nxyzf;
  
    int indices[ngqexq * nxyzet * nxyzft][4];
    xfp = nxyzft + 2;
    int indx, indy, indz, indb;
    k = 0;
    for (xf = 0; xf <= shellq; ++xf)
    {
        xfp = xfp + xf - 2;
        xfmax = xf * shellq;
        yfend = shellq - xf;
        xep = nxyzet + 2;
        for (xe = 0; xe <= shellp; ++xe)
        {
            xep = xep + xe - 2;
            xemax = xe * shellp;
            yeend = shellp - xe;
            indx = (xe + xf * (shellp + 1)) * ngqexq;

/*             ...middle loops over y,y-pairs. Skip multiplication */
/*                of y,y-contributions, if we have a 0,0-pair, as */
/*                then the 2DY integral is equal to 1. */
            xyfp = xfp - xfmax;
            for (yf = 0; yf <= yfend; ++yf)
            {
                xyf = xf + yf;
                --xyfp;
                sfend = MAX (shellc, xyf);
                xyep = xep - xemax;
                for (ye = 0; ye <= yeend; ++ye)
                {
                    xye = xe + ye;
                    --xyep;
                    seend = MAX (shella, xye);
                    indy = (ye + yf * (shellp + 1)) * ngqexq;

                    j = xyfp;
                    nxyzf = nxyzq;
                    for (sf = shellq; sf >= sfend; --sf)
                    {
                        zf = sf - xyf;
                        i = xyep;
                        nxyze = nxyzp;
                        for (se = shellp; se >= seend; --se)
                        {
                            ze = se - xye;
                            indz = (ze + zf * (shellp + 1)) * ngqexq;

                            indices[k][0] = indx;
                            indices[k][1] = indy;
                            indices[k][2] = indz;
                            indices[k][3] = i + j * nxyzet;
                            k++;
                             
                            i = i - nxyze + xe;
                            nxyze = nxyze - se - 1;
                        }
                        j = j - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
                    }
                }
            }
        }
    }

    int k1, m1;
    for(k1 = 0; k1 < k; k1++)
    {
        indx = indices[k1][0];
        indy = indices[k1][1];
        indz = indices[k1][2];
        indb = indices[k1][3];

        double sum = 0;
        for(m = 0; m < ngqexq; m+=SIMD_WIDTH)
        {
#pragma vector aligned
            for(m1 = 0; m1 < SIMD_WIDTH; m1++)
            {
                sum += scale[m + m1]
                    * int2dx[m + m1 + indx]
                    * int2dy[m + m1 + indy]
                    * int2dz[m + m1 + indz];
            }
        }
        batch[indb] = sum;
         
    }

    return 0;
}

#pragma offload_attribute(pop)
