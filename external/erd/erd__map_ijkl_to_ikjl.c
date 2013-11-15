#include <stdio.h>
#include <stdlib.h>


int erd__map_ijkl_to_ikjl (int ni, int nj, int nk, int nl,
                           double *x, double *y)
{
    int i;
    int j;
    int k;
    int l;

    for (l = 0; l < nl; l++)
    {
        for (k = 0; k < nk; k++)
        {
            for (j = 0; j < nj; j++)
            {
                for (i = 0; i < ni; i++)
                {
                    // Y(I,K,J,L) = X(I,J,K,L)
                    y[l * nj * nk * ni + 
                      j * nk * ni +
                      k * ni + i] = x[l * nk * nj * ni + 
                                      k * nj * ni +
                                      j * ni + i];
                }
            }
        }
    }

    return 0;
}


int erd__map_ijkl_to_ikjl_ (int *ni, int *nj, int *nk, int *nl,
                            int *tile, double *x, double *y)
{
    erd__map_ijkl_to_ikjl (*ni, *nj, *nk, *nl, x, y);
    
    return 0;
}