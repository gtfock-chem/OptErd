#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#include "../include/CInt.h"
#include "screening.h"


#pragma offload_attribute(push, target(mic))
ERD_t erd_mic;
double integrals[10000];
#pragma offload_attribute(pop)


int main (int argc, char **argv)
{
    int M;
    int N;
    int P;
    int Q;
    int nints;
    int ns;
    int nnz;
    double totalcalls = 0;
    int *shellptr;
    int *shellid;
    int *shellrid;
    double *shellvalue;
    struct timeval tv1;
    struct timeval tv2;
    double timepass;   
    double totalintls = 0;
    BasisSet_t basis;

    if (argc != 3)
    {
        printf ("Usage: %s <basisset> <xyz>\n", argv[0]);
        return 0;
    }

    FILE *ref_data_file = fopen ("ivalues.ref", "r");
    if (ref_data_file == NULL)
    {
        fprintf (stderr, "ivalues.ref does not exist\n");
        exit (0);
    }

    int errcount; 
    errcount = 0;
    // load basis set   
    CInt_createBasisSet (&basis);
    CInt_loadBasisSet (basis, argv[1], argv[2]);
    schwartz_screening (basis, &shellptr, &shellid, &shellrid, &shellvalue, &nnz);

    printf ("Molecule info:\n");
    printf ("  #Atoms\t= %d\n", CInt_getNumAtoms (basis));
    printf ("  #Shells\t= %d\n", CInt_getNumShells (basis));
    printf ("  #Funcs\t= %d\n", CInt_getNumFuncs (basis));
    printf ("  #OccOrb\t= %d\n", CInt_getNumOccOrb (basis));

    // Find the number of MIC cards available
    int num_devices = 0;   
#ifdef __INTEL_OFFLOAD
    num_devices = _Offload_number_of_devices();
#endif

    if(num_devices == 0)
    {
        printf("No target devices available. Exiting\n");
        exit(0);
    }
    else
    {
        printf("Number of Target devices installed: %d\n\n",num_devices);
    }

    printf ("Computing integrals ...\n");
    ns = CInt_getNumShells (basis);
    timepass = 0.0;
    int i;
    int j;
    int start1;
    int end1;
    int start2;
    int end2;
    double value1;
    double value2;

    for (i = 0; i < num_devices; i++)
    {
        #pragma offload target(mic:i) nocopy(basis_mic, erd_mic)
        {
            CInt_createERD (basis_mic, &erd_mic);
        }
    }

    int dimMax = CInt_getMaxShellDim (basis);
    int nints0;
    int k;
    double * integrals0 =
        (double *)malloc (sizeof(double) * dimMax * dimMax * dimMax * dimMax);
    assert (integrals0 != NULL);   
    
    for (M = 0; M < ns; M++)
    {
        start1 = shellptr[M];
        end1 = shellptr[M + 1];       
        for (i = start1; i < end1; i++)
        {
            N = shellid[i];
            value1 = shellvalue[i];
            for (P = 0; P < ns; P++)
            {
                start2 = shellptr[P];
                end2 = shellptr[P + 1];              
                for (j = start2; j < end2; j++)
                {
                    Q = shellid[j];
                    value2 = shellvalue[j];
                    if (M > N || P > Q || (M + N) > (P + Q))
                        continue;
                    if (fabs(value1 * value2) < TOLSRC * TOLSRC)
                        continue;                        
                    totalcalls = totalcalls + 1;
                    gettimeofday (&tv1, NULL);

                    #pragma offload target(mic:0) \
                            in(M, N, P, Q) \
                            nocopy(basis_mic, erd_mic) \
                            out(integrals, nints)
                    {
                        double *integrals_mic;
                        CInt_computeShellQuartet (basis_mic, erd_mic, M, N, P, Q,
                                &integrals_mic, &nints);
                        memcpy(integrals, integrals_mic, nints * sizeof(double));
                    }
                    gettimeofday (&tv2, NULL);
                    timepass += (tv2.tv_sec - tv1.tv_sec) +
                        (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;
                    totalintls = totalintls + nints;
                    
                    fread (&nints0, sizeof(int), 1, ref_data_file);
                    if (nints0 != 0)
                    {
                        fread (integrals0, sizeof(double), nints0, ref_data_file);
                    }
                    // compare results
                    if (nints == 0 && nints0 == 0)
                    {
                        continue;
                    }
                    else if (nints != 0 && nints0 == 0)
                    {
                        for (k = 0; k < nints; k++)
                        {
                            if (integrals[k] > 1e-10)
                            {
                                printf ("ERROR: %d %d %d %d: %le %le\n",
                                    M, N, P, Q, 0.0, integrals[k]);
                                errcount++;
                            }
                        }
                    }
                    else if (nints == 0 && nints != 0)
                    {
                        for (k = 0; k < nints0; k++)
                        {
                            if (integrals0[k] > 1e-10)
                            {
                                printf ("ERROR: %d %d %d %d: %le %le\n",
                                    M, N, P, Q, integrals0[k], 0.0);
                                errcount++;
                            }
                        }
                    
                    }
                    else if (nints == nints0 && nints != 0)
                    {
                     
                        for (k = 0; k < nints0; k++)
                        {
                            if (fabs(integrals0[k]) < 1e-6 ||
                                fabs(integrals[k]) < 1e-6)
                            {
                                if (fabs(integrals0[k] - integrals[k]) > 1e-10)
                                {
                                    printf ("* ERROR: %d %d %d %d: %le %le\n",
                                        M, N, P, Q, integrals0[k], integrals[k]);
                                    errcount++;
                                }
                            }
                            else
                            {
                                if (fabs(integrals0[k] - integrals[k])/fabs(integrals0[k]) >
                                    1e-6 && errcount < 10)
                                {
                                    printf ("* ERROR: %d %d %d %d: %le %le: %le\n",
                                        M, N, P, Q, integrals0[k], integrals[k],
                                        fabs(integrals0[k] - integrals[k])/fabs(integrals0[k]));
                                    errcount++;
                                }

                            }
                        }   
                    }
                    else
                    {
                        printf ("ERROR: nints0 %d nints %d\n", nints0, nints);
                    }

                    if (errcount > 10)
                    {
                        goto end;
                    }
                } /* for (j = start2; j < end2; j++) */
            } /* for (P = 0; P < ns; P++) */
        } /* for (i = start1; i < end1; i++) */
    } /* for (M = 0; M < ns; M++) */
    
    printf ("Done\n");
    printf ("\n");
    printf ("Number of calls: %.6le, Number of integrals: %.6le\n",
            totalcalls, totalintls);
    printf ("Total time: %.4lf secs\n", timepass);
    printf ("Average time per call: %.3le us\n",
            1000.0 * 1000.0 * timepass / totalcalls);

end:
    free (integrals0);
    for (i = 0; i < num_devices; i++)
    {
        #pragma offload target(mic:i) nocopy(erd_mic)
        {
            CInt_destroyERD (erd_mic);
        }
    }    
    CInt_destroyBasisSet (basis);
    free (shellptr);
    free (shellid);
    free (shellvalue);
    free (shellrid);
    
    fclose (ref_data_file);
    return 0;
}
