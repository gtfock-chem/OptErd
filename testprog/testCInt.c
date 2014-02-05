#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#include <screening.h>

#if (defined(OPTERD_TEST_REFERENCE) + defined(OPTERD_TEST_OPTIMIZED)) != 1
#error Either OPTERD_TEST_REFERENCE or OPTERD_TEST_OPTIMIZED must be defined
#endif


int main (int argc, char **argv)
{
    BasisSet_t basis;
    ERD_t erd;
    int M;
    int N;
    int P;
    int Q;
    int nints;
    int ns;
    int nnz;
    double totalcalls = 0;
    double *integrals;
    int *shellptr;
    int *shellid;
    int *shellrid;
    double *shellvalue;
    struct timeval tv1;
    struct timeval tv2;
    double timepass;   
    double totalintls = 0;

    if (argc != 3)
    {
        printf ("Usage: %s <basisset> <xyz>\n", argv[0]);
        return 0;
    }
#if defined(OPTERD_TEST_REFERENCE)
    FILE *ref_data_file = fopen ("ivalues.ref", "w+");
#elif defined(OPTERD_TEST_OPTIMIZED)
    FILE *ref_data_file = fopen ("ivalues.ref", "r");
    if (ref_data_file == NULL)
    {
        fprintf (stderr, "ivalues.ref does not exist\n");
        exit (0);
    }
#endif
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

    // compute intergrals
    CInt_createERD (basis, &erd);
    //printf ("max memory footprint per thread = %lf KB\n",
    //    CInt_getMaxMemory (erd)/1024.0);
    
    // Free basis data used in initialization but no longer required.
    CInt_freeInitDataBasisSet (basis);

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

#if defined(OPTERD_TEST_OPTIMIZED)
    int dimMax = CInt_getMaxShellDim (basis);
    int nints0;
    int k;
    double * integrals0 = (double *)malloc (sizeof(double) * dimMax * dimMax * dimMax * dimMax);
    assert (integrals0 != NULL);
#endif    
    
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
                    CInt_computeShellQuartet (basis, erd, M, N, P, Q,
                                              &integrals, &nints);

                    gettimeofday (&tv2, NULL);
                    timepass += (tv2.tv_sec - tv1.tv_sec) +
                        (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;
                    totalintls = totalintls + nints;
                    
#if defined(OPTERD_TEST_REFERENCE)
                    fwrite (&nints, sizeof(int), 1, ref_data_file);
                    if (nints != 0)
                    {
                        fwrite (integrals, sizeof(double), nints, ref_data_file);
                    }
#endif

#if defined(OPTERD_TEST_OPTIMIZED)
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
#endif /* defined(OPTERD_TEST_OPTIMIZED) */

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
#if defined(OPTERD_TEST_OPTIMIZED)
    free (integrals0);
#endif
    CInt_destroyERD (erd);
    CInt_destroyBasisSet (basis);
    free (shellptr);
    free (shellid);
    free (shellvalue);
    free (shellrid);
    
    fclose (ref_data_file);
    return 0;
}
