#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "CInt.h"

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
    double totalcalls = 0;
    double *integrals;
    int startM;
    int startN;
    int startP;
    int startQ;
    int dimM;
    int dimN;
    int dimP;
    int dimQ;
    int iM;
    int iN;
    int iP;
    int iQ;
    double ivalue;
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
    char line[1024];
    int A;
    int B;
    int C;
    int D;
    double ivalue0;
    int errcount; 
    errcount = 0;
    FILE *ref_data_file = fopen ("ivalues.ref", "r");
    if (ref_data_file == NULL)
    {
        fprintf (stderr, "ivalues.ref does not exist\n");
    }
#endif
    // load basis set   
    CInt_createBasisSet (&basis);
    CInt_loadBasisSet (basis, argv[1], argv[2]);

    printf ("Molecule info:\n");
    printf ("  #Atoms\t= %d\n", CInt_getNumAtoms (basis));
    printf ("  #Shells\t= %d\n", CInt_getNumShells (basis));
    printf ("  #Funcs\t= %d\n", CInt_getNumFuncs (basis));
    printf ("  #OccOrb\t= %d\n", CInt_getNumOccOrb (basis));

    // compute intergrals
    CInt_createERD (basis, &erd);

    printf ("Computing integrals ...\n");
    ns = CInt_getNumShells (basis);
    timepass = 0.0;
    for (M = 0; M < ns; M++)
    {
        for (N = 0; N < ns; N++)
        {
            for (P = 0; P < ns; P++)
            {
                for (Q = 0; Q < ns; Q++)
                {
                    if (M > N || P > Q || (M + N) > (P + Q))
                        continue;
                    totalcalls = totalcalls + 1;
                    gettimeofday (&tv1, NULL);
                    CInt_computeShellQuartet (basis, erd, M, N, P, Q,
                                              &integrals, &nints);

                    gettimeofday (&tv2, NULL);
                    timepass += (tv2.tv_sec - tv1.tv_sec) +
                        (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;
                    totalintls = totalintls + nints;
                    if (nints != 0)
                    {
                        startM = CInt_getFuncStartInd (basis, M);
                        startN = CInt_getFuncStartInd (basis, N);
                        startP = CInt_getFuncStartInd (basis, P);
                        startQ = CInt_getFuncStartInd (basis, Q);
                        dimM = CInt_getShellDim (basis, M);
                        dimN = CInt_getShellDim (basis, N);
                        dimP = CInt_getShellDim (basis, P);
                        dimQ = CInt_getShellDim (basis, Q);
                        for (iM = 0; iM < dimM; iM++)
                        {
                            for (iN = 0; iN < dimN; iN++)
                            {
                                for (iP = 0; iP < dimP; iP++)
                                {
                                    for (iQ = 0; iQ < dimQ; iQ++)
                                    {
                                        ivalue =
                                            integrals[iM +
                                                      dimM * (iN +
                                                              dimN * (iP +
                                                                      dimP *
                                                                      iQ))];
#if defined(OPTERD_TEST_REFERENCE)
                                        fprintf (ref_data_file,
                                                 "%d %d %d %d %le\n",
                                                 startM + iM, startN + iN,
                                                 startP + iP, startQ + iQ,
                                                 ivalue);
#endif
#if defined(OPTERD_TEST_OPTIMIZED)
                                      hehe:
                                        assert (fgets
                                                (line, 1000,
                                                 ref_data_file) != NULL);
                                        sscanf (line, "%d %d %d %d %le\n", &A,
                                                &B, &C, &D, &ivalue0);
                                        if ((A != startM + iM)
                                            || (B != startN + iN)
                                            || (C != startP + iP)
                                            || (D != startQ + iQ))
                                        {
                                            if (ivalue0 < 1e-10)
                                                goto hehe;
                                        }
                                        assert (A == startM + iM);
                                        assert (B == startN + iN);
                                        assert (C == startP + iP);
                                        assert (D == startQ + iQ);
                                        if (fabs(ivalue) < 1e-6 || fabs(ivalue0) < 1e-6)
                                        {
                                            if (fabs (ivalue0 - ivalue) >
                                                1e-10 && errcount < 10)
                                            {
                                                printf
                                                    ("** ERROR: %d %d %d %d %le %le: %le\n",
                                                     A, B, C, D, ivalue0,
                                                     ivalue,
                                                     fabs (ivalue0 - ivalue));
                                                errcount++;
                                            }
                                        }
                                        else
                                        {
                                            if (fabs (ivalue0 - ivalue)/fabs (ivalue) >
                                                1e-6 && errcount < 10)
                                            {
                                                printf
                                                    ("ERROR: %d %d %d %d %le %le: %le\n",
                                                     A, B, C, D, ivalue0,
                                                     ivalue,
                                                     fabs (ivalue0 -
                                                           ivalue) /
                                                     fabs (ivalue0));
                                                errcount++;
                                            }

                                        }
#endif /* defined(OPTERD_TEST_OPTIMIZED) */
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printf ("Done\n");
    printf ("\n");
    printf ("Number of calls: %.6le, Number of integrals: %.6le\n",
            totalcalls, totalintls);
    printf ("Total time: %.4lf secs\n", timepass);
    printf ("Average time per call: %.3le us\n",
            1000.0 * 1000.0 * timepass / totalcalls);

    CInt_destroyERD (erd);
    CInt_destroyBasisSet (basis);
    fclose (ref_data_file);
    return 0;
}
