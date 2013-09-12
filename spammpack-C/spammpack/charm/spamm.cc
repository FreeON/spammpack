/** @file
 *
 * Multiply two matrices with threshold.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 *
 *
 * @mainpage
 *
 * @section Introduction
 *
 * <a href="http://charm.cs.uiuc.edu/">Charm++</a> prototype implementation of
 * the SpAMM algorithm
 * @cite ChallacombeBock2010
 * @cite BockChallacombe2012
 * @cite BockSISC2013.
 * The SP2 algorithm is described in more detail in
 * @cite Niklasson2003.
 * The main program is documented as SpAMM::SpAMM. A more detailed description
 * of the API and components of this code can be found on the @link API API
 * @endlink page.
 *
 * @section example Example Use
 *
 * @code
 * ./configure.LB
 * make
 * ./charmrun +p3 spamm -N 1024 -b 16 --type decay --decay 8 --tolerance 1e-8 --verify --iterations 10
 * @endcode
 *
 * @section References
 *
 * A list of @link citelist references @endlink can be found here.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "spamm.h"
#include "blas_interface.h"
#include "messages.h"
#include "timer.h"
#include "logger.h"
#include "types.h"
#include "index.h"

#include <getopt.h>

/** Initialize the node. */
void initialize (void)
{
  /* Set the load balancing mode to manual. */
  TurnManualLBOn();

  /* Initialize the logger. */
  initializeLogger();
}

/** The main method.
 *
 * Currently known command line arguments are:
 *
 * - { -h | --help }           This help
 * - { -N | --N } N            Create NxN matrix (default: 1)
 * - { -b | --block } B        Create BxB dense blocks at leaf nodes (default: 1)
 * - { -i | --iterations } N   Iterate on the product N times (default:  1)
 * - { -t | --tolerance } T    Multiply with tolerance T (default: 0.00e+00)
 * - { -m | --type } TYPE      Use matrices of TYPE { full, decay, diagonal } (default: full)
 * - { -v | --verify }         Verify result
 * - { -d | --decay} GAMMA     Set matrix element decay, exp(-|i-j|/GAMMA)
 * - { -I | --intial-PE } PE   Put all chares initially on PE
 * - { -l | --load-balance }   Load balance after each iteration
 * - { -a | --align-PEs }      Align PEs for diagonal case
 * - { -p | --print-PEMap }    Print a PE map in each iteration
 * - { -o | --operation } OP   Test OP { multiply, add, trace, SP2 }
 *
 * @param msg The command line argument list.
 */
SpAMM::SpAMM (CkArgMsg *msg)
{
  int N = 1;
  int blocksize = 1;
  int numberIterations = 1;
  double tolerance = 0.0;
  enum matrix_t matrixType = full;
  bool verify = false;
  bool loadBalance = false;
  int initialPE = CK_PE_ANY;
  bool alignPEs = false;
  bool printPEMap = false;
  double verifyTolerance = 1.0e-10;
  double decayConstant = 0.1;
  enum operation_t operation = multiply;
  int Ne;
  double *PDense;

  int c;
  const char *short_options = "hN:b:i:t:m:vd:I:lapo:";
  const option long_options[] = {
    { "help",         no_argument,        NULL, 'h' },
    { "N",            required_argument,  NULL, 'N' },
    { "block",        required_argument,  NULL, 'b' },
    { "iterations",   required_argument,  NULL, 'i' },
    { "tolerance",    required_argument,  NULL, 't' },
    { "type",         required_argument,  NULL, 'm' },
    { "verify",       no_argument,        NULL, 'v' },
    { "decay",        required_argument,  NULL, 'd' },
    { "initial-PE",   required_argument,  NULL, 'I' },
    { "load-balance", no_argument,        NULL, 'l' },
    { "align-PEs",    no_argument,        NULL, 'a' },
    { "print-PEMap",  no_argument,        NULL, 'p' },
    { "operation",    required_argument,  NULL, 'o' },
    { NULL, 0, NULL, 0 }
  };

  while((c = getopt_long(msg->argc, msg->argv, short_options,
          long_options, NULL)) != -1)
  {
    switch(c)
    {
      case 'h':
        CkPrintf("\n");
        CkPrintf("Usage of spamm version %s\n", PACKAGE_VERSION);
        CkPrintf("\n");
        CkPrintf("{ -h | --help }           This help\n");
        CkPrintf("{ -N | --N } N            Create NxN matrix (default: %d)\n", N);
        CkPrintf("{ -b | --block } B        Create BxB dense blocks at leaf "
            "nodes (default: %d)\n", blocksize);
        CkPrintf("{ -i | --iterations } N   Iterate on the product N times (default: "
            " %d)\n", numberIterations);
        CkPrintf("{ -t | --tolerance } T    Multiply with tolerance T (default: %1.2e)\n",
            tolerance);
        CkPrintf("{ -m | --type } TYPE      Use matrices of TYPE { full, decay, diagonal } "
            "(default: full)\n");
        CkPrintf("{ -v | --verify }         Verify result\n");
        CkPrintf("{ -d | --decay} GAMMA     Set matrix element decay, exp(-|i-j|/GAMMA)\n");
        CkPrintf("{ -I | --intial-PE } PE   Put all chares initially on PE\n");
        CkPrintf("{ -l | --load-balance }   Load balance after each iteration\n");
        CkPrintf("{ -a | --align-PEs }      Align PEs for diagonal case\n");
        CkPrintf("{ -p | --print-PEMap }    Print a PE map in each iteration\n");
        CkPrintf("{ -o | --operation } OP   Test OP { multiply, add, trace, SP2 }\n");
        CkPrintf("\n");
        CkExit();
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'b':
        blocksize = strtol(optarg, NULL, 10);
        break;

      case 'i':
        numberIterations = strtol(optarg, NULL, 10);
        break;

      case 't':
        tolerance = strtod(optarg, NULL);
        break;

      case 'm':
        if(strcasecmp(optarg, "full") == 0)
        {
          matrixType = full;
        }
        else if (strcasecmp(optarg, "decay") == 0)
        {
          matrixType = decay;
        }
        else if (strcasecmp(optarg, "diagonal") == 0)
        {
          matrixType = diagonal;
        }
        else
        {
          ABORT("unknown matrix type\n");
        }
        break;

      case 'v':
        verify = true;
        break;

      case 'd':
        decayConstant = strtod(optarg, NULL);
        break;

      case 'I':
        initialPE = strtol(optarg, NULL, 10);
        break;

      case 'l':
        loadBalance = true;
        break;

      case 'a':
        alignPEs = true;
        break;

      case 'p':
        printPEMap = true;
        break;

      case 'o':
        if(strcasecmp(optarg, "multiply") == 0)
        {
          operation = multiply;
        }

        else if(strcasecmp(optarg, "add") == 0)
        {
          operation = add;
        }

        else if(strcasecmp(optarg, "trace") == 0)
        {
          operation = trace;
        }

        else if(strcasecmp(optarg, "SP2") == 0)
        {
          operation = SP2;
        }

        else
        {
          ABORT("unknown operation\n");
        }
        break;

      default:
        CkExit();
        break;
    }
  }

  CkPrintf("SpAMM version %s\n", PACKAGE_VERSION);

  switch(operation)
  {
    case SP2:
      DEBUG("calling runSP2() on this proxy\n");
      thisProxy.runSP2(N, PDense, Ne);
      break;

    default:
      DEBUG("calling run() on this proxy\n");
      thisProxy.run(N, blocksize, numberIterations, tolerance, matrixType,
          decayConstant, operation, verify, verifyTolerance, loadBalance,
          initialPE, alignPEs, printPEMap);
      break;
  }
}

/** The main method.
 *
 * @param N The size of the matrix.
 * @param blocksize The blocksize of the matrix.
 * @param numberIterations Total number of iterations for multiply.
 * @param tolerance The SpAMM tolerance.
 * @param matrixType The matrix type.
 * @param decayConstant The decay constant for matrices with exponential
 * decay.
 * @param operation The operation of type operation_t.
 * @param verify Whether to verify the matrix product.
 * @param verifyTolerance The absolute tolerance in the matrix
 * verification.
 * @param loadBalance Whether to load balance.
 * @param initialPE The initial PE to put chares on.
 * @param alignPEs Align PEs in the diagonal matrix case.
 * @param printPEMap Whether to print a PE map in each iteration.
 */
void SpAMM::run (int N, int blocksize, int numberIterations, double tolerance,
    int matrixType, double decayConstant, int operation, bool verify,
    double verifyTolerance, bool loadBalance, int initialPE, bool alignPEs,
    bool printPEMap)
{
  LBDatabase *db = LBDatabaseObj();

  CProxy_Matrix A = CProxy_Matrix::ckNew(initialPE, alignPEs, N, blocksize, 2, "A");
  CProxy_Matrix C = CProxy_Matrix::ckNew(initialPE, alignPEs, N, blocksize, 2, "C");

  /* Initialize the matrices. */
  double *ADense = new double[N*N];
  memset(ADense, 0, sizeof(double)*N*N);

  switch(matrixType)
  {
    case full:
      for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++)
        {
          ADense[BLOCK_INDEX(i, j, 0, 0, N)] = rand()/(double) RAND_MAX;
        }
      }
      break;

    case diagonal:
      for(int i = 0; i < N; i++)
      {
        ADense[BLOCK_INDEX(i, i, 0, 0, N)] = 0.3*(rand()/(double) RAND_MAX - 0.5)+1;
      }
      break;

    case decay:
      for(int i = 0; i < N; i++)
      {
        ADense[BLOCK_INDEX(i, i, 0, 0, N)] = 0.3*(rand()/(double) RAND_MAX - 0.5)+1;
      }

      for(int i = 0; i < N; i++) {
        for(int j = i+1; j < N; j++)
        {
          ADense[BLOCK_INDEX(i, j, 0, 0, N)] = exp(-fabs(i-j)/decayConstant)
            *ADense[BLOCK_INDEX(i, i, 0, 0, N)];
          ADense[BLOCK_INDEX(j, i, 0, 0, N)] = ADense[BLOCK_INDEX(i, j, 0, 0, N)];
        }
      }
      break;

    default:
      ABORT("unknown matrix type\n");
      break;
  }

#ifdef PRINT_MATRICES
  printDense(N, ADense, "ADense:");
#endif

  A.set(N, ADense, CkCallbackResumeThread());

  double *CExact = NULL;
  if(verify)
  {
    CExact = new double[N*N];
    memset(CExact, 0, sizeof(double)*N*N);
  }

  MatrixInfoMsg *AInfo = A.info();
  MatrixInfoMsg *CInfo = C.info();

  MatrixNodeMsg *ANodes = A.getNodes(AInfo->depth);
  MatrixNodeMsg *CNodes = C.getNodes(CInfo->depth);

  CProxy_Multiply M;

  if(operation == multiply)
  {
    M = CProxy_Multiply::ckNew(initialPE, alignPEs, A, A, C,
        AInfo->blocksize, AInfo->depth, ANodes->nodes, ANodes->nodes,
        CNodes->nodes);
  }

  delete ANodes;
  delete CNodes;

  CkPrintf("running %d iterations\n", numberIterations);
  for(int iteration = 0; iteration < numberIterations; iteration++)
  {
    switch(operation)
    {
      case multiply:
        {
          Timer t("iteration %d on %d PEs, multiplying C = A*A, tolerance = %e",
              iteration+1, CkNumPes(), tolerance);
          M.multiply(tolerance, CkCallbackResumeThread());
          t.stop();
          CkPrintf(t.to_str());
        }
        break;

      case add:
        {
          Timer t("iteration %d on %d PEs, adding C = A+B", iteration+1, CkNumPes());
          C.add(0.0, 1.0, A, CkCallbackResumeThread());
          C.add(1.0, 1.0, A, CkCallbackResumeThread());
          t.stop();
          CkPrintf(t.to_str());
        }
        break;

      case trace:
        {
          Timer t("iteration %d on %d PEs, trace(A)", iteration+1, CkNumPes());
          A.updateTrace(CkCallbackResumeThread());
          t.stop();
          CkPrintf(t.to_str());
        }
        break;

      default:
        ABORT("unknow operation\n");
        break;
    }

    if(printPEMap)
    {
      int NTier = 1 << AInfo->depth;

      DEBUG("NTier = %d\n", NTier);

      CkPrintf("PE map for A\n");
      A.updatePEMap(CkCallbackResumeThread());
      PEMapMsg *PEMap = A.getPEMap();

      CkPrintf("PEMap for matrix A:\n");
      for(int i = 0; i < NTier; i++) {
        for(int j = 0; j < NTier; j++)
        {
          int matrix_offset = BLOCK_INDEX(i, j, 0, 0, NTier);
          CkPrintf("PEMap(%d,%d) = %d (norm = %e)\n", i, j,
              PEMap->PEMap[matrix_offset],
              PEMap->PEMap_norm[matrix_offset]);
        }
      }
      CkPrintf("end of PEMap for matrix A\n");
      delete PEMap;

      CkPrintf("PE map for C\n");
      C.updatePEMap(CkCallbackResumeThread());
      PEMap = C.getPEMap();

      CkPrintf("PEMap for matrix C:\n");
      for(int i = 0; i < NTier; i++) {
        for(int j = 0; j < NTier; j++)
        {
          int matrix_offset = BLOCK_INDEX(i, j, 0, 0, NTier);
          CkPrintf("PEMap(%d,%d) = %d (norm = %e)\n", i, j,
              PEMap->PEMap[matrix_offset],
              PEMap->PEMap_norm[matrix_offset]);
        }
      }
      CkPrintf("end of PEMap for matrix C\n");
      delete PEMap;

      if(operation == multiply)
      {
        CkPrintf("PE map for convolution\n");
        M.updatePEMap(CkCallbackResumeThread());
        PEMap = M.getPEMap();

        CkPrintf("PEMap for convolution:\n");
        for(int i = 0; i < NTier; i++) {
          for(int j = 0; j < NTier; j++) {
            for(int k = 0; k < NTier; k++)
            {
              int matrix_offset = BLOCK_INDEX_3(i, j, k, NTier);
              CkPrintf("PEMap(%d,%d,%d) = %d (norm = %e)\n", i, j, k,
                  PEMap->PEMap[matrix_offset],
                  PEMap->PEMap_norm[matrix_offset]);
            }
          }
        }
        CkPrintf("end of PEMap for convolution\n");
        delete PEMap;
      }
    }

    /* Load balance. */
    if(loadBalance)
    {
      INFO("load balancing\n");
      db->StartLB();
      CkWaitQD();
    }

    if(verify)
    {
      /* Calculate the reference matrix. */
      Timer t("calculating reference result");
      switch(operation)
      {
        case multiply:
          {
#ifdef DGEMM
            double alpha = 1;
            double beta = 1;
            DGEMM("N", "N", &N, &N, &N, &alpha, ADense, &N, ADense, &N, &beta, CExact, &N);
#else
            for(int i = 0; i < N; i++) {
              for(int j = 0; j < N; j++) {
                for(int k = 0; k < N; k++)
                {
                  CExact[BLOCK_INDEX(i, j, 0, 0, N)] += ADense[BLOCK_INDEX(i, k, 0, 0, N)]
                    *ADense[BLOCK_INDEX(k, j, 0, 0, N)];
                }
              }
            }
#endif
            t.stop();
            CkPrintf(t.to_str());
          }
          break;

        case add:
          {
            for(int i = 0; i < N*N; i++)
            {
              CExact[i] = 2*ADense[i];
            }
          }
          break;

        case trace:
          {
            double trace = 0;
            for(int i = 0; i < N; i++)
            {
              trace += ADense[BLOCK_INDEX(i, i, 0, 0, N)];
            }
            DoubleMsg *spammTrace = A.getTrace();
            double absDiff = fabs(trace-spammTrace->x);
            if(absDiff > verifyTolerance)
            {
              ABORT("trace mismatch (abs. tolerance = %e, "
                  "abs. diff = %e), %e (reference) vs. %e (SpAMM)\n",
                  verifyTolerance, absDiff, trace, spammTrace->x);
            }
            CkPrintf("trace verified\n");
          }
          break;

        default:
          ABORT("unknown operation\n");
          break;
      }
    }
  }

#ifdef PRINT_MATRICES
  {
    if(verify)
    {
      printDense(N, CExact, "CExact");
    }
    DenseMatrixMsg *CDense = C.toDense();
    printDense(N, CDense->A, "C");
    delete CDense;
  }
#endif

  if(verify && operation != trace)
  {
    CkPrintf("verifying result\n");

    DenseMatrixMsg *CDense = C.toDense();

    int maxDiffRow;
    int maxDiffColumn;
    double maxAbsDiff = 0;

    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++)
      {
        double absDiff = fabs(CExact[BLOCK_INDEX(i, j, 0, 0, N)]
            -CDense->A[BLOCK_INDEX(i, j, 0, 0, N)]);

        if(absDiff > maxAbsDiff)
        {
          maxDiffRow = i;
          maxDiffColumn = j;
          maxAbsDiff = absDiff;
        }
      }
    }

    if(maxAbsDiff > verifyTolerance)
    {
      double relDiff = (CExact[BLOCK_INDEX(maxDiffRow, maxDiffColumn, 0, 0, N)] != 0
          ? maxAbsDiff/CExact[BLOCK_INDEX(maxDiffRow, maxDiffColumn, 0, 0, N)]
          : 0);
      ABORT("result mismatch (abs. tolerance = %e, "
          "abs. diff = %e, rel. diff = %e), "
          "C(%d,%d): %e (reference) vs. %e (SpAMM)\n",
          verifyTolerance, maxAbsDiff, relDiff,
          maxDiffRow, maxDiffColumn,
          CExact[BLOCK_INDEX(maxDiffRow, maxDiffColumn, 0, 0, N)],
          CDense->A[BLOCK_INDEX(maxDiffRow, maxDiffColumn, 0, 0, N)]);
    }
    CkPrintf("result verified\n");

    delete[] ADense;
    delete[] CExact;

    delete CDense;
  }

  delete AInfo;
  delete CInfo;

  INFO("done\n");
  CkExit();
}

/** The main method for running an SP2 calculation.
 *
 * @param N The matrix size.
 * @param PDense The initial density matrix.
 * @param Ne The total number of electrons.
 */
void SpAMM::runSP2 (int N, double *PDense, int Ne)
{
}

#include "spamm.def.h"
