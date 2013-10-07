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
 * ./charmrun +p3 spamm --operation SP2 --Ne 100 --density something.OrthoF --tolerance 1e-8
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
#include "lapack_interface.h"
#include "messages.h"
#include "timer.h"
#include "logger.h"
#include "types.h"
#include "index.h"
#include "bcsr.h"
#include "utilities.h"

#include <assert.h>
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
 * - { -P | --density } FILE   Load the density matrix from FILE
 * - { -T | --Ne } NE          The total number of electrons
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
  int Ne = -1;
  char *densityFilename = NULL;

  int c;
  const char *short_options = "hN:b:i:t:m:P:T:vd:I:lapo:";
  const struct option long_options[] = {
    { "help",         no_argument,        NULL, 'h' },
    { "N",            required_argument,  NULL, 'N' },
    { "block",        required_argument,  NULL, 'b' },
    { "iterations",   required_argument,  NULL, 'i' },
    { "tolerance",    required_argument,  NULL, 't' },
    { "type",         required_argument,  NULL, 'm' },
    { "density",      required_argument,  NULL, 'P' },
    { "Ne",           required_argument,  NULL, 'T' },
    { "verify",       no_argument,        NULL, 'v' },
    { "decay",        required_argument,  NULL, 'd' },
    { "initial-PE",   required_argument,  NULL, 'I' },
    { "load-balance", no_argument,        NULL, 'l' },
    { "align-PEs",    no_argument,        NULL, 'a' },
    { "print-PEMap",  no_argument,        NULL, 'p' },
    { "operation",    required_argument,  NULL, 'o' },
    { NULL, 0, NULL, 0 }
  };

  while((c = getopt_long(msg->argc, msg->argv, short_options, long_options,
          NULL)) != -1)
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
        CkPrintf("{ -P | --density } FILE   Load the density matrix from FILE\n");
        CkPrintf("{ -T | --Ne } NE          The total number of electrons\n");
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

      case 'P':
        densityFilename = strdup(optarg);
        break;

      case 'T':
        Ne = strtol(optarg, NULL, 10);
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
      {
        if(Ne < 0)
        {
          ABORT("missing number of electrons, --Ne\n");
        }

        if(densityFilename == NULL)
        {
          ABORT("missing density file\n");
        }

        DEBUG("calling runSP2() on this proxy\n");
        thisProxy.runSP2(strlen(densityFilename), densityFilename, Ne,
            blocksize, numberIterations, tolerance, loadBalance, initialPE,
            alignPEs, printPEMap);
      }
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

  CProxy_Matrix A = CProxy_Matrix::ckNew(initialPE, alignPEs, N, blocksize, strlen("A"), "A");
  CProxy_Matrix C = CProxy_Matrix::ckNew(initialPE, alignPEs, N, blocksize, strlen("B"), "C");

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
  printDense(N, ADense, "ADense");
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
          M.multiply(tolerance, 1.0, 1.0, CkCallbackResumeThread());
          t.stop();
          CkPrintf("%s\n", t.to_str());
        }
        break;

      case add:
        {
          Timer t("iteration %d on %d PEs, adding C = A+B", iteration+1, CkNumPes());
          C.add(0.0, 1.0, A, CkCallbackResumeThread());
          C.add(1.0, 1.0, A, CkCallbackResumeThread());
          t.stop();
          CkPrintf("%s\n", t.to_str());
        }
        break;

      case trace:
        {
          Timer t("iteration %d on %d PEs, trace(A)", iteration+1, CkNumPes());
          A.updateTrace(CkCallbackResumeThread());
          t.stop();
          CkPrintf("%s\n", t.to_str());
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
            CkPrintf("%s\n", t.to_str());
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

    int maxDiffRow = -1;
    int maxDiffColumn = -1;
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
 * @param length The length of the filename string.
 * @param filename The filename of the densit file.
 * @param Ne The total number of electrons.
 * @param blocksize The blocksize of the matrix.
 * @param maxIterations The maximum number of iterations.
 * @param tolerance The SpAMM tolerance.
 * @param loadBalance Whether to load balance.
 * @param initialPE The initial PE to put chares on.
 * @param alignPEs Align PEs in the diagonal matrix case.
 * @param printPEMap Whether to print a PE map in each iteration.
 */
void SpAMM::runSP2 (int length, char *filename, int Ne, int blocksize,
    int maxIterations, double tolerance, bool loadBalance, int initialPE,
    bool alignPEs, bool printPEMap)
{
  Timer *t;
  double *PDense;
  int NRows, NColumns;
  double F_min, F_max;

  LBDatabase *db = LBDatabaseObj();

  if(maxIterations == 1)
  {
    maxIterations = 100;
  }

  BCSR F(filename);

  F.getSpectralBounds(0, &F_min, &F_max);
  INFO("spectral bounds: [ % e, % e ], dF = %e\n", F_min, F_max, F_max-F_min);

  F.toDense(&NRows, &NColumns, &PDense);

  assert(NRows == NColumns);

  //printDense(NRows, PDense, "F");

  CProxy_Matrix P = CProxy_Matrix::ckNew(initialPE, alignPEs, NRows,
      blocksize, strlen("P"), (char*) "P");

  P.set(NRows, PDense, CkCallbackResumeThread());

  Timer total_time("total time");

  /* Scale Fockian to get initial density matrix guess.
   *
   * P_0 = (F_max*I-F)/(F_max-F_min)
   */
  P.addIdentity(-1, F_max, CkCallbackResumeThread());
  P.scale(1/(F_max-F_min), CkCallbackResumeThread());

  delete[] PDense;

  //DenseMatrixMsg *P0Dense = P.toDense();
  //printDense(P0Dense->N, P0Dense->A, "P0");
  //delete P0Dense;

  CProxy_Matrix P2 = CProxy_Matrix::ckNew(initialPE, alignPEs, NRows,
      blocksize, strlen("P2"), (char*) "P2");

  MatrixInfoMsg *PInfo = P.info();
  MatrixNodeMsg *PNodes = P.getNodes(PInfo->depth);
  MatrixNodeMsg *P2Nodes = P2.getNodes(PInfo->depth);

  CProxy_Multiply M = CProxy_Multiply::ckNew(initialPE, alignPEs, P, P, P2,
      blocksize, PInfo->depth, PNodes->nodes, PNodes->nodes, P2Nodes->nodes);

  delete PInfo;
  delete PNodes;
  delete P2Nodes;

  /* Start SP2 iterations. */
  double occupation[4] = { 0, 0, 0, 0 };
  P.updateTrace(CkCallbackResumeThread());
  DoubleMsg *trace_P = P.getTrace();
  INFO("iteration  0: trace(P) = %e (Ne/2 = %e)\n", trace_P->x, Ne/2.0);
  bool converged = false;
  for(int iteration = 0; iteration < maxIterations; iteration++)
  {
    t = new Timer("iteration %2d", iteration+1);

    M.multiply(tolerance, 1.0, 0.0, CkCallbackResumeThread()); /* P2 <- P*P */

    M.updateComplexity(CkCallbackResumeThread());
    IntMsg *complexity = M.getComplexity();

    //P0Dense = P.toDense();
    //printDense(P0Dense->N, P0Dense->A, "P%d^2", iteration+1);
    //delete P0Dense;

    P2.updateTrace(CkCallbackResumeThread());
    DoubleMsg *trace_P2 = P2.getTrace();
    DEBUG("trace(P%d^2) = %e (Ne/2 = %e)\n", iteration, trace_P2->x, Ne/2.0);

    if(fabs(trace_P2->x-Ne/2.0) < fabs(2*trace_P->x-trace_P2->x-Ne/2.0))
    {
      DEBUG("P%d <- P%d^2\n", iteration+1, iteration);
      P.setEqual(P2, CkCallbackResumeThread());
      delete trace_P;
      trace_P = trace_P2;
    }

    else
    {
      DEBUG("P%d <- 2P%d - P%d^2\n", iteration+1, iteration, iteration);
      P.add(2, -1, P2, CkCallbackResumeThread());
      trace_P->x = 2*trace_P->x-trace_P2->x;
      delete trace_P2;
    }

    t->stop();
    INFO("%s: trace(P) = %e (Ne/2 = %e, trace(P)-Ne/2 = % e) complexity %d\n", t->to_str(),
        trace_P->x, Ne/2.0, trace_P->x-Ne/2.0, complexity->i);
    delete t;
    delete complexity;

    occupation[0] = trace_P->x;
    for(int i = 3; i >= 1; i--)
    {
      occupation[i] = occupation[i-1];
    }

    if(iteration > MIN(20, maxIterations))
    {
      double idempotencyErrorNow = fabs(occupation[1]-occupation[2]);
      if(idempotencyErrorNow < 1.0e-2)
      {
        double idempotencyErrorLast = fabs(occupation[3]-occupation[2]);
        if(idempotencyErrorNow >= idempotencyErrorLast)
        {
          INFO("SP2 converged in %d steps\n", iteration+1);
          converged = true;
          break;
        }
      }
    }

    /* Print the PE map. */
    if(printPEMap)
    {
      MatrixInfoMsg *PInfo = P.info();
      int NTier = 1 << PInfo->depth;

      DEBUG("NTier = %d\n", NTier);

      P.updatePEMap(CkCallbackResumeThread());
      PEMapMsg *PEMap = P.getPEMap();

      CkPrintf("PEMap for matrix P:\n");
      for(int i = 0; i < NTier; i++) {
        for(int j = 0; j < NTier; j++)
        {
          int matrix_offset = BLOCK_INDEX(i, j, 0, 0, NTier);
          CkPrintf("PEMap(%d,%d) = %d (norm = %e)\n", i, j,
              PEMap->PEMap[matrix_offset],
              PEMap->PEMap_norm[matrix_offset]);
        }
      }
      CkPrintf("end of PEMap for matrix P\n");
      delete PEMap;

      P2.updatePEMap(CkCallbackResumeThread());
      PEMap = P2.getPEMap();

      CkPrintf("PEMap for matrix P2:\n");
      for(int i = 0; i < NTier; i++) {
        for(int j = 0; j < NTier; j++)
        {
          int matrix_offset = BLOCK_INDEX(i, j, 0, 0, NTier);
          CkPrintf("PEMap(%d,%d) = %d (norm = %e)\n", i, j,
              PEMap->PEMap[matrix_offset],
              PEMap->PEMap_norm[matrix_offset]);
        }
      }
      CkPrintf("end of PEMap for matrix P2\n");
      delete PEMap;

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

    /* Load balance. */
    if(loadBalance)
    {
      DEBUG("load balancing\n");
      db->StartLB();
      CkWaitQD();
    }
  }

  if(!converged)
  {
    INFO("SP2 did not converge in %d steps\n", maxIterations);
  }

  INFO("idempotency error          = %e\n", fabs(occupation[0]-occupation[1]));
  INFO("previous idempotency error = %e\n", fabs(occupation[2]-occupation[3]));

  //DenseMatrixMsg *PFinal = P.toDense();
  //printDense(NRows, PFinal->A, "PFinal");
  //delete PFinal;

  total_time.stop();
  INFO("%s\n", total_time.to_str());

  INFO("done\n");
  CkExit();
}

#include "spamm.def.h"
