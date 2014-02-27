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
 * The main program is documented as SpAMM_Charm::SpAMM_Charm. A more detailed
 * description of the API and components of this code can be found on the
 * @link API API @endlink page.
 *
 * @section example Example Use
 *
 * @code
 * ./configure.LB
 * make
 * ./charmrun +p3 spamm_charm -N 1024 -b 16 --type decay --decay 8 --tolerance 1e-8 --verify --iterations 10
 * ./charmrun +p3 spamm_charm --operation SP2 --Ne 100 --fockian something.OrthoF --tolerance 1e-8 --verify
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

#include "spamm_charm.h"

#include "backtrace.h"
#include "bcsr.h"
#include "index.h"
#include "lapack_interface.h"
#include "logger.h"
#include "messages.h"
#include "timer.h"
#include "types.h"
#include "utilities.h"

#include <assert.h>
#include <getopt.h>
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _OPENMP
#include <omp.h>
#endif

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
 * - { -h | --help }              This help
 * - { -N | --N } N               Create NxN matrix (default: 1)
 * - { -b | --block } B           Create BxB dense blocks at leaf nodes (default: 1)
 * - { -s | --basic } S           The basic blocksize (default: 1)
 * - { -i | --iterations } N      Iterate on the product N times (default:  1)
 * - { -t | --tolerance } T       Multiply with tolerance T (default: 0.00e+00)
 * - { -m | --type } TYPE         Use matrices of TYPE { full, decay, diagonal } (default: full)
 * - { -F | --fockian } FILE      Load the Fockian matrix from FILE
 * - { -P | --density } FILE      Load the density matrix from FILE
 * - { -T | --Ne } NE             The total number of electrons
 * - { -v | --verify }            Verify result
 * - { -d | --decay} GAMMA        Set matrix element decay, exp(-|i-j|/GAMMA)
 * - { -I | --intial-PE } PE      Put all chares initially on PE
 * - { -l | --load-balance }      Load balance after each iteration
 * - { -a | --align-PEs }         Align PEs for diagonal case
 * - { -p | --print-PEMap } FILE  Print a PE map in each iteration
 * - { -o | --operation } OP      Test OP { multiply, add, trace, scale, addIdentity, SP2 }
 * - { -1 | --F-min } MIN         The lower bound on the eigenspectrum of F\n");
 * - { -2 | --F-max } MAX         The upper bound on the eigenspectrum of F\n");
 *
 * @param msg The command line argument list.
 */
SpAMM_Charm::SpAMM_Charm (CkArgMsg *msg)
{
  int N = 1;
  int blocksize = 1;
  int N_basic = 1;
  int numberIterations = 0;
  double tolerance = 0.0;
  enum matrix_t matrixType = full;
  bool verify = false;
  bool loadBalance = false;
  int initialPE = CK_PE_ANY;
  bool alignPEs = false;
  bool printPEMap = false;
  char *filenamePEMap = strdup("");
  double verifyTolerance = 1.0e-10;
  double decayConstant = 0.1;
  enum operation_t operation = multiply;
  int Ne = -1;
  char *fockianFilename = NULL;
  char *densityFilename = NULL;
  double F_min = 0;
  double F_max = 0;

  int c;
  const char *short_options = "hN:b:s:i:t:m:F:P:T:vd:I:lap:o:1:2:";
  const struct option long_options[] = {
    { "help",         no_argument,        NULL, 'h' },
    { "N",            required_argument,  NULL, 'N' },
    { "block",        required_argument,  NULL, 'b' },
    { "basic",        required_argument,  NULL, 's' },
    { "iterations",   required_argument,  NULL, 'i' },
    { "tolerance",    required_argument,  NULL, 't' },
    { "type",         required_argument,  NULL, 'm' },
    { "fockian",      required_argument,  NULL, 'F' },
    { "density",      required_argument,  NULL, 'P' },
    { "Ne",           required_argument,  NULL, 'T' },
    { "verify",       no_argument,        NULL, 'v' },
    { "decay",        required_argument,  NULL, 'd' },
    { "initial-PE",   required_argument,  NULL, 'I' },
    { "load-balance", no_argument,        NULL, 'l' },
    { "align-PEs",    no_argument,        NULL, 'a' },
    { "print-PEMap",  required_argument,  NULL, 'p' },
    { "operation",    required_argument,  NULL, 'o' },
    { "F-min",        required_argument,  NULL, '1' },
    { "F-max",        required_argument,  NULL, '2' },
    { NULL, 0, NULL, 0 }
  };

  while((c = getopt_long(msg->argc, msg->argv, short_options, long_options,
          NULL)) != -1)
  {
    switch(c)
    {
      case 'h':
        CkPrintf("\n");
        CkPrintf("Usage of spamm_charm version %s\n", PACKAGE_VERSION);
        CkPrintf("\n");
        CkPrintf("{ -h | --help }               This help\n");
        CkPrintf("{ -N | --N } N                Create NxN matrix (default: %d)\n", N);
        CkPrintf("{ -b | --block } B            Create BxB dense blocks at leaf "
            "nodes (default: %d)\n", blocksize);
        CkPrintf("{ -s | --basic } S            Create SxS basic blocks (default: %d)\n", N_basic);
        CkPrintf("{ -i | --iterations } N       Iterate on the product N times (default: "
            " %d)\n", numberIterations);
        CkPrintf("{ -t | --tolerance } T        Multiply with tolerance T (default: %1.2e)\n",
            tolerance);
        CkPrintf("{ -m | --type } TYPE          Use matrices of TYPE { full, decay, diagonal } "
            "(default: full)\n");
        CkPrintf("{ -F | --fockian } FILE       Load the fockian matrix from FILE\n");
        CkPrintf("{ -P | --density } FILE       Load the density matrix from FILE\n");
        CkPrintf("{ -T | --Ne } NE              The total number of electrons\n");
        CkPrintf("{ -v | --verify }             Verify result\n");
        CkPrintf("{ -d | --decay} GAMMA         Set matrix element decay, exp(-|i-j|/GAMMA)\n");
        CkPrintf("{ -I | --intial-PE } PE       Put all chares initially on PE\n");
        CkPrintf("{ -l | --load-balance }       Load balance after each iteration\n");
        CkPrintf("{ -a | --align-PEs }          Align PEs for diagonal case\n");
        CkPrintf("{ -p | --print-PEMap } FILE   Print a PE map in each iteration\n");
        CkPrintf("{ -o | --operation } OP       Test OP { multiply, add, trace, addIdentity, scale, SP2 }\n");
        CkPrintf("{ -1 | --F-min } MIN          The lower bound on the eigenspectrum of F\n");
        CkPrintf("{ -2 | --F-max } MAX          The upper bound on the eigenspectrum of F\n");
        CkPrintf("\n");
        CkExit();
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'b':
        blocksize = strtol(optarg, NULL, 10);
        break;

      case 's':
        N_basic = strtol(optarg, NULL, 10);
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

      case 'F':
        fockianFilename = strdup(optarg);
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
        free(filenamePEMap);
        filenamePEMap = strdup(optarg);
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

        else if(strcasecmp(optarg, "addIdentity") == 0)
        {
          operation = addIdentity;
        }

        else if(strcasecmp(optarg, "scale") == 0)
        {
          operation = scale;
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

      case '1':
        F_min = strtod(optarg, NULL);
        break;

      case '2':
        F_max = strtod(optarg, NULL);
        break;

      default:
        CkExit();
        break;
    }
  }

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp master
    CkPrintf("SpAMM version %s using %d OpenMP thread(s)\n", PACKAGE_VERSION, omp_get_num_threads());
  }
#else
  CkPrintf("SpAMM version %s (serial)\n", PACKAGE_VERSION);
#endif

  /* Register backtrace handler. */
  struct sigaction act;
  memset(&act, 0, sizeof(struct sigaction));
  act.sa_handler = spamm_backtrace_handler;

  sigaction(SIGABRT, &act, NULL);
  sigaction(SIGSEGV, &act, NULL);

  switch(operation)
  {
    case SP2:
      {
        if(Ne < 0)
        {
          ABORT("missing number of electrons, --Ne\n");
        }

        if(fockianFilename == NULL)
        {
          INFO("missing fockian file, will generate random one\n");
          fockianFilename = strdup("");
        }

        if(densityFilename == NULL)
        {
          INFO("missing density file, will not compare result\n");
          densityFilename = strdup("");
        }

        if(numberIterations <= 0)
        {
          numberIterations = 100;
        }

        DEBUG("calling runSP2() on this proxy\n");
        thisProxy.runSP2(strlen(fockianFilename), fockianFilename,
            strlen(densityFilename), densityFilename, Ne, N, blocksize,
            N_basic, numberIterations, tolerance, loadBalance, initialPE,
            alignPEs, printPEMap, strlen(filenamePEMap), filenamePEMap, F_min,
            F_max);
      }
      break;

    default:

      if(numberIterations <= 0)
      {
        numberIterations = 1;
      }

      DEBUG("calling run() on this proxy\n");
      thisProxy.run(N, blocksize, N_basic, numberIterations, tolerance,
          matrixType, decayConstant, operation, verify, verifyTolerance,
          loadBalance, initialPE, alignPEs, printPEMap);
      break;
  }
}

/** The main method.
 *
 * @param N The size of the matrix.
 * @param blocksize The blocksize of the matrix.
 * @param N_basic The size of the basic sub-matrices.
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
void SpAMM_Charm::run (int N, int blocksize, int N_basic,
    int numberIterations, double tolerance, int matrixType,
    double decayConstant, int operation, bool verify, double verifyTolerance,
    bool loadBalance, int initialPE, bool alignPEs, bool printPEMap)
{
  double alpha = 0.8;
  double beta = 1.12;

  LBDatabase *db = LBDatabaseObj();

  CProxy_Matrix A = CProxy_Matrix::ckNew(initialPE, alignPEs, N, blocksize, N_basic, strlen("A"), "A");
  CProxy_Matrix C = CProxy_Matrix::ckNew(initialPE, alignPEs, N, blocksize, N_basic, strlen("C"), "C");

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

  MatrixNodeMsg *ANodes = A.getNodes();
  MatrixNodeMsg *CNodes = C.getNodes();

  CProxy_Multiply M;

  if(operation == multiply)
  {
    M = CProxy_Multiply::ckNew(A, A, C);
    M.init(initialPE, alignPEs, CkCallbackResumeThread());
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
          t.start();
          M.multiply(tolerance, 1.0, 1.0, CkCallbackResumeThread());
          t.stop();
          CkPrintf("%s\n", t.to_str());
        }
        break;

      case add:
        {
          Timer t("iteration %d on %d PEs, adding C = A+B", iteration+1, CkNumPes());
          t.start();
          C.add(0.0, 1.0, A, CkCallbackResumeThread());
          C.add(1.0, 1.0, A, CkCallbackResumeThread());
          t.stop();
          CkPrintf("%s\n", t.to_str());
        }
        break;

      case trace:
        {
          Timer t("iteration %d on %d PEs, trace(A)", iteration+1, CkNumPes());
          t.start();
          A.updateTrace(CkCallbackResumeThread());
          t.stop();
          CkPrintf("%s\n", t.to_str());
        }
        break;

      case scale:
        {
          Timer t("iteration %d on %d PEs, scale A by %e", iteration+1, CkNumPes(), alpha);
          t.start();
          C.setEqual(A, CkCallbackResumeThread());
          C.scale(alpha, CkCallbackResumeThread());
          t.stop();
          CkPrintf("%s\n", t.to_str());
        }
        break;

      case addIdentity:
        {
          Timer t("iteration %d on %d PEs, alpha * A + beta * I (alpha = %e, beta = %e)",
              iteration+1, CkNumPes(), alpha, beta);
          t.start();
          C.setEqual(A, CkCallbackResumeThread());
#ifdef PRINT_MATRICES
          DenseMatrixMsg *msg = C.toDense();
          printDense(N, msg->A, "C (after setEqual)");
          delete msg;
#endif
          C.addIdentity(alpha, beta, CkCallbackResumeThread());
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
      t.start();
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
          }
          break;

        case add:
          for(int i = 0; i < N*N; i++)
          {
            CExact[i] = 2*ADense[i];
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

        case scale:
          for(int i = 0; i < N*N; i++)
          {
            CExact[i] = alpha*ADense[i];
          }
          break;

        case addIdentity:
          for(int i = 0; i < N*N; i++)
          {
            CExact[i] = alpha*ADense[i];
          }
          for(int i = 0; i < N; i++)
          {
            CExact[BLOCK_INDEX(i, i, 0, 0, N)] += beta;
          }
          break;

        default:
          ABORT("unknown operation\n");
          break;
      }
      t.stop();
      CkPrintf("%s\n", t.to_str());
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

    double trace = 0;

    for(int i = 0; i < N; i++)
    {
      trace += CExact[BLOCK_INDEX(i, i, 0, 0, N)];

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

    C.updateTrace(CkCallbackResumeThread());
    DoubleMsg *sTrace = C.getTrace();

    INFO("trace(CExact) = % e\n", trace);
    INFO("trace(C)      = % e\n", sTrace->x);

    if(fabs(sTrace->x - trace) > verifyTolerance)
    {
      ABORT("trace mismatch (abs. tolerance = %e), "
          "trace(C) = %e, trace(CExact) = %e, "
          "abs. diff = %e\n", verifyTolerance, sTrace->x, trace,
          fabs(sTrace->x - trace));
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
 * @param lengthFockianFilename The length of the filename string.
 * @param fockianFilename The filename of the Fockian matrix file.
 * @param lengthDensityFilename The length of the filename string.
 * @param densityFilename The filename of the density matrix file.
 * @param Ne The total number of electrons.
 * @param N The matrix size in case the filename is empty.
 * @param blocksize The blocksize of the matrix.
 * @param N_basic The size of the basic sub-matrices.
 * @param maxIterations The maximum number of iterations.
 * @param tolerance The SpAMM tolerance.
 * @param loadBalance Whether to load balance.
 * @param initialPE The initial PE to put chares on.
 * @param alignPEs Align PEs in the diagonal matrix case.
 * @param printPEMap Whether to print a PE map in each iteration.
 * @param lengthPEMap The length of filenamePEMap.
 * @param filenamePEMap The base name of the PEMap files.
 * @param F_min The lower bound on the eigenspectrum.
 * @param F_max The upper bound on the eigenspectrum.
 */
void SpAMM_Charm::runSP2 (int lengthFockianFilename, char *fockianFilename,
    int lengthDensityFilename, char *densityFilename, int Ne, int N,
    int blocksize, int N_basic, int maxIterations, double tolerance,
    bool loadBalance, int initialPE, bool alignPEs, bool printPEMap,
    int lengthPEMap, char *filenamePEMap, double F_min, double F_max)
{
  double *PDense;
  int NRows, NColumns;

  LBDatabase *db = LBDatabaseObj();

  if(lengthFockianFilename > 0)
  {
    BCSR F(fockianFilename);

    if(F_min >= F_max)
    {
      F.getSpectralBounds(0, &F_min, &F_max);
    }

    /* Convert BCSR to dense. */
    F.toDense(&NRows, &NColumns, &PDense);
    assert(NRows == NColumns);
  }

  else
  {
    INFO("creating random, symmetric %d x %d matrix\n", N, N);
    NRows = N;
    NColumns = N;
    PDense = new double[N*N];
    for(int i = 0; i < N; i++) {
      for(int j = i; j < N; j++)
      {
        PDense[BLOCK_INDEX(i, j, 0, 0, N)] = rand()/(double) RAND_MAX;
        PDense[BLOCK_INDEX(j, i, 0, 0, N)] = PDense[BLOCK_INDEX(i, j, 0, 0, N)];
      }
    }
    getSpectralBounds(0, &F_min, &F_max, N, PDense);
  }

  INFO("spectral bounds: [ % e, % e ], dF = %e\n", F_min, F_max, F_max-F_min);

#ifdef PRINT_MATRICES
  printDense(NRows, PDense, "F");
#endif

  CProxy_Matrix P = CProxy_Matrix::ckNew(initialPE, alignPEs, NRows,
      blocksize, N_basic, strlen("P"), (char*) "P");

  P.set(NRows, PDense, CkCallbackResumeThread());
  //P.init(CkCallbackResumeThread());

  Timer total_time("total time");
  total_time.start();

  /* Scale Fockian to get initial density matrix guess.
   *
   * P_0 = (F_max*I-F)/(F_max-F_min)
   */
  P.addIdentity(-1, F_max, CkCallbackResumeThread());
  P.scale(1/(F_max-F_min), CkCallbackResumeThread());

  delete[] PDense;

#ifdef PRINT_MATRICES
  DenseMatrixMsg *P0Dense = P.toDense();
  printDense(P0Dense->N, P0Dense->A, "P0");
  delete P0Dense;
#endif

  CProxy_Matrix P2 = CProxy_Matrix::ckNew(initialPE, alignPEs, NRows,
      blocksize, N_basic, strlen("P2"), (char*) "P2");
  //P2.init(CkCallbackResumeThread());

  MatrixInfoMsg *PInfo = P.info();
  MatrixNodeMsg *PNodes = P.getNodes();
  MatrixNodeMsg *P2Nodes = P2.getNodes();

  CProxy_Multiply M = CProxy_Multiply::ckNew(P, P, P2);
  M.init(initialPE, alignPEs, CkCallbackResumeThread());

  int full_complexity = (int) ceil(NRows/(double) blocksize);
  full_complexity = full_complexity*full_complexity*full_complexity;

  delete PInfo;
  delete PNodes;
  delete P2Nodes;

  /* Start SP2 iterations. */
  double occupation[4] = { 0, 0, 0, 0 };
  P.updateTrace(CkCallbackResumeThread());
  DoubleMsg *trace_P = P.getTrace();
  INFO("tolerance = %e\n", tolerance);
  INFO("iteration  0: trace(P) = %e (Ne/2 = %e)\n", trace_P->x, Ne/2.0);
  bool converged = false;
  for(int iteration = 0; iteration < maxIterations; iteration++)
  {
    Timer tMultiply("multiply");
    Timer tSetEqual("setEq");
    Timer tAdd("add");

    tMultiply.start();
    M.multiply(tolerance, 1.0, 0.0, CkCallbackResumeThread()); /* P2 <- P*P */
    tMultiply.stop();

    M.updateComplexity(CkCallbackResumeThread());
    IntMsg *complexity = M.getComplexity();

    P2.updateTrace(CkCallbackResumeThread());
    DoubleMsg *trace_P2 = P2.getTrace();
    DEBUG("trace(P%d^2) = %e (Ne/2 = %e)\n", iteration, trace_P2->x, Ne/2.0);

    if(fabs(trace_P2->x-Ne/2.0) < fabs(2*trace_P->x-trace_P2->x-Ne/2.0))
    {
      DEBUG("P%d <- P%d^2\n", iteration+1, iteration);
      tSetEqual.start();
      P.setEqual(P2, CkCallbackResumeThread());
      tSetEqual.stop();
      delete trace_P;
      trace_P = trace_P2;
    }

    else
    {
      DEBUG("P%d <- 2P%d - P%d^2\n", iteration+1, iteration, iteration);
      tAdd.start();
      P.add(2, -1, P2, CkCallbackResumeThread());
      tAdd.stop();
      trace_P->x = 2*trace_P->x-trace_P2->x;
      delete trace_P2;
    }

    INFO("iteration %2d, %s, %s, %s: trace(P) = %e "
        "(Ne = %d, 2*trace(P)-Ne = % e) "
        "complexity %d (out of %d)\n",
        iteration+1,
        tMultiply.to_str(), tSetEqual.to_str(), tAdd.to_str(),
        trace_P->x, Ne, 2*trace_P->x-Ne, complexity->i, full_complexity);
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
      char filename[2000];

      snprintf(filename, 2000, "%s.%03d.PEMap", filenamePEMap, iteration+1);

      int fd = open(filename, O_WRONLY | O_EXCL | O_CREAT, 00644);
      if(fd == -1)
      {
        if(errno == EEXIST)
        {
          ABORT("PEMap file %s already exists\n", filename);
        }
        ABORT("error opening PEMap file %s: %s\n", filename, strerror(errno));
      }

      FILE *fs = fdopen(fd, "w");

      MatrixInfoMsg *PInfo = P.info();
      int NTier = 1 << PInfo->depth;

      DEBUG("NTier = %d\n", NTier);

      P.updatePEMap(CkCallbackResumeThread());
      PEMapMsg *PEMap = P.getPEMap();

      fprintf(fs, "PEMap for matrix P:\n");
      for(int i = 0; i < NTier; i++) {
        for(int j = 0; j < NTier; j++)
        {
          int matrix_offset = BLOCK_INDEX(i, j, 0, 0, NTier);
          fprintf(fs, "PEMap(%d,%d) = %d (norm = %e)\n", i, j,
              PEMap->PEMap[matrix_offset],
              PEMap->PEMap_norm[matrix_offset]);
        }
      }
      fprintf(fs, "end of PEMap for matrix P\n");
      delete PEMap;

      P2.updatePEMap(CkCallbackResumeThread());
      PEMap = P2.getPEMap();

      fprintf(fs, "PEMap for matrix P2:\n");
      for(int i = 0; i < NTier; i++) {
        for(int j = 0; j < NTier; j++)
        {
          int matrix_offset = BLOCK_INDEX(i, j, 0, 0, NTier);
          fprintf(fs, "PEMap(%d,%d) = %d (norm = %e)\n", i, j,
              PEMap->PEMap[matrix_offset],
              PEMap->PEMap_norm[matrix_offset]);
        }
      }
      fprintf(fs, "end of PEMap for matrix P2\n");
      delete PEMap;

      M.updatePEMap(CkCallbackResumeThread());
      PEMap = M.getPEMap();

      fprintf(fs, "PEMap for convolution:\n");
      for(int i = 0; i < NTier; i++) {
        for(int j = 0; j < NTier; j++) {
          for(int k = 0; k < NTier; k++)
          {
            int matrix_offset = BLOCK_INDEX_3(i, j, k, NTier);
            fprintf(fs, "PEMap(%d,%d,%d) = %d (norm = %e)\n", i, j, k,
                PEMap->PEMap[matrix_offset],
                PEMap->PEMap_norm[matrix_offset]);
          }
        }
      }
      fprintf(fs, "end of PEMap for convolution\n");
      delete PEMap;

      fclose(fs);
    }

    /* Load balance. */
    if(loadBalance)
    {
      DEBUG("load balancing\n");
      db->StartLB();
      CkWaitQD();
    }
  }

  INFO("idempotency error          = %e\n", fabs(occupation[0]-occupation[1]));
  INFO("previous idempotency error = %e\n", fabs(occupation[2]-occupation[3]));

  if(!converged)
  {
    ABORT("SP2 did not converge in %d steps\n", maxIterations);
  }

  //DenseMatrixMsg *PFinal = P.toDense();
  //printDense(NRows, PFinal->A, "PFinal");
  //delete PFinal;

  total_time.stop();
  INFO("%s\n", total_time.to_str());

  INFO("done\n");
  CkExit();
}

#include "spamm_charm.def.h"
