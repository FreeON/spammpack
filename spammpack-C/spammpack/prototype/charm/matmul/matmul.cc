/** @file
 *
 * Multiply two matrices with threshold.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
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
 * The main program is documented as Main.
 *
 * @section example Example Use
 *
 * @code
 * ./configure.LB
 * make
 * ./charmrun +p3 matmul -N 1024 -b 16 --type decay --decay 8 --tolerance 1e-8 --verify --iterations 10
 * @endcode
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "matmul.h"
#include "blas_interface.h"
#include "messages.h"
#include "timer.h"
#include "logger.h"
#include "types.h"
#include "index.h"

#include <getopt.h>

/** Set the load balancing mode to manual. */
void setLBMode (void)
{
  TurnManualLBOn();
}

/** The main method.
 *
 * @param msg The command line argument list.
 */
Main::Main (CkArgMsg *msg)
{
  int N = 1;
  int blocksize = 1;
  int numberIterations = 1;
  double tolerance = 0.0;
  enum matrix_t matrixType = full;
  bool verify = false;
  double verifyTolerance = 1.0e-10;
  double decayConstant = 0.1;

  int c;
  const char *short_options = "hN:b:i:t:m:vd:";
  const option long_options[] = {
    { "help",       no_argument,        NULL, 'h' },
    { "N",          required_argument,  NULL, 'N' },
    { "block",      required_argument,  NULL, 'b' },
    { "iterations", required_argument,  NULL, 'i' },
    { "tolerance",  required_argument,  NULL, 't' },
    { "type",       required_argument,  NULL, 'm' },
    { "verify",     no_argument,        NULL, 'v' },
    { "decay",      required_argument,  NULL, 'd' },
    { NULL, 0, NULL, 0 }
  };

  while((c = getopt_long(msg->argc, msg->argv, short_options,
          long_options, NULL)) != -1)
  {
    switch(c)
    {
      case 'h':
        CkPrintf("Usage:\n");
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
        CkPrintf("{ -v | --verify }         Verify matmul product\n");
        CkPrintf("{ -d | --decay} GAMMA     Set matrix element decay, exp(-|i-j|/GAMMA)\n");
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
        verify = !verify;
        break;

      case 'd':
        decayConstant = strtod(optarg, NULL);
        break;

      default:
        CkExit();
        break;
    }
  }

  CkPrintf("matmul version %s\n", PACKAGE_VERSION);

  DEBUG("calling run() on this proxy\n");
  thisProxy.run(N, blocksize, numberIterations, tolerance, matrixType,
      decayConstant, verify, verifyTolerance);
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
 * @param verify Whether to verify the matrix product.
 * @param verifyTolerance The absolute tolerance in the matrix
 * verification.
 */
void Main::run (int N, int blocksize, int numberIterations, double tolerance,
    int matrixType, double decayConstant, bool verify, double verifyTolerance)
{
  LBDatabase *db = LBDatabaseObj();

  CProxy_Matrix A = CProxy_Matrix::ckNew(N, blocksize);
  CProxy_Matrix C = CProxy_Matrix::ckNew(N, blocksize);

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

  MatrixInfoMsg *AInfo = A.info();
  MatrixInfoMsg *CInfo = C.info();

  MatrixNodeMsg *ANodes = A.getNodes(AInfo->depth);
  MatrixNodeMsg *CNodes = C.getNodes(AInfo->depth);

  /* Set the A matrix. */
  {
    double *block = new double[blocksize*blocksize];

    for(int i = 0; i < AInfo->NPadded/blocksize; i++) {
      for(int j = 0; j < AInfo->NPadded/blocksize; j++)
      {
        memset(block, 0, sizeof(double)*blocksize*blocksize);

        for(int l = i*blocksize; l < (i+1)*blocksize && l < N; l++) {
          for(int m = j*blocksize; m < (j+1)*blocksize && m < N; m++)
          {
            block[BLOCK_INDEX(l, m, i*blocksize, j*blocksize, blocksize)] =
              ADense[BLOCK_INDEX(l, m, 0, 0, N)];
          }
        }

        ANodes->nodes(i, j).set(blocksize, block);
      }
    }

    delete[] block;
  }

  double *CExact = NULL;
  if(verify)
  {
    CExact = new double[N*N];
    memset(CExact, 0, sizeof(double)*N*N);
  }

  CProxy_Multiply M = CProxy_Multiply::ckNew(A, A, C, AInfo->blocksize,
      AInfo->depth, ANodes->nodes, ANodes->nodes, CNodes->nodes);

  delete AInfo;
  delete CInfo;

  delete ANodes;
  delete CNodes;

  CkPrintf("running %d iterations\n", numberIterations);
  for(int iteration = 0; iteration < numberIterations; iteration++)
  {
    Timer t("iteration %d on %d PEs, multiplying C = A*A, tolerance = %e",
        iteration+1, CkNumPes(), tolerance);
    M.multiply(tolerance, CkCallbackResumeThread());
    t.stop();
    CkPrintf(t.to_str());

#ifdef PRINT_PE
    INFO("PE map for A\n");
    A.printPE(CkCallbackResumeThread());

    INFO("PE map for C\n");
    C.printPE(CkCallbackResumeThread());

    INFO("PE map for convolution\n");
    M.printPE(CkCallbackResumeThread());
#endif

    /* Load balance. */
    db->StartLB();
    CkWaitQD();

    if(verify)
    {
      /* Calculate the reference matrix. */
      Timer t("calculating reference result");
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

  if(verify)
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
          "C(%d,%d): %e (reference) vs. %e (matmul)\n",
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

  DEBUG("done\n");
  CkExit();
}

#include "matmul.def.h"
