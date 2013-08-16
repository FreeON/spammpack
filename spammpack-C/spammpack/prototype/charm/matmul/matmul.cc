/** @file
 *
 * Multiply two matrices with threshold.
 */

#include "config.h"
#include "matmul.decl.h"
#include "messages.h"
#include "timer.h"
#include "logger.h"
#include "types.h"
#include "index.h"

#include <getopt.h>

class Main : public CBase_Main
{
  public:

    Main (CkArgMsg *msg)
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
            CkPrintf("{ -m | --type } TYPE      Use matrices of TYPE { full, decay } "
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

      DEBUG("calling run() on this proxy\n");
      thisProxy.run(N, blocksize, numberIterations, tolerance, matrixType,
          decayConstant, verify, verifyTolerance);
    }

    void run (int N, int blocksize, int numberIterations, double tolerance,
        int matrixType, double decayConstant, bool verify,
        double verifyTolerance)
    {
      CProxy_Matrix A = CProxy_Matrix::ckNew(N, blocksize);
      CProxy_Matrix C = CProxy_Matrix::ckNew(N, blocksize);

      DenseMatrixMsg *ADense = NULL;
      DenseMatrixMsg *CExact = NULL;

      if(verify)
      {
        ADense = A.toDense();
        CExact = C.toDense();
      }

      MatrixInfoMsg *AInfo = A.info();
      MatrixInfoMsg *CInfo = C.info();

      CProxy_Multiply M = CProxy_Multiply::ckNew(A, A, C, AInfo->blocksize,
          AInfo->depth, AInfo->nodes, AInfo->nodes, CInfo->nodes);

      delete AInfo;
      delete CInfo;

      CkPrintf("running %d iterations\n", numberIterations);
      for(int iteration = 0; iteration < numberIterations; iteration++)
      {
        Timer t("iteration %d on %d PEs, multiplying C = A*A, tolerance = %e",
            iteration+1, CkNumPes(), tolerance);
        M.multiply(tolerance, CkCallbackResumeThread());
        t.stop();
        CkPrintf(t.to_str());

        if(verify)
        {
          CkPrintf("calculating reference result\n");
          for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
              for(int k = 0; k < N; k++)
              {
                CExact->A[BLOCK_INDEX(i, j, 0, 0, N)] += ADense->A[BLOCK_INDEX(i, k, 0, 0, N)]
                  *ADense->A[BLOCK_INDEX(k, j, 0, 0, N)];
              }
            }
          }
        }
      }

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
            double absDiff = fabs(CExact->A[BLOCK_INDEX(i, j, 0, 0, N)]
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
          double relDiff = (CExact->A[BLOCK_INDEX(maxDiffRow, maxDiffColumn, 0, 0, N)] != 0
              ? maxAbsDiff/CExact->A[BLOCK_INDEX(maxDiffRow, maxDiffColumn, 0, 0, N)]
              : 0);
          ABORT("result mismatch (abs. tolerance = %e, "
              "abs. diff = %e, rel. diff = %e), "
              "C(%d,%d): %e (reference) vs. %e (matmul)\n",
              verifyTolerance, maxAbsDiff, relDiff,
              maxDiffRow, maxDiffColumn,
              CExact->A[BLOCK_INDEX(maxDiffRow, maxDiffColumn, 0, 0, N)],
              CDense->A[BLOCK_INDEX(maxDiffRow, maxDiffColumn, 0, 0, N)]);
        }
        CkPrintf("result verified\n");

        delete ADense;
        delete CDense;
        delete CExact;
      }

      DEBUG("done\n");
      CkExit();
    }
};

#include "matmul.def.h"
