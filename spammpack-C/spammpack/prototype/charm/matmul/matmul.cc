#include "config.h"
#include "matmul.decl.h"
#include "logger.h"
#include "messages.h"
#include "timer.h"
#include "multiply.h"
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

      int c;
      const char *short_options = "hN:b:i:t:m:v";
      const option long_options[] = {
        { "help",       no_argument,        NULL, 'h' },
        { "N",          required_argument,  NULL, 'N' },
        { "block",      required_argument,  NULL, 'b' },
        { "iterations", required_argument,  NULL, 'i' },
        { "tolerance",  required_argument,  NULL, 't' },
        { "type",       required_argument,  NULL, 'm' },
        { "verify",     no_argument,        NULL, 'v' },
        { NULL, 0, NULL, 0 }
      };

      while((c = getopt_long(msg->argc, msg->argv, short_options,
              long_options, NULL)) != -1)
      {
        switch(c)
        {
          case 'h':
            CkPrintf("Usage:\n");
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

          default:
            CkExit();
            break;
        }
      }

      DEBUG("calling run() on this proxy\n");
      thisProxy.run(N, blocksize, numberIterations, tolerance, matrixType,
          verify, verifyTolerance);
    }

    void run (int N, int blocksize, int numberIterations, double tolerance,
        int matrixType, bool verify, double verifyTolerance)
    {
      CProxy_Matrix A = CProxy_Matrix::ckNew(N, blocksize);
      CProxy_Matrix C = CProxy_Matrix::ckNew(N, blocksize);

      switch(matrixType)
      {
        case full:
          DEBUG("generating random matrix\n");
          A.random(CkCallbackResumeThread());
          break;

        case decay:
          DEBUG("generating matrix with decay\n");
          A.decay(CkCallbackResumeThread());
          break;

        default:
          ABORT("unknown matrix type\n");
          break;
      }

      DEBUG("setting C to zero\n");
      C.zero(CkCallbackResumeThread());

#ifdef DEBUG_OUTPUT
      A.print(CkCallbackResumeThread());
      C.print(CkCallbackResumeThread());
#endif

      CkPrintf("running %d iterations\n", numberIterations);
      CProxy_Multiply M = CProxy_Multiply::ckNew();
      for(int iteration = 0; iteration < numberIterations; iteration++)
      {
        Timer t("iteration %d on %d PEs, multiplying C = A*A, tolerance = %e",
            iteration+1, CkNumPes(), tolerance);
        M.multiply(tolerance, A, A, C, CkCallbackResumeThread());
        t.stop();
        CkPrintf(t.to_str());
#ifdef DEBUG_OUTPUT
        C.printLeafPes(CkCallbackResumeThread());
#endif
      }

      M.getComplexity(CkCallbackResumeThread());

#ifdef DEBUG_OUTPUT
      C.print(CkCallbackResumeThread());
#endif

      if(verify)
      {
        CkPrintf("verifying result...\n");

        DenseMatrixMsg *ADense = A.getDense();
        DenseMatrixMsg *CDense = C.getDense();

        int maxDiffRow;
        int maxDiffColumn;
        double maxAbsDiff = 0;

        double *CExact = new double[N*N];
        memset(CExact, 0, sizeof(double)*N*N);

        for(int iteration = 0; iteration < numberIterations; iteration++) {
          for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
              for(int k = 0; k < N; k++)
              {
                CExact[BLOCK_INDEX(i, j, 0, 0, N)] += ADense->A[BLOCK_INDEX(i, k, 0, 0, N)]
                  *ADense->A[BLOCK_INDEX(k, j, 0, 0, N)];
              }
            }
          }
        }

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

        delete CExact;
        delete ADense;
        delete CDense;
      }

      DEBUG("done\n");
      CkExit();
    }
};

#include "matmul.def.h"
