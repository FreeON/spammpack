#include "config.h"
#include "matmul.decl.h"
#include "logger.h"
#include "messages.h"
#include "timer.h"
#include <getopt.h>

#ifdef VERIFY_MULTIPLY
#include "index.h"
#endif

void matmulInit()
{
  //DEBUG("on PE %d, turning manual LB on\n", CkMyPe());
  //TurnManualLBOn();
}

class Main : public CBase_Main
{
  public:

    Main (CkArgMsg *msg)
    {
      int N = 1;
      int blocksize = 1;

      int c;
      const char *short_options = "hN:b:";
      const option long_options[] = {
        { "help",  no_argument,       NULL, 0 },
        { "N",     required_argument, NULL, 'N' },
        { "block", required_argument, NULL, 'b' }
      };

      while((c = getopt_long(msg->argc, msg->argv, short_options,
              long_options, NULL)) != -1)
      {
        switch(c)
        {
          case 'h':
            CkPrintf("Usage:\n");
            CkPrintf("{ -h | --help }       This help\n");
            CkPrintf("{ -N | --N } N        Create NxN matrix (default: %d)\n", N);
            CkPrintf("{ -b | --block } B    Create BxB dense blocks at leaf "
                "nodes (default: %d)\n", blocksize);
            CkExit();
            break;

          case 'N':
            N = strtol(optarg, NULL, 10);
            break;

          case 'b':
            blocksize = strtol(optarg, NULL, 10);
            break;

          default:
            CkExit();
            break;
        }
      }

      DEBUG("calling run on this proxy\n");
      thisProxy.run(N, blocksize);
    }

    void run (int N, int blocksize)
    {
      CProxy_Matrix A = CProxy_Matrix::ckNew(N, blocksize);
      CProxy_Matrix C = CProxy_Matrix::ckNew(N, blocksize);

      DEBUG("generating random matrix\n");
      A.random(CkCallbackResumeThread());

#ifdef VERIFY_MULTIPLY
      DenseMatrixMsg *ADense = A.getDense();
#endif

      DEBUG("setting C to zero\n");
      C.zero(CkCallbackResumeThread());

#ifdef DEBUG_OUTPUT
      A.print(CkCallbackResumeThread());
#endif

      Timer t("%d PEs, multiplying C = A*A", CkNumPes());
      C.multiply(A, A, CkCallbackResumeThread());
      t.stop();
      CkPrintf(t.to_str());

#ifdef DEBUG_OUTPUT
      C.print(CkCallbackResumeThread());
#endif

#ifdef VERIFY_MULTIPLY
      DenseMatrixMsg *CDense = C.getDense();

      for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++)
        {
          double CExact = 0;
          for(int k = 0; k < N; k++)
          {
            CExact += ADense->A[BLOCK_INDEX(i, k, 0, 0, N)]
              *ADense->A[BLOCK_INDEX(k, j, 0, 0, N)];
          }

          if(fabs(CExact-CDense->A[BLOCK_INDEX(i, j, 0, 0, N)]) > VERIFY_TOLERANCE)
          {
            ABORT("result mismatch (abs. tolerance = %e), "
                "C(%d,%d): %e vs. %e\n", VERIFY_TOLERANCE, i, j,
                CExact, CDense->A[BLOCK_INDEX(i, j, 0, 0, N)]);
          }
        }
      }

      CkPrintf("result verified\n");

      delete ADense;
      delete CDense;
#endif

      DEBUG("done\n");
      CkExit();
    }
};

#include "matmul.def.h"
