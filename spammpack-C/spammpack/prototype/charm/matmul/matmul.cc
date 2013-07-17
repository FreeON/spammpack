#include "matmul.decl.h"

#include "config.h"
#include "logger.h"
#include "messages.h"

#include <getopt.h>

void matmulInit()
{
  LOG("on PE %d, turning manual LB on\n", CkMyPe());
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

      CkPrintf("running on %d PEs\n", CkNumPes());

      LOG("calling run on this proxy\n");
      thisProxy.run(N, blocksize);
    }

    void run (int N, int blocksize)
    {
      CProxy_Matrix A = CProxy_Matrix::ckNew(N, blocksize);
      CProxy_Matrix C = CProxy_Matrix::ckNew(N, blocksize);

      LOG("generating random matrix\n");
      A.random(CkCallbackResumeThread());

      LOG("setting C to zero\n");
      C.zero(CkCallbackResumeThread());

#ifdef DEBUG
      CkPrintf("A\n");
      A.print(CkCallbackResumeThread());
#endif

      LOG("multiplying C = A*A\n");
      C.multiply(A, A, CkCallbackResumeThread());

#ifdef DEBUG
      CkPrintf("C\n");
      C.print(CkCallbackResumeThread());
#endif

      LOG("done\n");
      CkExit();
    }
};

#include "matmul.def.h"
