#include "spammmultiply.decl.h"
#include "spamm.decl.h"

#include <getopt.h>

class SpAMMMultiply : public CBase_SpAMMMultiply
{
  public:

    SpAMMMultiply (CkArgMsg *msg)
    {
      //SpAMMMatrix A;

      int N = 1;

      int c;
      char shortoptions[] = "hN:";
      option longoptions[] = {
        { "help", no_argument, NULL, 'h' },
        { "N", required_argument, NULL, 'N' },
        { NULL, 0, NULL, 0 }
      };

      while((c = getopt_long(msg->argc, msg->argv, shortoptions, longoptions, NULL)) != -1)
      {
        switch(c)
        {
          case 'h':
            CkPrintf("Usage:\n");
            CkPrintf("\n");
            CkPrintf("{ -h | --help }     This help\n");
            CkPrintf("-N                  The matrix size\n");
            CkExit();
            break;

          case 'N':
            N = strtol(optarg, NULL, 10);
            break;

          default:
            CkPrintf("unknown option\n");
            CkExit();
            break;
        }
      }

      /* Create random matrix of size N. */
      //A = CProxy_SpAMMMatrix::ckNew(N);
    }

    SpAMMMultiply (CkMigrateMessage *msg) {}
};

#include "spammmultiply.def.h"
