#include "reductions.decl.h"
#include "messages.h"
#include <getopt.h>

class Main : public CBase_Main
{
  public:

    Main (CkArgMsg *msg)
    {
      int depth = 0;
      int childsize = 2;

      int c;
      const char *short_options = "hd:c:";
      const option long_options[] = {
        { "help", no_argument, NULL, 0 },
        { "depth", required_argument, NULL, 'd' },
        { "child", required_argument, NULL, 'c' }
      };

      while((c = getopt_long(msg->argc, msg->argv, short_options,
              long_options, NULL)) != -1)
      {
        switch(c)
        {
          case 'h':
            CkPrintf("Usage:\n");
            CkPrintf("{ -h | --help }       This help\n");
            CkPrintf("{ -d | --depth } d    Create d tiers\n");
            CkPrintf("{ -c | --child } c    Create c x c children nodes on each tier\n");
            CkExit();
            break;

          case 'd':
            depth = strtol(optarg, NULL, 10);
            break;

          case 'c':
            childsize = strtol(optarg, NULL, 10);
            break;

          default:
            CkExit();
            break;
        }
      }

      CkPrintf("creating matrix of depth %d with %dx%d node grid\n", depth,
          childsize, childsize);

      CProxy_Matrix A = CProxy_Matrix::ckNew(depth, childsize);

      CkPrintf("getting norm\n");
      CkCallback cb = CkCallback(CkIndex_Main::normDone(NULL), thisProxy);
      A.norm(cb);
    }

    void normDone (DoubleMsg *norm)
    {
      CkPrintf("norm = %f\n", norm->x);
      CkExit();
    }
};

#include "reductions.def.h"
