#include "charmtest_02.decl.h"

#include <getopt.h>

class SpAMMNode : public CBase_SpAMMNode
{
  private:

    int ID;

  public:

    SpAMMNode ()
    {
      CkPrintf("initializing node\n");
    }
};

class SpAMM : public CBase_SpAMM
{
  private:

    int N;
    int N_padded;
    SpAMMNode child[4];

  public:

    SpAMM (int N)
    {
      CkPrintf("initializing SpAMM object, N = %i\n", N);
    }
};

class Main : public CBase_Main
{
  private:

    CProxy_SpAMM A;

  public:

    Main (CkArgMsg *msg)
    {
      CkPrintf("running...\n");

      int c;
      const char shortOptions[] = "hN:";
      const struct option longOptions[] = {
        { "help", no_argument,        NULL, 'h' },
        { "N",    required_argument,  NULL, 'N' },
        { NULL, 0, NULL, 0 }
      };

      int N;

      while((c = getopt_long(msg->argc, msg->argv, shortOptions, longOptions, NULL)) != -1)
      {
        switch(c)
        {
          case 'h':
            CkPrintf("Usage:\n");
            CkExit();
            break;

          case 'N':
            N = strtol(optarg, NULL, 10);
            break;

          default:
            CkPrintf("Unknown option\n");
            CkExit();
            break;
        }
      }

      CkPrintf("creating new matrix\n");
      A = CProxy_SpAMM::ckNew(N);
      CkPrintf("created new matrix\n");

      CkExit();
    }
};

#include "charmtest_02.def.h"
