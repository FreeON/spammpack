#include "charmtest_02.decl.h"

#include <getopt.h>

class SpAMMNode : public CBase_SpAMMNode
{
  private:

    int N_lower[2];
    int N_upper[2];

    CProxy_SpAMMNode child[4];

  public:

    SpAMMNode (const int N_lower[2],
        const int N_upper[2])
    {
      CkPrintf("initializing node\n");

      for(int i = 0; i < 2; i++)
      {
        this->N_lower[i] = N_lower[i];
        this->N_upper[i] = N_upper[i];
      }

      for(int i = 0; i < 4; i++)
      {
        child[i] = NULL;
      }
    }

    void set (const int i,
        const int j,
        const float aij)
    {
    }
};

class SpAMM : public CBase_SpAMM
{
  private:

    int N;
    int chunksize;
    int N_padded;
    CProxy_SpAMMNode root;

    void initialize ()
    {
      for(int i = 0; i < N_padded; i++) {
        for(int j = 0; j < N_padded; j++)
        {
          set(i, j, rand()/(double) RAND_MAX);
        }
      }
    }

  public:

    SpAMM (const int N,
        const int chunksize)
    {
      CkPrintf("initializing SpAMM object, N = %i, chunksize = %i\n",
          N, chunksize);

      root = NULL;

      this->N = N;
      this->chunksize = chunksize;

      int depth = int(ceil(log(N/chunksize)/log(2)));
      N_padded = int(chunksize*pow(2, depth));

      CkPrintf("N_padded = %i\n", this->N_padded);

      /* Load a random matrix. */
      initialize();
    }

    void set (const int i,
        const int j,
        const float aij)
    {
      if(root == NULL)
      {
        int N_lower[2] = { 0, N_padded };
        int N_upper[2] = { 0, N_padded };

        CProxy_SpAMMNode new_root = CProxy_SpAMMNode::ckNew(N_lower, N_upper);
        root = &new_root;
      }

      root->set(i, j, aij);
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
      const char shortOptions[] = "hN:c:";
      const struct option longOptions[] = {
        { "help",   no_argument,        NULL, 'h' },
        { "N",      required_argument,  NULL, 'N' },
        { "chunk",  required_argument,  NULL, 'c' },
        { NULL, 0, NULL, 0 }
      };

      int N = 1;
      int chunksize = 1;

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

          case 'c':
            chunksize = strtol(optarg, NULL, 10);
            break;

          default:
            CkPrintf("Unknown option\n");
            CkExit();
            break;
        }
      }

      /* Create new matrix. */
      A = CProxy_SpAMM::ckNew(N, chunksize);
    }

    void done ()
    {
      CkPrintf("done...\n");
    }
};

#include "charmtest_02.def.h"
