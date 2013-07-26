#include "config.h"
#include "matmul.decl.h"
#include "logger.h"
#include "messages.h"
#include "timer.h"
#include "multiply.h"
#include <getopt.h>

#ifdef VERIFY_MULTIPLY
#include "index.h"
#endif

class Main : public CBase_Main
{
  public:

    Main (CkArgMsg *msg)
    {
      int N = 1;
      int blocksize = 1;
      int numberIterations = 1;

      int c;
      const char *short_options = "hN:b:i:";
      const option long_options[] = {
        { "help",       no_argument,        NULL, 0 },
        { "N",          required_argument,  NULL, 'N' },
        { "block",      required_argument,  NULL, 'b' },
        { "iterations", required_argument,  NULL, 'i' },
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

          default:
            CkExit();
            break;
        }
      }

      DEBUG("calling run on this proxy\n");
      thisProxy.run(N, blocksize, numberIterations);
      //thisProxy.reduceTest(N);
    }

    void run (int N, int blocksize, int numberIterations)
    {
#ifdef DIRECT_MULTIPLY
      CProxy_Node A = CProxy_Node::ckNew();
      CProxy_Node C = CProxy_Node::ckNew();

      A(0, 0).insert(N, 0, blocksize, 0, 0, N, 0, N);
      C(0, 0).insert(N, 0, blocksize, 0, 0, N, 0, N);

      A.doneInserting();
      C.doneInserting();
#else
      CProxy_Matrix A = CProxy_Matrix::ckNew(N, blocksize);
      CProxy_Matrix C = CProxy_Matrix::ckNew(N, blocksize);
#endif

      DEBUG("generating random matrix\n");
      A.random(CkCallbackResumeThread());

#ifdef VERIFY_MULTIPLY
      DenseMatrixMsg *ADense = A.getDense();
#endif

      DEBUG("setting C to zero\n");
      C.zero(CkCallbackResumeThread());

#ifdef DEBUG_OUTPUT
      A.print(CkCallbackResumeThread());
      C.print(CkCallbackResumeThread());
#endif

      CkPrintf("running %d iterations\n", numberIterations);
#ifndef DIRECT_MULTIPLY
      CProxy_Multiply M = CProxy_Multiply::ckNew();
#endif
      for(int i = 0; i < numberIterations; i++)
      {
        Timer t("iteration %d on %d PEs, multiplying C = A*A", i, CkNumPes());
#ifdef DIRECT_MULTIPLY
        int depth = 0;
        CProxy_MultiplyElement convolution = CProxy_MultiplyElement::ckNew();
        for(int i = 0; i < (1 << depth); i++) {
          for(int j = 0; j < (1 << depth); j++) {
            for(int k = 0; k < (1 << depth); k++)
            {
              DEBUG("inserting convolution at C(%d,%d) <- A(%d,%d) * B(%d,%d)\n",
                  i, j, i, k, k, j);

              convolution(i, j, k).insert(blocksize, A, A, C);
            }
          }
        }
        convolution.doneInserting();
        INFO("multiplying...\n");
        convolution.multiply(CkCallbackResumeThread());
        convolution.storeBack(CkCallbackResumeThread());
        INFO("done\n");
#else
        M.multiply(A, A, C, CkCallbackResumeThread());
#endif
        t.stop();
        CkPrintf(t.to_str());
#ifdef DEBUG_OUTPUT
        C.printLeafPes(CkCallbackResumeThread());
#endif
      }

#ifdef DEBUG_OUTPUT
      C.print(CkCallbackResumeThread());
#endif

#ifdef VERIFY_MULTIPLY
      DenseMatrixMsg *CDense = C.getDense();

      CkPrintf("verifying result...\n");

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

    void reduceTest (int N)
    {
      INFO("reduce test, N = %d\n", N);

      CProxy_ReductionData d = CProxy_ReductionData::ckNew(N, N);
      CProxy_Reduction r = CProxy_Reduction::ckNew();

      for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
          for(int k = 0; k < N; k++)
          {
            r(i, j, k).insert(N, d);
          }
        }
      }
      r.doneInserting();

      r.reduce(CkCallbackResumeThread());

      INFO("done\n");
      CkExit();
    }
};

#include "matmul.def.h"
