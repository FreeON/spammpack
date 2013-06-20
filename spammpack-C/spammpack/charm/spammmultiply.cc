#include "spammmultiply.decl.h"

#include <getopt.h>

class SpAMMMultiply : public CBase_SpAMMMultiply
{
  public:

    SpAMMMultiply (CkArgMsg *msg)
    {
      run(msg->argc, msg->argv);
    }

    void run (int argc, char **argv)
    {
      CProxy_SpAMMMatrix A;
      CProxy_SpAMMMatrix C;

      int N = 1;
      int chunksize = 1;

      int c;
      char shortoptions[] = "hN:c:";
      option longoptions[] = {
        { "help",   no_argument,        NULL, 'h' },
        { "N",      required_argument,  NULL, 'N' },
        { "chunk",  required_argument,  NULL, 'c' },
        { NULL, 0, NULL, 0 }
      };

      while((c = getopt_long(argc, argv, shortoptions, longoptions, NULL)) != -1)
      {
        switch(c)
        {
          case 'h':
            CkPrintf("Usage:\n");
            CkPrintf("\n");
            CkPrintf("{ -h | --help }     This help\n");
            CkPrintf("-N N                Set the matrix size to NxN\n");
            CkPrintf("{ -c | --chunk } N  Set the chunk size to NxN\n");
            CkExit();
            break;

          case 'N':
            N = strtol(optarg, NULL, 10);
            break;

          case 'c':
            chunksize = strtol(optarg, NULL, 10);
            break;

          default:
            CkPrintf("unknown option\n");
            CkExit();
            break;
        }
      }

      /* Done reading command line. */
      CkPrintf("starting...\n");

      /* Create new matrix. */
      A = CProxy_SpAMMMatrix::ckNew(N, chunksize);

      /* Create random matrix of size N. */
      for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++)
        {
          A.set(i, j, rand()/(float) RAND_MAX);
        }
      }

      /* Print matrix for debugging. */
      A.print();

      /* Form matrix square. */
      multiply(A, A, C);

      /* Done. */
      //CkExit();
    }

    SpAMMMultiply (CkMigrateMessage *msg) {}

    /** Multiply two matrices and store the produce in a third one, @f$ C
     * \leftarrow A \times B @f$. The result matrix C is overwritten.
     *
     * @param A The matrix A.
     * @param B The matrix B.
     * @param C The matrix C.
     */
    void multiply (CProxy_SpAMMMatrix A,
        CProxy_SpAMMMatrix B,
        CProxy_SpAMMMatrix C)
    {
    }
};

#include "spammmultiply.def.h"
