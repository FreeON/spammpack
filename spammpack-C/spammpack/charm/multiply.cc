#include "multiply.decl.h"
#include "messages.h"
#include "utilities.h"

#include <getopt.h>

class Multiply : public CBase_Multiply
{
  public:

    Multiply (CkArgMsg *msg)
    {
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

      while((c = getopt_long(msg->argc, msg->argv, shortoptions, longoptions, NULL)) != -1)
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

      run(N, chunksize);
    }

    void run (int N, int chunksize)
    {
      CProxy_Matrix A;
      CProxy_Matrix C;

      /* Done reading command line. */
      LOG_INFO("starting...\n");

      /* Create new matrix. */
      A = CProxy_Matrix::ckNew(N, chunksize);

      /* Create random matrix of size N. */
      for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++)
        {
          float aij = rand()/(float) RAND_MAX;
          LOG_INFO("setting matrix element A(%d,%d) <- %f\n", i, j, aij);
          EmptyMsg *m = A.set(i, j, aij);
          delete m;
        }
      }
      LOG_INFO("done setting matrix\n");

      /* Print matrix for debugging. */
      A.print();

      /* Form matrix square. */
      //multiply(A, A, C);

      /* Done. */
      CkExit();
    }

    /** Multiply two matrices and store the produce in a third one, @f$ C
     * \leftarrow A \times B @f$. The result matrix C is overwritten.
     *
     * @param A The matrix A.
     * @param B The matrix B.
     * @param C The matrix C.
     */
    void multiply (CProxy_Matrix A,
        CProxy_Matrix B,
        CProxy_Matrix C)
    {
    }
};

#include "multiply.def.h"
