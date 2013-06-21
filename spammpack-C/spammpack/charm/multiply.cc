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

      thisProxy.run(N, chunksize);
    }

    void run (int N, int chunksize)
    {
      EmptyMsg *m;
      CProxy_Matrix A;
      CProxy_Matrix C;

      /* Done reading command line. */
      LOG_INFO("starting...\n");

      /* Create new matrix. */
      A = CProxy_Matrix::ckNew(N, chunksize);
      C = CProxy_Matrix::ckNew();

      /* Create random matrix of size N. */
      for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++)
        {
          float aij = rand()/(float) RAND_MAX;
          LOG_INFO("setting matrix element A(%d,%d) <- %f\n", i, j, aij);
          m = A.set(i, j, aij); delete m;
        }
      }
      LOG_INFO("done setting matrix\n");

      /* Print matrix for debugging. */
      LOG_INFO("printing A\n");
      m = A.print(); delete m;

      /* Form matrix square. */
      m = thisProxy.multiply(A, A, C); delete m;

      /* Print result. */
      LOG_INFO("printing C <- A*A\n");
      m = C.print(); delete m;

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
    EmptyMsg * multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C)
    {
      EmptyMsg *m;

      IntMsg *N_A = A.getN();
      IntMsg *N_B = B.getN();

      if(N_A->i != N_B->i)
      {
        LOG_ERROR("mismatch in matrix dimensions\n");
        return new EmptyMsg;
      }

      IntMsg *chunksize_A = A.getChunksize();

      /* Delete C matrix. */
      LOG_INFO("Setting C to zero\n");
      m = C.remove(); delete m;

      /* Initialize C. */
      LOG_INFO("initializing C with N = %d, and chunksize = %d\n", N_A->i, chunksize_A->i);
      m = C.initialize(N_A->i, chunksize_A->i); delete m;

      /* Multiply matrices. */
      m = C.multiply(A, B); delete m;

      return new EmptyMsg;
    }
};

#include "multiply.def.h"
