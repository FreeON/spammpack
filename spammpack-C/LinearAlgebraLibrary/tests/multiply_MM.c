#include "config.h"
#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  struct spamm_t A;
  struct spamm_t A2;

  double *A_dense;
  double *A2_dense;
  double alpha, beta;

  struct spamm_tree_stats_t stats;

  int i, j, k;

  spamm_log("loading matrix\n", __FILE__, __LINE__);
  spamm_read_MM("e20r0000.mtx", 10, 10, 2, 2, 1e-10, &A);
  spamm_tree_stats(&stats, &A);
  spamm_log("read %ix%i matrix, %i nodes, %i dense blocks, tree = %i bytes, blocks = %i bytes (%1.1f%%)\n",
      __FILE__, __LINE__, A.M, A.N, stats.number_nodes, stats.number_dense_blocks,
      stats.memory_tree, stats.memory_dense_blocks, (stats.memory_tree+stats.memory_dense_blocks)/(double) (A.M*A.N*sizeof(double))*100);

  spamm_new(A.M, A.N, A.M_block, A.N_block, A.M_child, A.N_child, A.threshold, &A2);

  spamm_log("converting tree to dense\n", __FILE__, __LINE__);
  spamm_spamm_to_dense(&A, &A_dense);
  A2_dense = (double*) malloc(sizeof(double)*A.M*A.N);
  for (i = 0; i < A.M; ++i) {
    for (j = 0; j < A.N; ++j)
    {
      A2_dense[spamm_dense_index(i, j, A.N)] = 0.0;
    }
  }

#ifdef DGEMM
  spamm_log("multiplying matrix with BLAS to get reference product\n", __FILE__, __LINE__);
  alpha = 1.0;
  beta = 0.0;
  DGEMM("N", "N", &A.M, &A.N, &A.M, &alpha, A_dense, &A.M, A_dense, &A.N, &beta, A2_dense, &A.M);
#else
  spamm_log("multiplying matrix directly to get reference product\n", __FILE__, __LINE__);
  for (i = 0; i < A.M; ++i) {
    for (j = 0; j < A.N; ++j) {
      for (k = 0; k < A.M; ++k)
      {
        A2_dense[spamm_dense_index(i, j, A.N)] += A_dense[spamm_dense_index(i, k, A.N)]*A_dense[spamm_dense_index(k, j, A.N)];
      }
    }
  }
#endif

  spamm_log("multiplying matrix with spamm\n", __FILE__, __LINE__);
  spamm_multiply(1.0, &A, &A, 0.0, &A2);

  spamm_log("comparing matrices\n", __FILE__, __LINE__);
  for (i = 0; i < A.M; ++i) {
    for (j = 0; j < A.N; ++j)
    {
      if (spamm_get(i, j, &A2) != A2_dense[spamm_dense_index(i, j, A.N)])
      {
        printf("[multiply_spamm] mismatch: (A2[%i][%i] = %e) != (A2_dense[%i][%i] = %e)\n", i, j, spamm_get(i, j, &A2), i, j, A2_dense[spamm_dense_index(i, j, A.N)]);
        exit(1);
      }
    }
  }

  return 0;
}
