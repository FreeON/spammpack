#include "config.h"
#include <spamm.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define THRESHOLD 1e-14

int
main ()
{
  struct spamm_t A;
  struct spamm_t A2;

  float_t *A_dense;
  float_t *A2_dense;
  float_t alpha, beta;

  struct spamm_tree_stats_t stats;

  double max_diff;
  int max_diff_i, max_diff_j;

  int i, j, k;

  int result = 0;

  spamm_log("loading matrix\n", __FILE__, __LINE__);
  spamm_read_MM("e20r0000.mtx", 32, 32, 2, 2, 1e-10, &A);
  //spamm_read_MM("sparseHamiltonian.mtx", 32, 32, 2, 2, 1e-10, &A);
  spamm_tree_stats(&stats, &A);
  spamm_log("read %ix%i matrix, %i nodes, %i dense blocks, tree = %i bytes, blocks = %i bytes (%1.1f%%)\n",
      __FILE__, __LINE__, A.M, A.N, stats.number_nodes, stats.number_dense_blocks,
      stats.memory_tree, stats.memory_dense_blocks, (stats.memory_tree+stats.memory_dense_blocks)/(double) (A.M*A.N*sizeof(double))*100);

  spamm_log("converting tree to dense\n", __FILE__, __LINE__);
  spamm_spamm_to_dense(&A, &A_dense);
  A2_dense = (float_t*) malloc(sizeof(float_t)*A.M*A.N);
  for (i = 0; i < A.M; ++i) {
    for (j = 0; j < A.N; ++j)
    {
      A2_dense[spamm_dense_index(i, j, A.M, A.N)] = 0.0;
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
        A2_dense[spamm_dense_index(i, j, A.M, A.N)] += A_dense[spamm_dense_index(i, k, A.M, A.N)]*A_dense[spamm_dense_index(k, j, A.M, A.N)];
      }
    }
  }
#endif

  spamm_log("multiplying matrix with spamm\n", __FILE__, __LINE__);
  spamm_new(A.M, A.N, A.M_block, A.N_block, A.M_child, A.N_child, A.threshold, &A2);
  spamm_multiply(1.0, &A, &A, 1.0, &A2);

  spamm_log("comparing matrices\n", __FILE__, __LINE__);
  max_diff = 0;
  for (i = 0; i < A.M; ++i) {
    for (j = 0; j < A.N; ++j)
    {
      if (fabs(spamm_get(i, j, &A2)-A2_dense[spamm_dense_index(i, j, A.M, A.N)]) > max_diff)
      {
        max_diff = fabs(spamm_get(i, j, &A2)-A2_dense[spamm_dense_index(i, j, A.M, A.N)]);
        max_diff_i = i;
        max_diff_j = j;
      }
    }
  }

  if (max_diff > THRESHOLD)
  {
    printf("[multiply_spamm] biggest mismatch above threshold of %e: (A2[%i][%i] = %e) != (A2_dense[%i][%i] = %e), |diff| = %e\n",
        THRESHOLD,
        max_diff_i, max_diff_j, spamm_get(max_diff_i, max_diff_j, &A2),
        max_diff_i, max_diff_j, A2_dense[spamm_dense_index(max_diff_i, max_diff_j, A.M, A.N)],
        fabs(spamm_get(max_diff_i, max_diff_j, &A2)-A2_dense[spamm_dense_index(max_diff_i, max_diff_j, A.M, A.N)]));
    result = 1;
  }

  return result;
}
