#include "config.h"
#include <spamm.h>
#include <sparsekit.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define THRESHOLD 1e-14

int
main (int argc, char **argv)
{
  struct spamm_t A;
  struct spamm_t A2;

  floating_point_t *A_dense;
  floating_point_t *A2_dense;

#ifdef DGEMM
  floating_point_t alpha, beta;
#endif

  floating_point_t *A_CSR;
  int *A_i_CSR, *A_j_CSR;

  floating_point_t *C_CSR;
  int *C_i_CSR, *C_j_CSR, *work, *degree;

  struct spamm_tree_stats_t stats;

  double max_diff;
  int max_diff_i = 0;
  int max_diff_j = 0;

  int i, j;

#if ! defined(DGEMM)
  int k;
#endif

  int nonzero, C_nonzero;

  int ierr, job;
  int result = 0;

  struct timeval start, stop;


  if (argc != 2)
  {
    LOG2_FATAL("what file should I load?\n");
    exit(1);
  }

  LOG2_INFO("loading matrix\n");
  spamm_read_MM(argv[1], 8, 8, 2, 2, 1e-10, &A);
  spamm_tree_stats(&stats, &A);
  LOG_INFO("read %ix%i matrix, %ix%i blocks, %i nodes, %i dense blocks, depth = %i, tree = %i bytes, blocks = %i bytes (%1.1f%%), ",
      A.M, A.N, A.M_block, A.N_block, stats.number_nodes,
      stats.number_dense_blocks, A.tree_depth, stats.memory_tree,
      stats.memory_dense_blocks,
      (stats.memory_tree+stats.memory_dense_blocks)/(double) (A.M*A.N*sizeof(double))*100);
  printf("avg. sparsity = %1.1f%%\n", stats.average_sparsity*100);

  LOG2_INFO("converting tree to dense... ");
  gettimeofday(&start, NULL);
  spamm_spamm_to_dense(&A, &A_dense);
  gettimeofday(&stop, NULL);
  printf("time elapsed: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

  LOG2_INFO("converting dense to CSR\n");
  A_CSR = (floating_point_t*) malloc(sizeof(floating_point_t)*spamm_number_nonzero(&A));
  A_j_CSR = (int*) malloc(sizeof(int)*spamm_number_nonzero(&A));
  A_i_CSR = (int*) malloc(sizeof(int)*A.M);
  LOG_INFO("found %u nonzero elements\n", spamm_number_nonzero(&A));
  nonzero = spamm_number_nonzero(&A);
  dnscsr_single_(&A.M, &A.N, &nonzero, A_dense, &A.M, A_CSR, A_j_CSR, A_i_CSR, &ierr);

  LOG2_INFO("allocating A2_dense\n");
  A2_dense = (floating_point_t*) malloc(sizeof(floating_point_t)*A.M*A.N);
  for (i = 0; i < A.M; ++i) {
    for (j = 0; j < A.N; ++j)
    {
      A2_dense[spamm_dense_index(i, j, A.M, A.N)] = 0.0;
    }
  }

#ifdef DGEMM
  LOG2_INFO("multiplying matrix with BLAS to get reference product... ");
  alpha = 1.0;
  beta = 0.0;
  gettimeofday(&start, NULL);
  DGEMM("N", "N", &A.M, &A.N, &A.M, &alpha, A_dense, &A.M, A_dense, &A.N, &beta, A2_dense, &A.M);
  gettimeofday(&stop, NULL);
  printf("time elapsed: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);
#else
  LOG2_INFO("multiplying matrix directly to get reference product... ");
  gettimeofday(&start, NULL);
  for (i = 0; i < A.M; ++i) {
    for (j = 0; j < A.N; ++j) {
      for (k = 0; k < A.M; ++k)
      {
        A2_dense[spamm_dense_index(i, j, A.M, A.N)] += A_dense[spamm_dense_index(i, k, A.M, A.N)]*A_dense[spamm_dense_index(k, j, A.M, A.N)];
      }
    }
  }
  gettimeofday(&stop, NULL);
  printf("time elapsed: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);
#endif

  LOG2_INFO("multiplying matrix with sparsekit... ");
  gettimeofday(&start, NULL);
  degree = (int*) malloc(sizeof(int)*A.M);
  work = (int*) malloc(sizeof(int)*A.N);
  amubdg_(&A.M, &A.N, &A.N, A_j_CSR, A_i_CSR, A_j_CSR, A_i_CSR, degree, &C_nonzero, work);
  printf("need %i nonzeros in product... ", C_nonzero);
  C_CSR = (floating_point_t*) malloc(sizeof(floating_point_t)*C_nonzero);
  C_j_CSR = (int*) malloc(sizeof(int)*C_nonzero);
  C_i_CSR = (int*) malloc(sizeof(int)*A.M);
  job = 1;
  amub_single_(&A.M, &A.N, &job, A_CSR, A_j_CSR, A_i_CSR, A_CSR, A_j_CSR, A_i_CSR, C_CSR, C_j_CSR, C_i_CSR, &C_nonzero, work, &ierr);
  gettimeofday(&stop, NULL);
  printf("time elapsed: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

  LOG2_INFO("multiplying matrix with spamm\n");
  spamm_new(A.M, A.N, A.M_block, A.N_block, A.M_child, A.N_child, A.threshold, &A2);
  gettimeofday(&start, NULL);
  spamm_multiply(tree, 1.0, &A, &A, 1.0, &A2);
  LOG_INFO("product has %u nonzero elements\n", spamm_number_nonzero(&A2));
  gettimeofday(&stop, NULL);
  LOG_INFO("total spamm time elapsed: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

  LOG2_INFO("comparing matrices\n");
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
