#include <spamm.h>
#include <stdlib.h>
#include <sys/time.h>

int
main (int argc, char **argv)
{
  struct spamm_t A;
  struct spamm_t A2;

  floating_point_t *A_dense;
  floating_point_t *A2_dense;

  struct spamm_tree_stats_t stats;
  struct timeval start, stop;

  int result = 0;

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

  /* Purify using SP2. */
  LOG2_INFO("SP2 on SpAMM... ");
  gettimeofday(&start, NULL);

  gettimeofday(&stop, NULL);
  printf("time elapsed: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

  /* Release memory. */
  spamm_delete(&A);
  free(A_dense);

  return result;
}
