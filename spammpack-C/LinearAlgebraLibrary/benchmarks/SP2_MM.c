#include <spamm.h>
#include <stdlib.h>

int
main (int argc, char **argv)
{
  struct spamm_t *A;
  struct spamm_t *A2;

  float_t *A_dense;
  float_t *A2_dense;

  int result = 0;

  if (argc != 2)
  {
    spamm_log("what file should I load?\n", __FILE__, __LINE__);
    exit(1);
  }

  spamm_log("loading matrix\n", __FILE__, __LINE__);
  spamm_read_MM(argv[1], 8, 8, 2, 2, 1e-10, &A);
  spamm_tree_stats(&stats, &A);
  spamm_log("read %ix%i matrix, %ix%i blocks, %i nodes, %i dense blocks, depth = %i, tree = %i bytes, blocks = %i bytes (%1.1f%%), ",
      __FILE__, __LINE__, A.M, A.N, A.M_block, A.N_block, stats.number_nodes, stats.number_dense_blocks, A.tree_depth,
      stats.memory_tree, stats.memory_dense_blocks, (stats.memory_tree+stats.memory_dense_blocks)/(double) (A.M*A.N*sizeof(double))*100);
  printf("avg. sparsity = %1.1f%%\n", stats.average_sparsity*100);

  spamm_log("converting tree to dense... ", __FILE__, __LINE__);
  gettimeofday(&start, NULL);
  spamm_spamm_to_dense(&A, &A_dense);
  gettimeofday(&stop, NULL);
  printf("time elapsed: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

  return result;
}
