#include "spamm.h"
#include "spamm_types_private.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Initialize a new matrix object.
 *
 * @param M The number of rows of the matrix.
 * @param N The number of columns of the matrix.
 * @param N_block The size of a basic matrix block as NxN single matrix elements.
 * @param number_contiguous_tiers The number of tiers that are allocated
 * contiguously for the computational kernel.
 * @param number_hashed_tiers The number of tiers from the bottom of the
 * matrix tree that are stored in hashed format. This number can be zero in
 * which case the whole tree is stored in hierarchical format.
 * @param kernel_tier The tier at which to store the dense kernel matrix
 * blocks.
 *
 * @return The newly allocated matrix. This matrix has to be freed by calling
 * spamm_delete().
 */
struct spamm_matrix_t *
spamm_new (const unsigned int M, const unsigned int N,
    const unsigned int N_block,
    const unsigned int number_contiguous_tiers,
    const unsigned int number_hashed_tiers,
    const unsigned int kernel_tier,
    const enum spamm_layout_t layout)
{
  struct spamm_matrix_t *A;
  double x, x_M, x_N;

  if (kernel_tier == 0)
  {
    spamm_error_fatal(__FILE__, __LINE__, "kernel_tier has to be greater 0\n");
  }

  A = calloc(1, sizeof(struct spamm_matrix_t*));
  A->number_hashed_tiers = number_hashed_tiers;
  A->kernel_tier = kernel_tier;
  A->blocksize = (unsigned int) pow(2, kernel_tier);

  /* Pad to powers of M_child x N_child. */
  x_M = (log(M) > log(N_block) ? log(M) - log(N_block) : 0)/log(2);
  x_N = (log(N) > log(N_block) ? log(N) - log(N_block) : 0)/log(2);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  /* The ceil() function can lead to a depth that is one tier too large
   * because of numerical errors in the calculation of x. We need to check
   * whether the depth is appropriate.
   */
  A->depth = (unsigned int) ceil(x);

  /* Double check depth. */
  if (A->depth >= 1 && ((int) (N_block*pow(2, A->depth-1)) >= M && (int) (N_block*pow(2, A->depth-1)) >= N))
  {
    (A->depth)--;
  }

  /* Adjust tree to kernel depth. */
  if (A->depth < number_contiguous_tiers) { A->depth = number_contiguous_tiers; }

  /* Set matrix size. */
  A->M = M;
  A->N = N;

  /* Set padded matrix size. */
  A->N_padded = (int) (N_block*pow(2, A->depth));

  /* Set the kernel tier. */
  A->kernel_tier = A->depth-number_contiguous_tiers;

  /* Set the layout. */
  switch (layout)
  {
    case row_major:
    case column_major:
    case Z_curve:
    case dense_column_major:
      A->layout = layout;
      break;

    default:
      spamm_error_fatal(__FILE__, __LINE__, "unknown layout (%i)\n", layout);
      exit(1);
      break;
  }

  return A;
}

/** Initialize new matrix object.
 *
 * @param M Number of rows of dense input matrix.
 * @param N Number of columns of dense input matrix.
 *
 * @return A pointer to the matrix.
 */
struct spamm_hashed_t *
spamm_hashed_new (const unsigned int M, const unsigned int N, const enum spamm_layout_t layout)
{
  struct spamm_hashed_t *A;
  double x, x_M, x_N;
  unsigned int tier;

  if (M <= 0)
  {
    fprintf(stderr, "M <= 0\n");
    exit(1);
  }

  if (N <= 0)
  {
    fprintf(stderr, "N <= 0\n");
    exit(1);
  }

  /* Allocate memory. */
  A = (struct spamm_hashed_t*) malloc(sizeof(struct spamm_hashed_t));

  /* Set the layout. */
  switch (layout)
  {
    case row_major:
    case column_major:
    case Z_curve:
    case dense_column_major:
      A->layout = layout;
      break;

    default:
      fprintf(stderr, "[spamm new] unknown layout (%i)\n", layout);
      exit(1);
      break;
  }

  /* Pad to powers of M_child x N_child. */
  x_M = (log(M) > log(SPAMM_N_BLOCK) ? log(M) - log(SPAMM_N_BLOCK) : 0)/log(2);
  x_N = (log(N) > log(SPAMM_N_BLOCK) ? log(N) - log(SPAMM_N_BLOCK) : 0)/log(2);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  /* The ceil() function can lead to a depth that is one tier too large
   * because of numerical errors in the calculation of x. We need to check
   * whether the depth is appropriate.
   */
  A->depth = (unsigned int) ceil(x);

  /* Double check depth. */
  if (A->depth >= 1 && ((int) (SPAMM_N_BLOCK*pow(2, A->depth-1)) >= M && (int) (SPAMM_N_BLOCK*pow(2, A->depth-1)) >= N))
  {
    (A->depth)--;
  }

  /* Adjust tree to kernel depth. */
  if (A->depth < SPAMM_KERNEL_DEPTH) { A->depth = SPAMM_KERNEL_DEPTH; }

  /* Set matrix size. */
  A->M = M;
  A->N = N;

  /* Set padded matrix size. */
  A->N_padded = (int) (SPAMM_N_BLOCK*pow(2, A->depth));

  /* Set the kernel tier. */
  A->kernel_tier = A->depth-SPAMM_KERNEL_DEPTH;

  /* Create the tier hash tables. */
  A->tier_hashtable = (struct spamm_hashtable_t**) malloc(sizeof(struct spamm_hashtable_t*)*(A->kernel_tier+1));
  for (tier = 0; tier <= A->kernel_tier; tier++)
  {
    A->tier_hashtable[tier] = spamm_hashtable_new();
  }

  return A;
}

/** Create a new recursive matrix object.
 *
 * @param M Number of rows of dense input matrix.
 * @param N Number of columns of dense input matrix.
 * @param blocksize The size of the dense matrix blocks.
 *
 * @return A pointer to the matrix.
 */
struct spamm_recursive_t *
spamm_recursive_new (const unsigned int M, const unsigned int N, const unsigned int blocksize)
{
  struct spamm_recursive_t *A = NULL;
  double x, x_M, x_N;

  if (M <= 0)
  {
    fprintf(stderr, "M <= 0\n");
    exit(1);
  }

  if (N <= 0)
  {
    fprintf(stderr, "N <= 0\n");
    exit(1);
  }

  /* Allocate memory. */
  A = calloc(1, sizeof(struct spamm_recursive_t));

  /* Store the blocksize. */
  A->blocksize = blocksize;

  /* Pad to powers of M_child x N_child. */
  x_M = (log(M) > log(blocksize) ? log(M) - log(blocksize) : 0)/log(2);
  x_N = (log(N) > log(blocksize) ? log(N) - log(blocksize) : 0)/log(2);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  /* The ceil() function can lead to a depth that is one tier too large
   * because of numerical errors in the calculation of x. We need to check
   * whether the depth is appropriate.
   */
  A->depth = (unsigned int) ceil(x);

  /* Double check depth. */
  if (A->depth >= 1 && ((int) (blocksize*pow(2, A->depth-1)) >= M && (int) (blocksize*pow(2, A->depth-1)) >= N))
  {
    (A->depth)--;
  }

  /* Set matrix size. */
  A->M = M;
  A->N = N;

  /* Set padded matrix size. */
  A->N_padded = (int) (blocksize*pow(2, A->depth));

  return A;
}

/** Allocate a new node of a matrix tree.
 *
 * @param tier The tier this node will be on.
 * @param index_2D The 2D linear matrix index of this node.
 *
 * @return A pointer to the newly allocated node.
 */
struct spamm_hashed_node_t *
spamm_hashed_new_node (const unsigned int tier, const unsigned int index_2D)
{
  struct spamm_hashed_node_t *node = (struct spamm_hashed_node_t*) malloc(sizeof(struct spamm_hashed_node_t));

  node->tier = tier;
  node->index_2D = index_2D;

  node->norm = 0.0;
  node->norm2 = 0.0;

  return node;
}

/** Allocate a new data node of a matrix tree.
 *
 * @param tier The tier this node will be on.
 * @param index_2D The 2D linear matrix index of this node.
 * @param layout The layout of the basic matrix blocks.
 *
 * @return A pointer to the newly allocated node.
 */
struct spamm_data_t *
spamm_new_block (const unsigned int tier, const unsigned int index_2D, const enum spamm_layout_t layout)
{
  int i, j;
  struct spamm_data_t *data;

#ifdef HAVE_POSIX_MEMALIGN
  int result;

  /* Allocate data. */
  if ((result = posix_memalign((void**) &data, SPAMM_PAGE_ALIGNMENT, sizeof(struct spamm_data_t))) != 0)
  {
    switch (result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        exit(1);
        break;

      case ENOMEM:
        printf("There was insufficient memory to fulfill the allocation request.\n");
        exit(1);
        break;

      default:
        printf("unknown error code: %i\n", result);
        exit(1);
        break;
    }
  }

  /* Set matrix elements to zero. */
  for (i = 0; i < SPAMM_N_KERNEL; i++) {
    for (j = 0; j < SPAMM_N_KERNEL; j++)
    {
      data->block_dense[i*SPAMM_N_KERNEL+j] = 0.0;
      data->block_dense_store[i*SPAMM_N_KERNEL+j] = 0.0;
      data->block_dense_transpose[i*SPAMM_N_KERNEL+j] = 0.0;

      data->block_dense_dilated[4*(i*SPAMM_N_KERNEL+j)+0] = 0.0;
      data->block_dense_dilated[4*(i*SPAMM_N_KERNEL+j)+1] = 0.0;
      data->block_dense_dilated[4*(i*SPAMM_N_KERNEL+j)+2] = 0.0;
      data->block_dense_dilated[4*(i*SPAMM_N_KERNEL+j)+3] = 0.0;
    }
  }

  for (i = 0; i < SPAMM_N_KERNEL_BLOCKED*SPAMM_N_KERNEL_BLOCKED; i++)
  {
    data->norm[i] = 0.0;
    data->norm2[i] = 0.0;
  }

  data->node_norm = 0.0;
  data->node_norm2 = 0.0;

  for (i = 0; i < 8; i++)
  {
    data->norm_upper[i] = 0.0;
    data->norm_upper_transpose[i] = 0.0;
  }
#else
  /* Allocate data (this is with unknown alignment, i.e. it is aligned to
   * whatever malloc() aligns it to. */
  data = (struct spamm_data_t*) calloc(1, sizeof(struct spamm_data_t));
#endif

  switch (layout)
  {
    case row_major:
    case column_major:
    case Z_curve:
    case dense_column_major:
      data->layout = layout;
      break;

    default:
      fprintf(stderr, "[spamm new block] unknown layout (%i)\n", layout);
      exit(1);
      break;
  }

  /* Set some information on the data block. */
  data->tier = tier;
  data->index_2D = index_2D;

  return data;
}

/** Allocate a new node of a recursive matrix tree.
 *
 * @param tier The tier this node will be on.
 *
 * @return A pointer to the newly allocated node.
 */
struct spamm_recursive_node_t *
spamm_recursive_new_node (const unsigned int tier, const unsigned int blocksize)
{
  struct spamm_recursive_node_t *node = NULL;

  node = calloc(1, sizeof(struct spamm_recursive_node_t));
  node->tier = tier;
  node->blocksize = blocksize;

  return node;
}
