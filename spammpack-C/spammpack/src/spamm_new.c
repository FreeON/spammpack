#include "config.h"
#include "spamm.h"
#include "spamm_types_private.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Initialize new matrix object.
 *
 * @param tier The tier.
 * @param depth The depth of the matrix tree.
 * @param N_padded The padded size of the matrix.
 * @param M_lower The lower row index of this submatrix.
 * @param M_upper The upper row index of this submatrix.
 * @param N_lower The lower column index of this submatrix.
 * @param N_upper The upper column index of this submatrix.
 *
 * @return A pointer to the matrix.
 */
struct spamm_hashed_t *
spamm_hashed_new (const unsigned int tier,
    const unsigned int kernel_tier,
    const unsigned int depth,
    const unsigned int N_padded,
    const unsigned int M_lower,
    const unsigned int M_upper,
    const unsigned int N_lower,
    const unsigned int N_upper)
{
  unsigned int i;
  struct spamm_hashed_t *A;

  /* Allocate memory. */
  A = calloc(1, sizeof(struct spamm_hashed_t));

  /* Store kernel_tier. */
  A->kernel_tier = kernel_tier;

  /* Store tier. */
  A->tier = tier;

  /* Store depth. */
  A->depth = depth;

  /* Store matrix size. */
  A->N_padded = N_padded;

  /* Store bounding box. */
  A->M_lower = M_lower;
  A->M_upper = M_upper;
  A->N_lower = N_lower;
  A->N_upper = N_upper;

  /* Create the tier hash tables. */
  A->tier_hashtable = (struct spamm_hashtable_t**) malloc(sizeof(struct spamm_hashtable_t*)*(A->kernel_tier-A->tier+1));
  for (i = A->tier; i <= A->kernel_tier; i++)
  {
    A->tier_hashtable[i-A->tier] = spamm_hashtable_new();
  }

  return A;
}

/** Create a new recursive matrix object.
 *
 * @param number_dimensions The number of dimensions.
 * @param N An array of the number of rows/columns of dense input matrix.
 * @param blocksize The size of the dense matrix blocks.
 *
 * @return A pointer to the matrix.
 */
struct spamm_recursive_t *
spamm_recursive_new (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int N_contiguous)
{
  unsigned int dim;
  struct spamm_recursive_t *A = NULL;
  double x, x_N;

  for (dim = 0; dim < number_dimensions; dim++)
  {
    if (N[dim] <= 0)
    {
      SPAMM_FATAL("N[%u] <= 0\n", dim);
    }
  }

  /* Allocate memory. */
  A = calloc(1, sizeof(struct spamm_recursive_t));

  /* Store the blocksize. */
  A->N_contiguous = N_contiguous;

  /* Pad to powers of M_child x N_child. */
  x = 0;
  for (dim = 0; dim < number_dimensions; dim++)
  {
    x_N = (log(N[dim]) > log(N_contiguous) ? log(N[dim]) - log(N_contiguous) : 0)/log(2);
    if (x_N > x)
    {
      x = x_N;
    }
  }

  /* The ceil() function can lead to a depth that is one tier too large
   * because of numerical errors in the calculation of x. We need to check
   * whether the depth is appropriate.
   */
  A->depth = (unsigned int) ceil(x);

  /* Double check depth. */
  if (A->depth >= 1)
  {
    for (dim = 0; dim < A->number_dimensions; dim++)
    {
      if ((int) (N_contiguous*ipow(2, A->depth-1)) < A->N[dim])
      {
        A->depth++;
        break;
      }
    }
    A->depth--;
  }

  /* Store the number of dimensions. */
  A->number_dimensions = number_dimensions;

  /* Set matrix size. */
  A->N = calloc(number_dimensions, sizeof(unsigned int));
  for (dim = 0; dim < number_dimensions; dim++)
  {
    A->N[dim] = N[dim];
  }

  /* Set padded matrix size. */
  A->N_padded = (int) (N_contiguous*ipow(2, A->depth));

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
struct spamm_hashed_data_t *
spamm_hashed_new_data (const unsigned int tier, const unsigned int index_2D, const enum spamm_layout_t layout)
{
  int i, j;
  struct spamm_hashed_data_t *data;

#ifdef HAVE_POSIX_MEMALIGN
  int result;

  /* Allocate data. */
  if ((result = posix_memalign((void**) &data, SPAMM_PAGE_ALIGNMENT, sizeof(struct spamm_hashed_data_t))) != 0)
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

#else
  /* Allocate data (this is with unknown alignment, i.e. it is aligned to
   * whatever malloc() aligns it to. */
  data = (struct spamm_hashed_data_t*) calloc(1, sizeof(struct spamm_hashed_data_t));
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
 * @param number_dimensions The number of dimensions.
 * @param N_contiguous The size of the contiguous submatrix block.
 * @param N_linear The size of the matrix when switching to a linear tree.
 * @param N_lower An array of the lowest row index of this submatrix node.
 * @param N_upper An array of the lowest row index of this submatrix node.
 *
 * @return A pointer to the newly allocated node.
 */
struct spamm_recursive_node_t *
spamm_recursive_new_node (const unsigned int tier,
    const unsigned int number_dimensions,
    const unsigned int N_contiguous,
    const unsigned int N_linear,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper)
{
  unsigned int dim;
  struct spamm_recursive_node_t *node = NULL;

  /* Allocate memory. */
  node = calloc(1, sizeof(struct spamm_recursive_node_t));

  node->tier = tier;
  node->number_dimensions = number_dimensions;
  node->N_contiguous = N_contiguous;
  node->N_linear = N_linear;

  node->N_lower = calloc(number_dimensions, sizeof(unsigned int));
  node->N_upper = calloc(number_dimensions, sizeof(unsigned int));

  /* Allocate child matrix. */
  node->child = calloc(ipow(2, number_dimensions), sizeof(struct spamm_recursive_node_t*));

  for (dim = 0; dim < number_dimensions; dim++)
  {
    node->N_lower[dim] = N_lower[dim];
    node->N_upper[dim] = N_upper[dim];
  }

  return node;
}

/** Initialize a new matrix object.
 *
 * Two settings determine when and if the tree is stored in linear or
 * hierarchical format. If N_linear == N_contiguous, then the tree format will
 * switch to linear tree format at N_contiguous. If N_linear < N_contiguous,
 * then the whole tree will be stored hierarchically and the spamm condition
 * is applied to the SpAMM chunks.
 *
 * @param number_dimensions The number of dimensions of this matrix.
 * @param N The number of rows/columns of the matrix. This array has to have
 * a size of number_dimensions.
 * @param N_linear The tier at which to switch from hierarchical to linear
 * tree.
 * @param contiguous_tier The tier at which to store contiguous submatrix
 * blocks in the hierarhical tree.
 * @param layout The storage layout of the matrix elements.
 *
 * @return The newly allocated matrix. This matrix has to be freed by calling
 * spamm_delete().
 */
struct spamm_matrix_t *
spamm_new (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int N_linear,
    const unsigned int contiguous_tier,
    const enum spamm_layout_t layout)
{
  int dim;
  unsigned int N_temp;
  struct spamm_matrix_t *A = NULL;
  double x, x_N;

  for (dim = 0; dim < number_dimensions; dim++)
  {
    if (N[dim] == 0)
    {
      SPAMM_FATAL("N[%u] == 0\n", dim);
    }
  }

  /* Allocate memory. */
  A = calloc(1, sizeof(struct spamm_matrix_t));

  /* Store the number of dimensions. */
  A->number_dimensions = number_dimensions;

  /* Store matrix dimensions. */
  A->N = calloc(number_dimensions, sizeof(unsigned int));
  for (dim = 0; dim < number_dimensions; dim++)
  {
    A->N[dim] = N[dim];
  }

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
  x = 0;
  for (dim = 0; dim < number_dimensions; dim++)
  {
    /* Make sure we pad at least to the extend that we can store that linear
     * kernel matrix. */
    if (N[dim] < SPAMM_N_KERNEL)
    {
      N_temp = SPAMM_N_KERNEL;
    }

    else
    {
      N_temp = N[dim];
    }

    x_N = (log(N_temp) > log(SPAMM_N_BLOCK) ? log(N_temp) - log(SPAMM_N_BLOCK) : 0)/log(2);
    if (x_N > x)
    {
      x = x_N;
    }
  }

  /* The ceil() function can lead to a depth that is one tier too large
   * because of numerical errors in the calculation of x. We need to check
   * whether the depth is appropriate.
   */
  A->depth = (unsigned int) ceil(x);

  /* Double check depth. */
  if (A->depth >= 1)
  {
    for (dim = 0; dim < A->number_dimensions; dim++)
    {
      if ((int) (SPAMM_N_BLOCK*ipow(2, A->depth-1)) < A->N[dim])
      {
        A->depth++;
        break;
      }
    }
    A->depth--;
  }

  /* Adjust tree to kernel depth. */
  if (A->depth < SPAMM_KERNEL_DEPTH)
  {
    /* We should have already made sure that the matrix is big enough to fit
     * the kernel matrix. */
    SPAMM_FATAL("there is some logic error in this function\n");
  }

  /* Adjust the linear depth. */
  if (linear_tier+SPAMM_KERNEL_DEPTH > A->depth)
  {
    SPAMM_WARN("linear tier (%u) + kernel depth (%u) is greater than depth (%u)\n", linear_tier, SPAMM_KERNEL_DEPTH, A->depth);
    A->linear_tier = A->depth+1;
  }

  else
  {
    A->linear_tier = linear_tier;
  }

  /* Adjust the contiguous tier depth. */
  if (contiguous_tier > A->depth)
  {
    SPAMM_WARN("contiguous tier (%u) is greater than depth (%u), adjusting to depth\n", contiguous_tier, A->depth);
    A->contiguous_tier = A->depth;
  }

  else
  {
    A->contiguous_tier = contiguous_tier;
  }

  if (A->contiguous_tier >= A->linear_tier)
  {
    /* Reset contiguous_tier, it will not get used. */
    A->contiguous_tier = A->depth-SPAMM_KERNEL_DEPTH;
  }

  /* Set padded matrix size. */
  A->N_padded = (int) (SPAMM_N_BLOCK*ipow(2, A->depth));

  /* Adjust hashed size. In case linear_tier is larger than the depth, we do
   * not adjust N_linear, i.e. we disable the linear tree. */
  if (A->linear_tier < A->depth)
  {
    A->N_linear = (int) (SPAMM_N_BLOCK*ipow(2, A->depth-A->linear_tier));
  }

  /* Adjust contiguous size. */
  A->N_contiguous = (int) (SPAMM_N_BLOCK*ipow(2, A->depth-A->contiguous_tier));

  if (A->N_linear > A->N_padded)
  {
    A->N_linear = A->N_padded;
  }

  /* Set the kernel tier. */
  A->kernel_tier = A->depth-SPAMM_KERNEL_DEPTH;

  /* Done. */
  return A;
}
