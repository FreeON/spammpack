#include "config.h"
#include "spamm.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>

/** Multiply a chunk with a scalar.
 *
 * @param alpha The factor.
 * @param chunk The chunk.
 *
 * @return The new squared norm.
 */
float
spamm_chunk_multiply_scalar (const float alpha,
    spamm_chunk_t *chunk)
{
  unsigned int i;
  unsigned int N_contiguous;
  unsigned int number_dimensions;

  float *A;
  float *norm;
  float *norm2;

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  N_contiguous = spamm_chunk_get_N_contiguous(chunk);
  A = spamm_chunk_get_matrix(chunk);

  for (i = 0; i < ipow(N_contiguous, number_dimensions); i++)
  {
    A[i] *= alpha;
  }

  norm = spamm_chunk_get_norm(chunk);
  norm2 = spamm_chunk_get_norm2(chunk);
  for (i = 0; i < 0; i++)
  {
    norm[i] *= alpha;
    norm2[i] *= alpha*alpha;
  }

  return norm2[0];
}

//#define PRINT_DEBUG
float
spamm_chunk_multiply (const float tolerance,
    const float alpha,
    spamm_chunk_t *chunk_A,
    spamm_chunk_t *chunk_B,
    spamm_chunk_t *chunk_C,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const unsigned int depth,
    const unsigned int N_block,
    const unsigned int linear_index_A,
    const unsigned int linear_index_B,
    const unsigned int linear_index_C,
    struct spamm_timer_t *timer,
    sgemm_func sgemm,
    const enum spamm_kernel_t kernel)
{
  unsigned int i, j, k;

  unsigned int number_dimensions_A;
  unsigned int number_dimensions_B;
  unsigned int number_dimensions_C;

  unsigned int new_linear_index_A;
  unsigned int new_linear_index_B;
  unsigned int new_linear_index_C;

#ifdef PRINT_DEBUG
  int dim;

  unsigned int *N_lower;
  unsigned int *N_upper;
#endif

  float *norm_A;
  float *norm_B;

  float *norm_C;
  float *norm2_C;

  float alpha_sgemm = alpha;
  float beta = 1.0;

  int N_contiguous;

  float *matrix_A;
  float *matrix_B;
  float *matrix_C;

  short use_linear_tree;

  use_linear_tree = *spamm_chunk_get_use_linear_tree(chunk_A);

  if (use_linear_tree)
  {
    spamm_linear_multiply(tolerance, alpha, chunk_A, chunk_B, chunk_A,
        chunk_B, beta, chunk_C, chunk_C, timer, kernel);
  }

  else
  {
    matrix_A = spamm_chunk_get_matrix(chunk_A);
    matrix_B = spamm_chunk_get_matrix(chunk_B);
    matrix_C = spamm_chunk_get_matrix(chunk_C);

    N_contiguous = spamm_chunk_get_N_contiguous(chunk_A);

    if (sgemm)
    {
      sgemm("N", "N", &N_contiguous, &N_contiguous, &N_contiguous,
          &alpha_sgemm, matrix_A, &N_contiguous, matrix_B, &N_contiguous,
          &beta, matrix_C, &N_contiguous);
    }

    else
    {
      /* Braindead multiply in nested loops. */
      for (i = 0; i < N_contiguous; i++) {
        for (j = 0; j < N_contiguous; j++) {
          for (k = 0; k < N_contiguous; k++)
          {
            matrix_C[spamm_index_column_major(i, j, N_contiguous, N_contiguous)] += alpha
              *matrix_A[spamm_index_column_major(i, k, N_contiguous, N_contiguous)]
              *matrix_B[spamm_index_column_major(k, j, N_contiguous, N_contiguous)];
          }
        }
      }
    }
  }
}

void
spamm_chunk_add (const float alpha,
    spamm_chunk_t **A,
    const float beta,
    spamm_chunk_t *B)
{
  float *A_matrix;
  float *B_matrix;

  unsigned int i;
  unsigned int N_contiguous;
  unsigned int number_dimensions;

  A_matrix = spamm_chunk_get_matrix(*A);
  B_matrix = spamm_chunk_get_matrix(B);

  N_contiguous = spamm_chunk_get_N_contiguous(B);
  number_dimensions = *spamm_chunk_get_number_dimensions(B);

  for (i = 0; i < ipow(N_contiguous, number_dimensions); i++)
  {
    A_matrix[i] = alpha*A_matrix[i]+beta*B_matrix[i];
  }
}

void
spamm_chunk_copy (spamm_chunk_t **A,
    const float beta,
    spamm_chunk_t *B)
{
  unsigned int *number_dimensions = spamm_chunk_get_number_dimensions(B);
  unsigned int *number_tiers = spamm_chunk_get_number_tiers(B);
  unsigned int *N = spamm_chunk_get_N(B);
  unsigned int *N_lower = spamm_chunk_get_N_lower(B);
  unsigned int *N_upper = spamm_chunk_get_N_upper(B);

  float *norm_A;
  float *norm_B;
  float *norm2_A;
  float *norm2_B;

  float *A_matrix;
  float *B_matrix;

  unsigned int N_contiguous;

  unsigned int i;

  spamm_delete_chunk(A);

  *A = spamm_new_chunk(*number_dimensions, (*number_tiers == 1 ? 0 : 1), N, N_lower, N_upper);

  norm_A = spamm_chunk_get_norm(*A);
  norm_B = spamm_chunk_get_norm(B);
  norm2_A = spamm_chunk_get_norm2(*A);
  norm2_B = spamm_chunk_get_norm2(B);

  for (i = 0; i < ipow(*number_dimensions, *number_tiers); i++)
  {
    norm_A[i] = beta*norm_B[i];
    norm2_A[i] = beta*beta*norm2_B[i];
  }

  A_matrix = spamm_chunk_get_matrix(*A);
  B_matrix = spamm_chunk_get_matrix(B);

  N_contiguous = spamm_chunk_get_N_contiguous(B);

  for (i = 0; i < N_contiguous; i++)
  {
    A_matrix[i] = beta*B_matrix[i];
  }
}

/** Set an element in a SpAMM chunk.
 *
 * @param i The row/column index array.
 * @param Aij The value of the matrix element.
 * @param chunk The SpAMM chunk.
 */
void
spamm_chunk_set (const unsigned int *const i,
    const float Aij,
    spamm_chunk_t *chunk)
{
  int dim;

  short use_linear_tree;

  unsigned int tier;
  unsigned int linear_index;

  unsigned int number_dimensions;
  unsigned int number_tiers;

  unsigned int *N_lower;
  unsigned int *N_upper;

  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  unsigned int *new_i;

  float *norm;
  float *norm2;

  float *A;

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  number_tiers = *spamm_chunk_get_number_tiers(chunk);
  use_linear_tree = *spamm_chunk_get_use_linear_tree(chunk);

  N_lower = spamm_chunk_get_N_lower(chunk);
  N_upper = spamm_chunk_get_N_upper(chunk);

  new_N_lower = calloc(number_dimensions, sizeof(unsigned int));
  new_N_upper = calloc(number_dimensions, sizeof(unsigned int));

  for (dim = 0; dim < number_dimensions; dim++)
  {
    new_N_lower[dim] = N_lower[dim];
    new_N_upper[dim] = N_upper[dim];
  }

  /* Z-curve ordering down to SPAMM_N_KERNEL. */
  for (tier = 0, linear_index = 0; tier < number_tiers; tier++)
  {
    norm = spamm_chunk_get_tier_norm(tier, chunk);
    norm2 = spamm_chunk_get_tier_norm2(tier, chunk);

    /* Update norm. */
    norm2[linear_index] += Aij*Aij;
    norm[linear_index] = sqrt(norm2[linear_index]);

    if (tier+1 < number_tiers)
    {
      /* Recurse. */
      linear_index <<= number_dimensions;

      for (dim = 0; dim < number_dimensions; dim++)
      {
        if (i[dim] < new_N_lower[dim]+(new_N_upper[dim]-new_N_lower[dim])/2)
        {
          new_N_lower[dim] = new_N_lower[dim];
          new_N_upper[dim] = new_N_lower[dim]+(new_N_upper[dim]-new_N_lower[dim])/2;
        }

        else
        {
          new_N_upper[dim] = new_N_upper[dim];
          new_N_lower[dim] = new_N_lower[dim]+(new_N_upper[dim]-new_N_lower[dim])/2;
          linear_index |= (1 << dim);
        }
      }
    }

    else
    {
      break;
    }
  }

  A = spamm_chunk_get_matrix(chunk);
  if (use_linear_tree)
  {
    A[linear_index*ipow(SPAMM_N_KERNEL, number_dimensions)*sizeof(float)
      +spamm_index_kernel_block(i[0]-new_N_lower[0], i[1]-new_N_lower[1], row_major)] = Aij;
  }

  else
  {
    new_i = calloc(number_dimensions, sizeof(unsigned int));
    for (dim = 0; dim < number_dimensions; dim++)
    {
      new_i[dim] = i[dim]-N_lower[dim];
    }
    A[spamm_index_column_major_2(number_dimensions, N_upper[0]-N_lower[0], new_i)] = Aij;
    free(new_i);
  }

  free(new_N_lower);
  free(new_N_upper);
}

/** Get an element from a SpAMM chunk.
 *
 * @param i The row/column index array.
 * @param chunk The chunk.
 *
 * @return The matrix element.
 */
float
spamm_chunk_get (const unsigned int *i,
    spamm_chunk_t *chunk)
{
  float Aij = 0;

  int dim;

  short use_linear_tree;

  unsigned int tier;
  unsigned int linear_index;

  unsigned int number_dimensions;
  unsigned int number_tiers;

  unsigned int *N_lower;
  unsigned int *N_upper;

  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  unsigned int *new_i;

  float *norm;
  float *norm2;

  float *A;

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  number_tiers = *spamm_chunk_get_number_tiers(chunk);
  use_linear_tree = *spamm_chunk_get_use_linear_tree(chunk);

  N_lower = spamm_chunk_get_N_lower(chunk);
  N_upper = spamm_chunk_get_N_upper(chunk);

  new_N_lower = calloc(number_dimensions, sizeof(unsigned int));
  new_N_upper = calloc(number_dimensions, sizeof(unsigned int));

  for (dim = 0; dim < number_dimensions; dim++)
  {
    new_N_lower[dim] = N_lower[dim];
    new_N_upper[dim] = N_upper[dim];
  }

  /* Z-curve ordering down to SPAMM_N_KERNEL. */
  for (tier = 0, linear_index = 0; tier < number_tiers; tier++)
  {
    norm = spamm_chunk_get_tier_norm(tier, chunk);
    norm2 = spamm_chunk_get_tier_norm2(tier, chunk);

    /* Update norm. */
    norm2[linear_index] += Aij*Aij;
    norm[linear_index] = sqrt(norm2[linear_index]);

    if (tier+1 < number_tiers)
    {
      /* Recurse. */
      linear_index <<= number_dimensions;

      for (dim = 0; dim < number_dimensions; dim++)
      {
        if (i[dim] < new_N_lower[dim]+(new_N_upper[dim]-new_N_lower[dim])/2)
        {
          new_N_lower[dim] = new_N_lower[dim];
          new_N_upper[dim] = new_N_lower[dim]+(new_N_upper[dim]-new_N_lower[dim])/2;
        }

        else
        {
          new_N_upper[dim] = new_N_upper[dim];
          new_N_lower[dim] = new_N_lower[dim]+(new_N_upper[dim]-new_N_lower[dim])/2;
          linear_index |= (1 << dim);
        }
      }
    }

    else
    {
      break;
    }
  }

  A = spamm_chunk_get_matrix(chunk);
  if (use_linear_tree)
  {
    Aij = A[linear_index*ipow(SPAMM_N_KERNEL, number_dimensions)*sizeof(float)
      +spamm_index_kernel_block(i[0]-new_N_lower[0], i[1]-new_N_lower[1], row_major)];
  }

  else
  {
    new_i = calloc(number_dimensions, sizeof(unsigned int));
    for (dim = 0; dim < number_dimensions; dim++)
    {
      new_i[dim] = i[dim]-N_lower[dim];
    }
    Aij = A[spamm_index_column_major_2(number_dimensions, N_upper[0]-N_lower[0], new_i)];
    free(new_i);
  }

  free(new_N_lower);
  free(new_N_upper);

  return Aij;
}

/** Pad memory address to some alignment.
 *
 * @param address Address to pad.
 * @param alignment The byte boundary to align to.
 *
 * @return The aligned address.
 */
size_t
spamm_chunk_pad (const size_t address,
    const size_t alignment)
{
  size_t new_address = address & (-alignment);

  if (new_address != address)
  {
    new_address += alignment;
  }

  return new_address;
}

/** Get the number_dimensions.
 *
 * @param chunk The chunk.
 *
 * @return The number_dimensions.
 */
unsigned int *
spamm_chunk_get_number_dimensions (spamm_chunk_t *chunk)
{
  return (unsigned int*) chunk;
}

/** Get N_block.
 *
 * @param chunk The chunk.
 *
 * @return number_tiers
 */
unsigned int *
spamm_chunk_get_number_tiers (spamm_chunk_t *chunk)
{
  return (unsigned int*) ((intptr_t) chunk + sizeof(unsigned int));
}

/** Get use_linear_tree.
 *
 * @param chunk The chunk.
 *
 * @return use_linear_tree.
 */
unsigned int *
spamm_chunk_get_use_linear_tree (spamm_chunk_t *chunk)
{
  return (unsigned int*) ((intptr_t) chunk + 2*sizeof(unsigned int));
}

/** Get the address of the N array.
 *
 * @param chunk The chunk.
 *
 * @return The N array.
 */
unsigned int *
spamm_chunk_get_N (spamm_chunk_t *chunk)
{
  void **chunk_pointer = (void*) ((intptr_t) chunk + 4*sizeof(unsigned int));
  return (unsigned int*) ((intptr_t) chunk + (intptr_t) chunk_pointer[0]);
}

/** Get the address of the N_lower array.
 *
 * @param chunk The chunk.
 *
 * @return The N_lower array.
 */
unsigned int *
spamm_chunk_get_N_lower (spamm_chunk_t *chunk)
{
  void **chunk_pointer = (void*) ((intptr_t) chunk + 4*sizeof(unsigned int));
  return (unsigned int*) ((intptr_t) chunk + (intptr_t) chunk_pointer[1]);
}

/** Get the address of the N_upper array.
 *
 * @param chunk The chunk.
 *
 * @return The N_upper array.
 */
unsigned int *
spamm_chunk_get_N_upper (spamm_chunk_t *chunk)
{
  void **chunk_pointer = (void*) ((intptr_t) chunk + 4*sizeof(unsigned int));
  return (unsigned int*) ((intptr_t) chunk + (intptr_t) chunk_pointer[2]);
}

/** Get the size of the contiguous matrix block.
 *
 * @param chunk The chunk.
 *
 * @return The size of the contiguous matrix.
 */
unsigned int
spamm_chunk_get_N_contiguous (spamm_chunk_t *chunk)
{
  unsigned int *N_lower;
  unsigned int *N_upper;

  N_lower = spamm_chunk_get_N_lower(chunk);
  N_upper = spamm_chunk_get_N_upper(chunk);

  return N_upper[0]-N_lower[0];
}

/** Get the address of matrix block.
 *
 * @param chunk The chunk.
 *
 * @return The address of the matrix.
 */
float *
spamm_chunk_get_matrix (spamm_chunk_t *chunk)
{
  void **chunk_pointer = (void*) ((intptr_t) chunk + 4*sizeof(unsigned int));
  return (float*) ((intptr_t) chunk + (intptr_t) chunk_pointer[3]);
}

/** Get the address of dilated matrix block.
 *
 * @param chunk The chunk.
 *
 * @return The address of the dilated matrix.
 */
float *
spamm_chunk_get_matrix_dilated (spamm_chunk_t *chunk)
{
  void **chunk_pointer = (void*) ((intptr_t) chunk + 4*sizeof(unsigned int));
  return (float*) ((intptr_t) chunk + (intptr_t) chunk_pointer[4]);
}

/** Get the total number of norm entries across all tiers stored in the SpAMM
 * chunk.
 *
 * @param chunk The chunk.
 *
 * @return The total number of entries.
 */
unsigned int
spamm_chunk_get_number_norm_entries (spamm_chunk_t *chunk)
{
  unsigned int result = 0;
  unsigned int number_dimensions;
  unsigned int tier;

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  for (tier = 0; tier < *spamm_chunk_get_number_tiers(chunk); tier++)
  {
    result += ipow(ipow(2, number_dimensions), tier);
  }

  return result;
}

/** Get a pointer to the norm arrays.
 *
 * @param chunk The chunk.
 *
 * @return The pointer to the norm arrays.
 */
float *
spamm_chunk_get_norm (spamm_chunk_t *chunk)
{
  void **chunk_pointer = (void*) ((intptr_t) chunk + 4*sizeof(unsigned int));
  return (float*) ((intptr_t) chunk + (intptr_t) chunk_pointer[5]);
}

/** Calculate the starting address of the norm arrays inside a SpAMM chunk.
 * The norm arrays start at tier == chunk_tier, with one entry, and
 * generally have pow(pow(2, number_dimensions), tier-chunk_tier)
 * entries.
 *
 * @param tier The tier.
 * @param chunk The chunk.
 *
 * @return A pointer to the start of the norm chunk at this tier.
 */
float *
spamm_chunk_get_tier_norm (const unsigned int tier,
    spamm_chunk_t *chunk)
{
  float *norm;
  unsigned int number_dimensions;
  unsigned int offset = 0;

  unsigned int i;

  norm = spamm_chunk_get_norm(chunk);
  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);

  for (i = 0; i < tier; i++)
  {
    offset += ipow(ipow(2, number_dimensions), i);
  }

  return &norm[offset];
}

/** Get a pointer to the square of the norm arrays.
 *
 * @param chunk The chunk.
 *
 * @return The pointer to the norm2 arrays.
 */
float *
spamm_chunk_get_norm2 (spamm_chunk_t *chunk)
{
  void **chunk_pointer = (void*) ((intptr_t) chunk + 4*sizeof(unsigned int));
  return (float*) ((intptr_t) chunk + (intptr_t) chunk_pointer[6]);
}

/** Calculate the starting address of the norm2 arrays inside a SpAMM chunk.
 * The norm arrays start at tier == chunk_tier, with one entry, and
 * generally have pow(pow(2, number_dimensions), tier-chunk_tier)
 * entries.
 *
 * @param tier The tier.
 * @param chunk The chunk.
 *
 * @return A pointer to the start of the norm chunk at this tier.
 */
float *
spamm_chunk_get_tier_norm2 (const unsigned int tier,
    spamm_chunk_t *chunk)
{
  float *norm2;
  unsigned int number_dimensions;
  unsigned int offset = 0;

  unsigned int i;

  norm2 = spamm_chunk_get_norm2(chunk);
  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);

  for (i = 0; i < tier; i++)
  {
    offset += ipow(ipow(2, number_dimensions), i);
  }

  return &norm2[offset];
}

/** Calculate a linear offset into a SpAMM chunk matrix. The matrix data is
 * stored in a matrix of width N_contiguous, of submatrices of width N_block.
 * Storage order is column-major (for reasons of compatibilty with Fortran).
 *
 * @param number_dimensions The number of dimensions.
 * @param N_block The block size.
 * @param N_lower The bounding box upper index array.
 * @param i The index array.
 *
 * @return The linear offset into the matrix data.
 */
unsigned int
spamm_chunk_matrix_index (const unsigned int number_dimensions,
    const unsigned int N_block,
    const unsigned int *const N_lower,
    const unsigned int *const i)
{
  unsigned int offset;
  unsigned int block_offset = 0;
  unsigned int *i_temp;

  int dim;

  i_temp = calloc(number_dimensions, sizeof(unsigned int));
  for (dim = 0; dim < number_dimensions; dim++)
  {
    i_temp[dim] = (i[dim]-N_lower[dim])/N_block;
  }
  offset = ipow(N_block, number_dimensions)*spamm_index_linear(number_dimensions, i_temp);
  free(i_temp);

  for (dim = number_dimensions-1; dim >= 1; dim--)
  {
    block_offset = N_block*(block_offset+(i[dim]-N_lower[dim])%N_block);
  }
  block_offset += (i[0]-N_lower[0])%N_block;
  offset += block_offset;

  return offset;
}

/** Get the size of a SpAMM data chunk.
 *
 * See the documentation for spamm_new_chunk() for a detailed description of
 * the contents of a SpAMM data chunk.
 *
 * @param number_dimensions [in] The number of dimensions.
 * @param use_linear_tree [in] Whether to use the linear tree code or not.
 * @param number_tiers [out] The number of tiers stored in this chunk.
 * @param N [in] The array of the original, unpadded matrix size.
 * @param N_lower [in] The array of the bounding box.
 * @param N_upper [in] The array of the bounding box.
 * @param N_lower_pointer [out] Pointer to N_lower[].
 * @param N_upper_pointer [out] Pointer to N_upper[].
 * @param A_pointer [out] Pointer to A.
 * @param A_dilated_pointer [out] Pointer to A_dilated.
 * @param norm_pointer [out] Pointer to norm[].
 * @param norm2_pointer [out] Pointer to norm2[].
 *
 * @return The size in bytes of the chunk.
 */
size_t
spamm_chunk_get_size (const unsigned int number_dimensions,
    const short use_linear_tree,
    unsigned int *number_tiers,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    unsigned int **N_pointer,
    unsigned int **N_lower_pointer,
    unsigned int **N_upper_pointer,
    float **A_pointer,
    float **A_dilated_pointer,
    float **norm_pointer,
    float **norm2_pointer)
{
  unsigned int tier;
  unsigned int N_contiguous;
  int dim;

  double tier_temp;

  size_t size = 0;

  /* Calculate sizes needed. */
  for (dim = 0, N_contiguous = 0; dim < number_dimensions; dim++)
  {
    if (N_contiguous == 0)
    {
      N_contiguous = N_upper[dim]-N_lower[dim];
    }

    if (N_upper[dim]-N_lower[dim] != N_contiguous)
    {
      SPAMM_FATAL("only square shaped regions are supported (N_upper[%u]-N_lower[%u] = %u, N_contiguous = %u)\n",
          dim, dim, N_upper[dim]-N_lower[dim], N_contiguous);
    }
  }

  /* Calculate the number of tiers stored here. */
  if (use_linear_tree)
  {
    if (N_contiguous < SPAMM_N_KERNEL)
    {
      SPAMM_FATAL("logic error, N_contiguous (%u) has to be at least %u\n", N_contiguous, SPAMM_N_KERNEL);
    }
    tier_temp = log(N_contiguous/SPAMM_N_KERNEL)/log(2);

    if (tier_temp-round(tier_temp) < 1e-10)
    {
      *number_tiers = (unsigned int) round(tier_temp);
    }

    else
    {
      *number_tiers = (unsigned int) ceil(tier_temp);
    }

    /* Add another tier for norms at SPAMM_N_BLOCK. */
    *number_tiers += 1;

    if (SPAMM_N_KERNEL*ipow(2, *number_tiers-1) != N_contiguous)
    {
      SPAMM_FATAL("logic error, number_tiers = %u, N_contiguous = %u, %u*2^%u = %u\n",
          *number_tiers, N_contiguous, SPAMM_N_KERNEL, *number_tiers-1,
          SPAMM_N_KERNEL*ipow(2, *number_tiers-1));
    }
  }

  else
  {
    *number_tiers = 1;
  }

  /* Pointers to fields. */
  size += sizeof(unsigned int);  /* number_dimensions */
  size += sizeof(unsigned int);  /* number_tiers */
  size += sizeof(unsigned int);  /* use_linear_tree */
  size += sizeof(unsigned int);  /* Padding. */

  size += sizeof(unsigned int*); /* N_pointer */
  size += sizeof(unsigned int*); /* N_lower_pointer */
  size += sizeof(unsigned int*); /* N_upper_pointer */
  size += sizeof(float*);        /* A_pointer */
  size += sizeof(float*);        /* A_dilated_pointer */
  size += sizeof(float*);        /* norm_pointer */
  size += sizeof(float*);        /* norm2_pointer */

  /* Fields. */
  *N_pointer       = (unsigned int*) size; size += number_dimensions*sizeof(unsigned int); /* N[number_dimensions] */
  *N_lower_pointer = (unsigned int*) size; size += number_dimensions*sizeof(unsigned int); /* N_lower[number_dimensions] */
  *N_upper_pointer = (unsigned int*) size; size += number_dimensions*sizeof(unsigned int); /* N_upper[number_dimensions] */

  /* Pad. */
  size = spamm_chunk_pad(size, SPAMM_ALIGNMENT);

  /* Matrix data. */
  *A_pointer = (float*) size; size += ipow(N_contiguous, number_dimensions)*sizeof(float); /* A[ipow(N_contiguous, number_dimensions)] */

  /* Pad. */
  size = spamm_chunk_pad(size, SPAMM_ALIGNMENT);

  /* Dilated matrix data. */
  *A_dilated_pointer = (float*) size; size += 4*ipow(N_contiguous, number_dimensions)*sizeof(float); /* A_dilated[4*N_contiguous, number_dimensions)] */

  /* Norm. */
  *norm_pointer = (float*) size;

  /* Add up all tiers. */
  for (tier = 0; tier < *number_tiers; tier++)
  {
    size += ipow(ipow(2, tier), number_dimensions)*sizeof(float);
  }

  /* Squared norm. */
  *norm2_pointer = (float*) size;

  /* Add up all tiers. */
  for (tier = 0; tier < *number_tiers; tier++)
  {
    size += ipow(ipow(2, tier), number_dimensions)*sizeof(float);
  }

  return size;
}

/** Print some information on a SpAMM chunk.
 *
 * @param chunk The chunk.
 */
void
spamm_chunk_print (spamm_chunk_t *chunk)
{
  int dim;
  unsigned int number_dimensions;
  unsigned int N_contiguous;
  unsigned int *N;
  unsigned int *N_lower;
  unsigned int *N_upper;
  float *A;
  float *norm;

  printf("chunk:\n");
  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  printf("number_dimensions: %u\n", number_dimensions);
  N_contiguous = spamm_chunk_get_N_contiguous(chunk);
  printf("N_contiguous: %u\n", N_contiguous);
  N = spamm_chunk_get_N(chunk);
  printf("N = [");
  for (dim = 0; dim < number_dimensions; dim++)
  {
    printf(" %u", N[dim]);
    if (dim+1 < number_dimensions)
    {
      printf(",");
    }
  }
  printf(" ]\n");
  N_lower = spamm_chunk_get_N_lower(chunk);
  N_upper = spamm_chunk_get_N_upper(chunk);
  for (dim = 0; dim < number_dimensions; dim++)
  {
    printf("N[%u] = [ %u, %u ]\n", dim, N_lower[dim], N_upper[dim]);
  }
  A = spamm_chunk_get_matrix(chunk);
  spamm_print_dense(N_contiguous, N_contiguous, column_major, A);
  norm = spamm_chunk_get_norm(chunk);
}

/** Allocate a SpAMM data chunk.
 *
 * The chunk contains the following data fields. In order to guarantee this
 * layout we allocate a larger chunk of memory and then manage the data inside
 * of it ourselves.  In order to simplify access to the fields, we start the
 * chunk with a pointer array that points to the field variables.
 *
 * \code
 * struct spamm_chunk_t
 * {
 *   unsigned int *number_dimensions_pointer;
 *   unsigned int *N_block_pointer;
 *   unsigned int *N_lower_pointer;
 *   unsigned int *N_upper_pointer;
 *   float        *A_pointer;
 *   float        *A_dilated_pointer;
 *   float        *norm_pointer;
 *   float        *norm2_pointer;
 *
 *   unsigned int number_dimensions;
 *   unsigned int N_block;
 *   unsigned int N_lower[number_dimensions];
 *   unsigned int N_upper[number_dimensions];
 *
 *   spamm_float_t *A;
 *
 *   spamm_float_t *A_dilated;
 *
 *   spamm_float_t norm[];
 *   spamm_float_t norm2[];
 * };
 * \endcode
 *
 * @param number_dimensions The number of dimensions.
 * @param use_linear_tree Whether to use the linear code for the chunk or not.
 * @param N The size of original matrix (unpadded).
 * @param N_lower The lower bounds of the bounding box.
 * @param N_lower The upper bounds of the bounding box.
 *
 * @return A pointer to the newly allocated chunk.
 */
spamm_chunk_t *
spamm_new_chunk (const unsigned int number_dimensions,
    const short use_linear_tree,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper)
{
  void **pointer_pointer;
  unsigned int *int_pointer;

  unsigned int number_tiers;

  unsigned int *N_pointer;
  unsigned int *N_lower_pointer;
  unsigned int *N_upper_pointer;
  float *A_pointer;
  float *A_dilated_pointer;
  float *norm_pointer;
  float *norm2_pointer;

  int dim;

  spamm_chunk_t *chunk;

  chunk = spamm_allocate(spamm_chunk_get_size(number_dimensions,
        use_linear_tree, &number_tiers, N, N_lower, N_upper, &N_pointer,
        &N_lower_pointer, &N_upper_pointer, &A_pointer, &A_dilated_pointer,
        &norm_pointer, &norm2_pointer), 1);

  int_pointer = chunk;
  pointer_pointer = (void**) ((intptr_t) chunk + 4*sizeof(unsigned int));

  int_pointer[0] = number_dimensions;
  int_pointer[1] = number_tiers;
  int_pointer[2] = use_linear_tree;

  pointer_pointer[0] = (void*) N_pointer;
  pointer_pointer[1] = (void*) N_lower_pointer;
  pointer_pointer[2] = (void*) N_upper_pointer;
  pointer_pointer[3] = (void*) A_pointer;
  pointer_pointer[4] = (void*) A_dilated_pointer;
  pointer_pointer[5] = (void*) norm_pointer;
  pointer_pointer[6] = (void*) norm2_pointer;

  /* Store bounding box. */
  N_pointer       = (unsigned int*) ((intptr_t) chunk + (intptr_t) N_pointer);
  N_lower_pointer = (unsigned int*) ((intptr_t) chunk + (intptr_t) N_lower_pointer);
  N_upper_pointer = (unsigned int*) ((intptr_t) chunk + (intptr_t) N_upper_pointer);

  for (dim = 0; dim < number_dimensions; dim++)
  {
    N_pointer[dim] = N[dim];
    N_lower_pointer[dim] = N_lower[dim];
    N_upper_pointer[dim] = N_upper[dim];
  }

  return chunk;
}

/** Delete a chunk. This function simply calls free().
 *
 * @param chunk The chunk to delete.
 */
void
spamm_delete_chunk (spamm_chunk_t **chunk)
{
  if (chunk == NULL) { return; }

  if (*chunk != NULL)
  {
    free(*chunk);
  }
  *chunk = NULL;
}
