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
    const unsigned int contiguous_tier,
    const unsigned int depth,
    const unsigned int N_block,
    const unsigned int linear_index_A,
    const unsigned int linear_index_B,
    const unsigned int linear_index_C,
    sgemm_func sgemm)
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

  float *matrix_A;
  float *matrix_B;
  float *matrix_C;

  number_dimensions_A = *spamm_chunk_get_number_dimensions(chunk_A);
  number_dimensions_B = *spamm_chunk_get_number_dimensions(chunk_B);
  number_dimensions_C = *spamm_chunk_get_number_dimensions(chunk_C);

#ifdef PRINT_DEBUG
  N_lower = spamm_chunk_get_N_lower(chunk_A);
  N_upper = spamm_chunk_get_N_upper(chunk_A);
  printf("chunk_A: {");
  for (dim = 0; dim < number_dimensions_A; dim++)
  {
    printf(" [ %u, %u )", N_lower[dim], N_upper[dim]);
    if (dim+1 < number_dimensions_A) { printf(","); }
  }
  printf(" }\n");

  N_lower = spamm_chunk_get_N_lower(chunk_B);
  N_upper = spamm_chunk_get_N_upper(chunk_B);
  printf("chunk_B: {");
  for (dim = 0; dim < number_dimensions_B; dim++)
  {
    printf(" [ %u, %u )", N_lower[dim], N_upper[dim]);
    if (dim+1 < number_dimensions_B) { printf(","); }
  }
  printf(" }\n");

  N_lower = spamm_chunk_get_N_lower(chunk_C);
  N_upper = spamm_chunk_get_N_upper(chunk_C);
  printf("chunk_C: {");
  for (dim = 0; dim < number_dimensions_C; dim++)
  {
    printf(" [ %u, %u )", N_lower[dim], N_upper[dim]);
    if (dim+1 < number_dimensions_C) { printf(","); }
  }
  printf(" }\n");
#endif

  norm_A = spamm_chunk_get_tier_norm(tier-contiguous_tier, chunk_A);
  norm_B = spamm_chunk_get_tier_norm(tier-contiguous_tier, chunk_B);

  norm_C = spamm_chunk_get_tier_norm(tier-contiguous_tier, chunk_C);
  norm2_C = spamm_chunk_get_tier_norm2(tier-contiguous_tier, chunk_C);

  if (norm_A[linear_index_A]*norm_B[linear_index_B] > tolerance)
  {
    if (number_dimensions_A == 2 &&
        number_dimensions_B == 2 &&
        number_dimensions_C == 2)
    {
      if (tier == depth)
      {
        if (sgemm)
        {
          SPAMM_FATAL("FIXME\n");
        }

        else
        {
          matrix_A = spamm_chunk_get_matrix(chunk_A);
          matrix_B = spamm_chunk_get_matrix(chunk_B);
          matrix_C = spamm_chunk_get_matrix(chunk_C);

          matrix_A += ipow(N_block, number_dimensions_A)*linear_index_A;
          matrix_B += ipow(N_block, number_dimensions_B)*linear_index_B;
          matrix_C += ipow(N_block, number_dimensions_C)*linear_index_C;

          norm2_C[linear_index_C] = 0;
          for (i = 0; i < N_block; i++) {
            for (j = 0; j < N_block; j++) {
              for (k = 0; k < N_block; k++)
              {
                matrix_C[i+N_block*j] += alpha*matrix_A[i+N_block*k]*matrix_B[k+N_block*j];
                norm2_C[linear_index_C] += ipow(matrix_C[i+N_block*j], 2);
              }
            }
          }
          norm_C[linear_index_C] = sqrt(norm2_C[linear_index_C]);
        }
      }

      else
      {
        norm2_C[linear_index_C] = 0;
        for (i = 0; i < 2; i++) {
          for (j = 0; j < 2; j++)
          {
            new_linear_index_C = linear_index_C << number_dimensions_C;
            new_linear_index_C |= i | (j << 1);

            for (k = 0; k < 2; k++)
            {
              new_linear_index_A = linear_index_A << number_dimensions_A;
              new_linear_index_B = linear_index_B << number_dimensions_B;

              new_linear_index_A |= i | (k << 1);
              new_linear_index_B |= k | (j << 1);

              norm2_C[linear_index_C] += spamm_chunk_multiply(tolerance,
                  alpha, chunk_A, chunk_B, chunk_C, tier+1, contiguous_tier,
                  depth, N_block, new_linear_index_A, new_linear_index_B,
                  new_linear_index_C, sgemm);
            }
          }
        }
        norm_C[linear_index_C] = sqrt(norm2_C[linear_index_C]);
      }
    }

    else
    {
      SPAMM_FATAL("not implemented\n");
    }
  }

  return norm2_C[linear_index_C];
}

void
spamm_chunk_add (const float alpha,
    spamm_chunk_t **A,
    const float beta,
    spamm_chunk_t *B)
{
  SPAMM_FATAL("FIXME\n");
}

void
spamm_chunk_copy (spamm_chunk_t **A,
    const float beta,
    spamm_chunk_t *B)
{
  unsigned int *number_dimensions = spamm_chunk_get_number_dimensions(B);
  unsigned int *N_block = spamm_chunk_get_N_block(B);
  unsigned int *N = spamm_chunk_get_N(B);
  unsigned int *N_lower = spamm_chunk_get_N_lower(B);
  unsigned int *N_upper = spamm_chunk_get_N_upper(B);

  float *norm_A;
  float *norm_B;
  float *norm2_A;
  float *norm2_B;

  float *A_matrix;
  float *B_matrix;

  unsigned int number_tiers;
  unsigned int N_contiguous;

  unsigned int i;

  spamm_delete_chunk(A);

  *A = spamm_new_chunk(*number_dimensions, *N_block, N, N_lower, N_upper);

  norm_A = spamm_chunk_get_norm(*A);
  norm_B = spamm_chunk_get_norm(B);
  norm2_A = spamm_chunk_get_norm2(*A);
  norm2_B = spamm_chunk_get_norm2(B);

  number_tiers = spamm_chunk_get_number_tiers(B);

  for (i = 0; i < ipow(*number_dimensions, number_tiers); i++)
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
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const unsigned int depth,
    const unsigned int linear_index,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    spamm_chunk_t *chunk)
{
  int dim;

  unsigned int number_dimensions;

  float *norm;
  float *norm2;

  unsigned int N_block;
  float *A;

  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  unsigned int new_linear_index = linear_index;

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  norm = spamm_chunk_get_tier_norm(tier-contiguous_tier, chunk);
  norm2 = spamm_chunk_get_tier_norm2(tier-contiguous_tier, chunk);

  /* Update norm. */
  norm2[linear_index] += Aij*Aij;
  norm[linear_index]   = sqrt(norm2[linear_index]);

  if (tier == depth)
  {
    N_block = *spamm_chunk_get_N_block(chunk);
    A = spamm_chunk_get_matrix(chunk);

    A[ipow(N_block, number_dimensions)*linear_index
      +spamm_index_column_major_2(number_dimensions, N_block, N_lower, i)] = Aij;
  }

  else
  {
    new_linear_index <<= number_dimensions;

    new_N_lower = calloc(number_dimensions, sizeof(unsigned int));
    new_N_upper = calloc(number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < number_dimensions; dim++)
    {
      if (i[dim] < N_lower[dim]+(N_upper[dim]-N_lower[dim])/2)
      {
        new_N_lower[dim] = N_lower[dim];
        new_N_upper[dim] = N_lower[dim]+(N_upper[dim]-N_lower[dim])/2;
      }

      else
      {
        new_N_lower[dim] = N_lower[dim]+(N_upper[dim]-N_lower[dim])/2;
        new_N_upper[dim] = N_upper[dim];
        new_linear_index |= (1 << dim);
      }
    }

    spamm_chunk_set(i, Aij, tier+1, contiguous_tier, depth, new_linear_index, new_N_lower, new_N_upper, chunk);

    free(new_N_lower);
    free(new_N_upper);
  }
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
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const unsigned int depth,
    const unsigned int linear_index,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    spamm_chunk_t *chunk)
{
  float Aij = 0;

  int dim;

  unsigned int number_dimensions;

  unsigned int N_block;
  float *A;

  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  unsigned int new_linear_index = linear_index;

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);

  if (tier == depth)
  {
    N_block = *spamm_chunk_get_N_block(chunk);
    A = spamm_chunk_get_matrix(chunk);

    Aij = A[ipow(N_block, number_dimensions)*linear_index
      +spamm_index_column_major_2(number_dimensions, N_block, N_lower, i)];
  }

  else
  {
    new_linear_index <<= number_dimensions;

    new_N_lower = calloc(number_dimensions, sizeof(unsigned int));
    new_N_upper = calloc(number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < number_dimensions; dim++)
    {
      if (i[dim] < N_lower[dim]+(N_upper[dim]-N_lower[dim])/2)
      {
        new_N_lower[dim] = N_lower[dim];
        new_N_upper[dim] = N_lower[dim]+(N_upper[dim]-N_lower[dim])/2;
      }

      else
      {
        new_N_lower[dim] = N_lower[dim]+(N_upper[dim]-N_lower[dim])/2;
        new_N_upper[dim] = N_upper[dim];
        new_linear_index |= (1 << dim);
      }
    }

    Aij = spamm_chunk_get(i, tier+1, contiguous_tier, depth, new_linear_index, new_N_lower, new_N_upper, chunk);

    free(new_N_lower);
    free(new_N_upper);
  }

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
 * @return N_block.
 */
unsigned int *
spamm_chunk_get_N_block (spamm_chunk_t *chunk)
{
  return (unsigned int*) ((intptr_t) chunk + sizeof(unsigned int));
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
  void **chunk_pointer = (void*) ((intptr_t) chunk + 2*sizeof(unsigned int));
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
  void **chunk_pointer = (void*) ((intptr_t) chunk + 2*sizeof(unsigned int));
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
  void **chunk_pointer = (void*) ((intptr_t) chunk + 2*sizeof(unsigned int));
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
  void **chunk_pointer = (void*) ((intptr_t) chunk + 2*sizeof(unsigned int));
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
  void **chunk_pointer = (void*) ((intptr_t) chunk + 2*sizeof(unsigned int));
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
  for (tier = 0; tier < spamm_chunk_get_number_tiers(chunk); tier++)
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
  void **chunk_pointer = (void*) ((intptr_t) chunk + 2*sizeof(unsigned int));
  return (float*) ((intptr_t) chunk + (intptr_t) chunk_pointer[5]);
}

/** Calculate the starting address of the norm arrays inside a SpAMM chunk.
 * The norm arrays start at tier == contiguous_tier, with one entry, and
 * generally have pow(pow(2, number_dimensions), tier-contiguous_tier)
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
  void **chunk_pointer = (void*) ((intptr_t) chunk + 2*sizeof(unsigned int));
  return (float*) ((intptr_t) chunk + (intptr_t) chunk_pointer[6]);
}

/** Calculate the starting address of the norm2 arrays inside a SpAMM chunk.
 * The norm arrays start at tier == contiguous_tier, with one entry, and
 * generally have pow(pow(2, number_dimensions), tier-contiguous_tier)
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

/** Get the number of tiers stored in the SpAMM chunk.
 *
 * @param chunk The chunk.
 *
 * @return The number of tiers.
 */
unsigned int
spamm_chunk_get_number_tiers (spamm_chunk_t *chunk)
{
  unsigned int N_contiguous;
  unsigned int N_block;
  unsigned int number_tiers;
  unsigned int N;

  N_contiguous = spamm_chunk_get_N_contiguous(chunk);
  N_block = *spamm_chunk_get_N_block(chunk);
  for (N = N_contiguous, number_tiers = 0; N >= N_block; N /= 2)
  {
    number_tiers++;
  }
  return number_tiers;
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
 * @param N_block [in] The size at which the SpAMM condition is applied.
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
    const unsigned int N_block,
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
  unsigned int N_temp;
  unsigned int N_contiguous;
  int dim;

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

  /* Pointers to fields. */
  size += sizeof(unsigned int);  /* number_dimensions */
  size += sizeof(unsigned int);  /* N_block */
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

  for (N_temp = N_contiguous; N_temp >= N_block; N_temp /= 2)
  {
    size += ipow(N_contiguous/N_temp, number_dimensions)*sizeof(float);
  }

  /* Squared norm. */
  *norm2_pointer = (float*) size;

  for (N_temp = N_contiguous; N_temp >= N_block; N_temp /= 2)
  {
    size += ipow(N_contiguous/N_temp, number_dimensions)*sizeof(float);
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
 * @param N_block The size of matrices to which the SpAMM condition is
 * applied.
 * @param N The size of original matrix (unpadded).
 * @param N_lower The lower bounds of the bounding box.
 * @param N_lower The upper bounds of the bounding box.
 *
 * @return A pointer to the newly allocated chunk.
 */
spamm_chunk_t *
spamm_new_chunk (const unsigned int number_dimensions,
    const unsigned int N_block,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper)
{
  void **chunk_pointer;
  unsigned int *int_pointer;

  unsigned int *N_pointer;
  unsigned int *N_lower_pointer;
  unsigned int *N_upper_pointer;
  float *A_pointer;
  float *A_dilated_pointer;
  float *norm_pointer;
  float *norm2_pointer;

  int dim;

  spamm_chunk_t *chunk;

  chunk = spamm_allocate(spamm_chunk_get_size(number_dimensions, N_block, N,
        N_lower, N_upper, &N_pointer, &N_lower_pointer, &N_upper_pointer,
        &A_pointer, &A_dilated_pointer, &norm_pointer, &norm2_pointer), 1);

  chunk_pointer = (void**) ((intptr_t) chunk + 2*sizeof(unsigned int));
  int_pointer = chunk;

  int_pointer[0] = number_dimensions;
  int_pointer[1] = N_block;

  chunk_pointer[0] = (void*) N_pointer;
  chunk_pointer[1] = (void*) N_lower_pointer;
  chunk_pointer[2] = (void*) N_upper_pointer;
  chunk_pointer[3] = (void*) A_pointer;
  chunk_pointer[4] = (void*) A_dilated_pointer;
  chunk_pointer[5] = (void*) norm_pointer;
  chunk_pointer[6] = (void*) norm2_pointer;

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
