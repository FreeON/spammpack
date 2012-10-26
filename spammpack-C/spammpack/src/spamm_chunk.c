#include "config.h"
#include "spamm.h"

#include <math.h>
#include <stdio.h>
#include <stdint.h>

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

/** Get the number of tiers stored in the SpAMM chunk.
 *
 * @param chunk
 *
 * @return The number of tiers.
 */
unsigned int
spamm_chunk_get_number_tiers (spamm_chunk_t *chunk)
{
  unsigned int N_contiguous;
  unsigned int number_tiers;
  unsigned int N;

  N_contiguous = spamm_chunk_get_N_contiguous(chunk);
  for (N = N_contiguous, number_tiers = 0; N >= SPAMM_N_BLOCK; N /= 2)
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
 * @param N_lower The bounding box lower index array.
 * @param N_upper The bounding box upper index array.
 * @param i The index array.
 *
 * @return The linear offset into the matrix data.
 */
unsigned int
spamm_chunk_matrix_index (const unsigned int number_dimensions,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    const unsigned int *const i)
{
  unsigned int offset = 0;
  unsigned int dim;

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

  for (N_temp = N_contiguous; N_temp >= SPAMM_N_BLOCK; N_temp /= 2)
  {
    size += ipow(N_contiguous/N_temp, 2)*sizeof(float);
  }

  /* Squared norm. */
  *norm2_pointer = (float*) size;

  for (N_temp = N_contiguous; N_temp >= SPAMM_N_BLOCK; N_temp /= 2)
  {
    size += ipow(N_contiguous/N_temp, 2)*sizeof(float);
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
  unsigned int i, j;
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
