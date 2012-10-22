#include "spamm.h"

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
spamm_chunk_align (const size_t address,
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
  spamm_chunk_t **chunk_pointer = chunk;
  return (unsigned int*) ((intptr_t) chunk + (intptr_t) chunk_pointer[0]);
}

/** Get N_contiguous.
 *
 * @param chunk The chunk.
 *
 * @return N_contiguous.
 */
unsigned int *
spamm_chunk_get_N_contiguous (spamm_chunk_t *chunk)
{
  spamm_chunk_t **chunk_pointer = chunk;
  return (unsigned int*) ((intptr_t) chunk + (intptr_t) chunk_pointer[1]);
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
  spamm_chunk_t **chunk_pointer = chunk;
  return (unsigned int*) ((intptr_t) chunk + (intptr_t) chunk_pointer[2]);
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
  spamm_chunk_t **chunk_pointer = chunk;
  return (unsigned int*) ((intptr_t) chunk + (intptr_t) chunk_pointer[3]);
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
  spamm_chunk_t **chunk_pointer = chunk;
  return (float*) ((intptr_t) chunk + (intptr_t) chunk_pointer[4]);
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
  spamm_chunk_t **chunk_pointer = chunk;
  return (float*) ((intptr_t) chunk + (intptr_t) chunk_pointer[5]);
}

/** Get the size of a SpAMM data chunk.
 *
 * See the documentation for spamm_new_chunk() for a detailed description of
 * the contents of a SpAMM data chunk.
 *
 * @param number_dimensions The number of dimensions.
 * @param N_contiguous The size of the contigous matrix in this chunk.
 * @param number_dimension_pointer [out] Pointer to number_dimension.
 * @param N_contiguous_pointer [out] Pointer to N_contiguous.
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
    const unsigned int N_contiguous,
    unsigned int **number_dimension_pointer,
    unsigned int **N_contiguous_pointer,
    unsigned int **N_lower_pointer,
    unsigned int **N_upper_pointer,
    float **A_pointer,
    float **A_dilated_pointer,
    float **norm_pointer,
    float **norm2_pointer)
{
  size_t size = 0;

  /* Pointers to fields. */
  size += sizeof(unsigned int*); /* number_dimensions_pointer */
  size += sizeof(unsigned int*); /* N_contiguous_pointer */
  size += sizeof(unsigned int*); /* N_lower_pointer */
  size += sizeof(unsigned int*); /* N_upper_pointer */
  size += sizeof(float*);    /* A_pointer */
  size += sizeof(float*);    /* A_dilated_pointer */
  size += sizeof(float*);    /* norm_pointer */
  size += sizeof(float*);    /* norm2_pointer */

  /* Fields. */
  *number_dimension_pointer = (unsigned int*) size; size += sizeof(unsigned int);  /* number_dimensions */
  *N_contiguous_pointer = (unsigned int*) size;     size += sizeof(unsigned int);  /* N_contiguous. */
  *N_lower_pointer = (unsigned int*) size;          size += number_dimensions*sizeof(unsigned int); /* N_lower[number_dimensions] */
  *N_upper_pointer = (unsigned int*) size;          size += number_dimensions*sizeof(unsigned int); /* N_upper[number_dimensions] */

  /* Pad. */
  size = spamm_chunk_align(size, SPAMM_ALIGNMENT);

  /* Matrix data. */
  *A_pointer = (float*) size; size += ipow(N_contiguous, number_dimensions)*sizeof(float); /* A[ipow(N_contiguous, number_dimensions)] */

  /* Pad. */
  size = spamm_chunk_align(size, SPAMM_ALIGNMENT);

  /* Dilated matrix data. */
  *A_dilated_pointer = (float*) size; size += 4*ipow(N_contiguous, number_dimensions)*sizeof(float); /* A_dilated[4*N_contiguous, number_dimensions)] */

  /* Norm. */

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
  unsigned int *N_lower;
  unsigned int *N_upper;
  float *A;

  printf("chunk:\n");
  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  printf("number_dimensions: %u\n", number_dimensions);
  N_contiguous = *spamm_chunk_get_N_contiguous(chunk);
  printf("N_contiguous: %u\n", N_contiguous);
  N_lower = spamm_chunk_get_N_lower(chunk);
  N_upper = spamm_chunk_get_N_upper(chunk);
  for (dim = 0; dim < number_dimensions; dim++)
  {
    printf("N[%u] = [ %u, %u ]\n", dim, N_lower[dim], N_upper[dim]);
  }
  A = spamm_chunk_get_matrix(chunk);
  spamm_print_dense(N_contiguous, N_contiguous, column_major, A);
}

/** Allocate a SpAMM data chunk.
 *
 * The chunk contains the following data fields. Since compilers don't
 * guarantee certain alignments and padding unless coerced, we manage the
 * position and size of the fields in the chunk ourselves. In order to
 * simplify access to the fields, we start the chunk with a pointer array that
 * points to the field variables.
 *
 * \code
 * struct spamm_chunk_t
 * {
 *   unsigned int  *number_dimensions_pointer;
 *   unsigned int  *N_contiguous_pointer;
 *   unsigned int  *N_lower_pointer;
 *   unsigned int  *N_upper_pointer;
 *   spamm_float_t *A_pointer;
 *   spamm_float_t *A_dilated_pointer;
 *   spamm_float_t *norm_pointer;
 *   spamm_float_t *norm2_pointer;
 *
 *   unsigned int number_dimensions;
 *   unsigned int N_contiguous;
 *   unsigned int N_lower[number_dimensions];
 *   unsigned int N_upper[number_dimensions];
 *
 *   char padding[];
 *
 *   spamm_float_t *A;
 *
 *   char padding[];
 *
 *   spamm_float_t *A_dilated;
 *
 *   char padding[];
 *
 *   spamm_float_t norm[];
 *   spamm_float_t norm2[];
 * };
 * \endcode
 *
 * @param number_dimensions The number of dimensions.
 * @param N_contiguous The size of the contigous matrix in this chunk.
 *
 * @return A pointer to the newly allocated chunk.
 */
spamm_chunk_t *
spamm_new_chunk (const unsigned int number_dimensions,
    const unsigned int N_contiguous)
{
  void **chunk_pointer;
  unsigned int *number_dimension_pointer;
  unsigned int *N_contiguous_pointer;
  unsigned int *N_lower_pointer;
  unsigned int *N_upper_pointer;
  float *A_pointer;
  float *A_dilated_pointer;
  float *norm_pointer;
  float *norm2_pointer;

  spamm_chunk_t *chunk;

  chunk = spamm_allocate(spamm_chunk_get_size(number_dimensions, N_contiguous,
        &number_dimension_pointer, &N_contiguous_pointer,
        &N_lower_pointer, &N_upper_pointer,
        &A_pointer, &A_dilated_pointer,
        &norm_pointer, &norm2_pointer), 1);

  chunk_pointer = chunk;

  chunk_pointer[0] = (void*) number_dimension_pointer;
  chunk_pointer[1] = (void*) N_contiguous_pointer;
  chunk_pointer[2] = (void*) N_lower_pointer;
  chunk_pointer[3] = (void*) N_upper_pointer;
  chunk_pointer[4] = (void*) A_pointer;
  chunk_pointer[5] = (void*) A_dilated_pointer;
  chunk_pointer[6] = (void*) norm_pointer;
  chunk_pointer[7] = (void*) norm2_pointer;

  number_dimension_pointer = (unsigned int*) ((intptr_t) chunk + (intptr_t) number_dimension_pointer);
  *number_dimension_pointer = number_dimensions;

  N_contiguous_pointer = (unsigned int*) ((intptr_t) chunk + (intptr_t) N_contiguous_pointer);
  *N_contiguous_pointer = N_contiguous;

  return chunk;
}
