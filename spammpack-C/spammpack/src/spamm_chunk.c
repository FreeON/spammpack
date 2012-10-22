#include "spamm.h"

#include <stdio.h>

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
uint32_t *
spamm_chunk_get_number_dimensions (spamm_chunk_t *chunk)
{
  spamm_chunk_t **chunk_pointer = chunk;
  return (uint32_t*) ((uint64_t) chunk + (uint64_t) chunk_pointer[0]);
}

/** Get N_contiguous.
 *
 * @param chunk The chunk.
 *
 * @return N_contiguous.
 */
uint32_t *
spamm_chunk_get_N_contiguous (spamm_chunk_t *chunk)
{
  spamm_chunk_t **chunk_pointer = chunk;
  return (uint32_t*) ((uint64_t) chunk + (uint64_t) chunk_pointer[1]);
}

/** Get the address of the N_lower array.
 *
 * @param chunk The chunk.
 *
 * @return The N_lower array.
 */
uint32_t *
spamm_chunk_get_N_lower (spamm_chunk_t *chunk)
{
  spamm_chunk_t **chunk_pointer = chunk;
  return (uint32_t*) ((uint64_t) chunk + (uint64_t) chunk_pointer[2]);
}

/** Get the address of the N_upper array.
 *
 * @param chunk The chunk.
 *
 * @return The N_upper array.
 */
uint32_t *
spamm_chunk_get_N_upper (spamm_chunk_t *chunk)
{
  spamm_chunk_t **chunk_pointer = chunk;
  return (uint32_t*) ((uint64_t) chunk + (uint64_t) chunk_pointer[3]);
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
  return (float*) ((uint64_t) chunk + (uint64_t) chunk_pointer[4]);
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
  return (float*) ((uint64_t) chunk + (uint64_t) chunk_pointer[5]);
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
    uint32_t **number_dimension_pointer,
    uint32_t **N_contiguous_pointer,
    uint32_t **N_lower_pointer,
    uint32_t **N_upper_pointer,
    float **A_pointer,
    float **A_dilated_pointer,
    float **norm_pointer,
    float **norm2_pointer)
{
  size_t size = 0;

  /* Pointers to fields. */
  size += sizeof(uint32_t*); /* number_dimensions_pointer */
  size += sizeof(uint32_t*); /* N_contiguous_pointer */
  size += sizeof(uint32_t*); /* N_lower_pointer */
  size += sizeof(uint32_t*); /* N_upper_pointer */
  size += sizeof(float*);    /* A_pointer */
  size += sizeof(float*);    /* A_dilated_pointer */
  size += sizeof(float*);    /* norm_pointer */
  size += sizeof(float*);    /* norm2_pointer */

  /* Fields. */
  *number_dimension_pointer = (uint32_t*) size; size += sizeof(uint32_t);  /* number_dimensions */
  *N_contiguous_pointer = (uint32_t*) size;     size += sizeof(uint32_t);  /* N_contiguous. */
  *N_lower_pointer = (uint32_t*) size;          size += number_dimensions*sizeof(uint32_t); /* N_lower[number_dimensions] */
  *N_upper_pointer = (uint32_t*) size;          size += number_dimensions*sizeof(uint32_t); /* N_upper[number_dimensions] */

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
  uint32_t number_dimensions;
  uint32_t *N_lower;
  uint32_t *N_upper;

  printf("chunk:\n");
  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  printf("number_dimensions: %u\n", number_dimensions);
  printf("N_contiguous: %u\n", *spamm_chunk_get_N_contiguous(chunk));
  N_lower = spamm_chunk_get_N_lower(chunk);
  N_upper = spamm_chunk_get_N_upper(chunk);
  for (dim = 0; dim < number_dimensions; dim++)
  {
    printf("N[%u] = [ %u, %u ]\n", dim, N_lower[dim], N_upper[dim]);
  }
}
