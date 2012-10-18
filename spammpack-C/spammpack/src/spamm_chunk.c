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
  return (uint32_t*) chunk;
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
  return (uint32_t*) ((void*) spamm_chunk_get_number_dimensions(chunk) + sizeof(uint32_t));
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
  return (uint32_t*) ((void*) spamm_chunk_get_N_contiguous(chunk) + sizeof(uint32_t));
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
  return (uint32_t*) ((void*) spamm_chunk_get_N_lower(chunk)
      + (*spamm_chunk_get_number_dimensions(chunk))*sizeof(uint32_t));
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
  return (float*) (spamm_chunk_align((size_t) ((void*) spamm_chunk_get_N_upper(chunk)+(*spamm_chunk_get_number_dimensions(chunk))*sizeof(uint32_t)), SPAMM_ALIGNMENT));
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
  return (float*) (spamm_chunk_align((size_t) ((void*) spamm_chunk_get_matrix(chunk)
          +ipow(*spamm_chunk_get_N_contiguous(chunk), *spamm_chunk_get_number_dimensions(chunk))*sizeof(float)), SPAMM_ALIGNMENT));
}

/** Get the size of a SpAMM data chunk.
 *
 * @param number_dimensions The number of dimensions.
 * @param N_contiguous The size of the contigous matrix in this chunk.
 *
 * @return The size in bytes of the chunk.
 */
size_t
spamm_get_chunk_size (const unsigned int number_dimensions,
    const unsigned int N_contiguous)
{
  size_t size;

  /* Dimensions and bounding box. */
  size = sizeof(uint32_t)                /* number_dimensions */
    +sizeof(uint32_t)                    /* N_contiguous. */
    +number_dimensions*sizeof(uint32_t)  /* N_lower[number_dimensions] */
    +number_dimensions*sizeof(uint32_t); /* N_upper[number_dimensions] */

  /* Pad. */
  size = spamm_chunk_align(size, SPAMM_ALIGNMENT);

  /* Matrix data. */
  size += ipow(N_contiguous, number_dimensions)*sizeof(float); /* A[ipow(N_contiguous, number_dimensions)] */

  /* Pad. */
  size = spamm_chunk_align(size, SPAMM_ALIGNMENT);

  /* Dilated matrix data. */
  size += 4*ipow(N_contiguous, number_dimensions)*sizeof(float); /* A_dilated[4*N_contiguous, number_dimensions)] */

  return size;
}
