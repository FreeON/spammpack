#include "config.h"
#include "spamm.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>

/** Get the total number of norm entries stored in a SpAMM chunk.
 *
 * @param number_tiers The number of tiers stored in the chunk.
 * @param number_dimensions The number of dimensions in the chunk.
 *
 * @return The total number of norms.
 */
unsigned int
spamm_chunk_get_total_number_norms (const unsigned int number_tiers,
    const unsigned int number_dimensions)
{
  unsigned int tier;
  unsigned int number_norms;

  for(tier = 1, number_norms = 1; tier < number_tiers; tier++)
  {
    number_norms += ipow(ipow(2, number_dimensions), tier);
  }

  return number_norms;
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

  if(new_address != address)
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
spamm_chunk_get_number_dimensions (const spamm_chunk_t *const chunk)
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
spamm_chunk_get_number_tiers (const spamm_chunk_t *const chunk)
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
spamm_chunk_get_use_linear_tree (const spamm_chunk_t *const chunk)
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
spamm_chunk_get_N (const spamm_chunk_t *const chunk)
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
spamm_chunk_get_N_lower (const spamm_chunk_t *const chunk)
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
spamm_chunk_get_N_upper (const spamm_chunk_t *const chunk)
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
spamm_chunk_get_N_contiguous (const spamm_chunk_t *const chunk)
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
spamm_chunk_get_matrix (const spamm_chunk_t *const chunk)
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
spamm_chunk_get_matrix_dilated (const spamm_chunk_t *const chunk)
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
spamm_chunk_get_number_norm_entries (const spamm_chunk_t *const chunk)
{
  unsigned int result = 0;
  unsigned int number_dimensions;
  unsigned int tier;

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  for(tier = 0; tier < *spamm_chunk_get_number_tiers(chunk); tier++)
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
spamm_norm_t *
spamm_chunk_get_norm (const spamm_chunk_t *const chunk)
{
  void **chunk_pointer = (void*) ((intptr_t) chunk + 4*sizeof(unsigned int));
  return (spamm_norm_t*) ((intptr_t) chunk + (intptr_t) chunk_pointer[5]);
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
spamm_norm_t *
spamm_chunk_get_tier_norm (const unsigned int tier,
    const spamm_chunk_t *const chunk)
{
  spamm_norm_t *norm;
  unsigned int number_dimensions;
  unsigned int offset = 0;

  unsigned int i;

  norm = spamm_chunk_get_norm(chunk);
  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);

  for(i = 0; i < tier; i++)
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
spamm_norm_t *
spamm_chunk_get_norm2 (const spamm_chunk_t *const chunk)
{
  void **chunk_pointer = (void*) ((intptr_t) chunk + 4*sizeof(unsigned int));
  return (spamm_norm_t*) ((intptr_t) chunk + (intptr_t) chunk_pointer[6]);
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
spamm_norm_t *
spamm_chunk_get_tier_norm2 (const unsigned int tier,
    spamm_chunk_t *chunk)
{
  spamm_norm_t *norm2;
  unsigned int number_dimensions;
  unsigned int offset = 0;

  unsigned int i;

  norm2 = spamm_chunk_get_norm2(chunk);
  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);

  for(i = 0; i < tier; i++)
  {
    offset += ipow(ipow(2, number_dimensions), i);
  }

  return &norm2[offset];
}

/** Calculate a linear offset into a SpAMM chunk matrix. When using a linear
 * tree, the matrix data is stored in a matrix of width N_contiguous, in
 * submatrices of width N_block.  The blocks are stored in Z-curve order, and
 * the matrix elements inside the blocks are stored in column-major order for
 * reasons of compatibilty with Fortran. When not using a linear tree, the
 * matrix elements are stored in a N_contiguous x N_contiguous matrix in
 * column-major storage order.
 *
 * @param number_dimensions The number of dimensions.
 * @param use_linear_tree Whether we are using the linear tree.
 * @param N_lower The bounding box lower index array.
 * @param N_upper The bounding box upper index array.
 * @param i The index array.
 *
 * @return The linear offset into the matrix data.
 */
unsigned int
spamm_chunk_matrix_index (const unsigned int number_dimensions,
    const short use_linear_tree,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    const unsigned int *const i)
{
  unsigned int offset;
  unsigned int block_offset;
  unsigned int *i_temp;

  int dim;

  i_temp = calloc(number_dimensions, sizeof(unsigned int));

  if(use_linear_tree)
  {
    for(dim = 0; dim < number_dimensions; dim++)
    {
      i_temp[dim] = (i[dim]-N_lower[dim])/SPAMM_N_BLOCK;
    }
    offset = ipow(SPAMM_N_BLOCK, number_dimensions)*spamm_index_linear(number_dimensions, i_temp);
    free(i_temp);

    for(dim = number_dimensions-1, block_offset = 0; dim >= 1; dim--)
    {
      block_offset = SPAMM_N_BLOCK*(block_offset+(i[dim]-N_lower[dim])%SPAMM_N_BLOCK);
    }
    block_offset += (i[0]-N_lower[0])%SPAMM_N_BLOCK;
    offset += block_offset;
  }

  else
  {
    i_temp = calloc(number_dimensions, sizeof(unsigned int));
    for(dim = 0; dim < number_dimensions; dim++)
    {
      i_temp[dim] = i[dim]-N_lower[dim];
    }
    offset = spamm_index_column_major_2(number_dimensions, N_upper[0]-N_lower[0], i_temp);
  }

  free(i_temp);

  return offset;
}

/** Get the size of a SpAMM data chunk.
 *
 * See the documentation for spamm_new_chunk() for a detailed description of
 * the contents of a SpAMM data chunk.
 *
 * @param number_dimensions [in] The number of dimensions.
 * @param use_linear_tree [in] Whether to use the linear tree code or not.
 * @param number_tiers [out] The number of tiers stored in this chunk. If
 * use_linear_tree then the chunk matrix size is >= SPAMM_N_KERNEL, and the
 * number of tiers is down to SPAMM_N_BLOCK.
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
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    unsigned int **N_pointer,
    unsigned int **N_lower_pointer,
    unsigned int **N_upper_pointer,
    float **A_pointer,
    float **A_dilated_pointer,
    spamm_norm_t **norm_pointer,
    spamm_norm_t **norm2_pointer)
{
  unsigned int N_contiguous;
  int dim;

  double tier_temp;

  size_t size = 0;

  /* Calculate sizes needed. */
  for(dim = 0, N_contiguous = 0; dim < number_dimensions; dim++)
  {
    if(N_contiguous == 0)
    {
      N_contiguous = N_upper[dim]-N_lower[dim];
    }

    if(N_upper[dim]-N_lower[dim] != N_contiguous)
    {
      SPAMM_FATAL("only square shaped regions are supported (N_upper[%u]-N_lower[%u] = %u, N_contiguous = %u)\n",
          dim, dim, N_upper[dim]-N_lower[dim], N_contiguous);
    }
  }

  /* Calculate the number of tiers stored here. */
  if(use_linear_tree)
  {
    if(N_contiguous < SPAMM_N_KERNEL)
    {
      SPAMM_FATAL("logic error, N_contiguous (%u) has to be at least %u\n", N_contiguous, SPAMM_N_KERNEL);
    }
    tier_temp = log(N_contiguous/SPAMM_N_KERNEL)/log(2)+1;

    if(tier_temp-round(tier_temp) < 1e-10)
    {
      *number_tiers = (unsigned int) round(tier_temp);
    }

    else
    {
      *number_tiers = (unsigned int) ceil(tier_temp);
    }

    /* Add more tiers for norms at SPAMM_N_BLOCK. */
    *number_tiers += SPAMM_KERNEL_DEPTH;

    if(SPAMM_N_KERNEL*ipow(2, *number_tiers-1-SPAMM_KERNEL_DEPTH) != N_contiguous)
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
  size += sizeof(spamm_norm_t*); /* norm_pointer */
  size += sizeof(spamm_norm_t*); /* norm2_pointer */

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
  *norm_pointer = (spamm_norm_t*) size;

  /* Add up all tiers. */
  size += spamm_chunk_get_total_number_norms(*number_tiers, number_dimensions)*sizeof(spamm_norm_t);

  /* Squared norm. */
  *norm2_pointer = (spamm_norm_t*) size;

  /* Add up all tiers. */
  size += spamm_chunk_get_total_number_norms(*number_tiers, number_dimensions)*sizeof(spamm_norm_t);

  return size;
}
