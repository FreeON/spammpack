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
    const spamm_chunk_t *const chunk)
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

/** Calculate a linear offset into a SpAMM chunk matrix.
 *
 * When using a linear tree, the matrix data is stored in a hierarchy of
 * submatrices. At the highest level, the chunk stores the matrix elements in
 * a matrix of size N_contiguous x N_contiguous. The next tier are submatrices
 * of size SPAMM_N_KERNEL x SPAMM_N_KERNEL which are Z-curve ordered. This is
 * the kernel tier. Inside a kernel submatrix, there are submatrices of size
 * SPAMM_N_BLOCK x SPAMM_N_BLOCK, which are in row-major order. The elements
 * inside each such basic block are also stored in row-major order.
 *
 * When not using a linear tree, the matrix elements are stored in a
 * N_contiguous x N_contiguous matrix in column-major storage order (so we are
 * compatible with Fortran).
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
  unsigned int *i_mapped;
  unsigned int *i_remainder;

  int dim;

  i_mapped = calloc(number_dimensions, sizeof(unsigned int));
  i_remainder = calloc(number_dimensions, sizeof(unsigned int));

  if(use_linear_tree)
  {
#ifdef INDEX_DEBUG
    SPAMM_INFO("i        = { % 3u, % 3u }\n", i[0], i[1]);
#endif

    /* Shift indices to lower corner and divide by SPAMM_N_KERNEL. */
    for(dim = 0; dim < 2; dim++)
    {
      i_mapped[dim] = (i[dim]-N_lower[dim])/SPAMM_N_KERNEL;
      i_remainder[dim] = (i[dim]-N_lower[dim])%SPAMM_N_KERNEL;
    }
#ifdef INDEX_DEBUG
    SPAMM_INFO("i_mapped = { % 3u, % 3u } -> i_remainder = { % 3u, % 3u }\n",
        i_mapped[0], i_mapped[1], i_remainder[0], i_remainder[1]);
#endif

    /* Calculate Z-curve offset for top tier. */
    offset = SPAMM_N_KERNEL*SPAMM_N_KERNEL*spamm_index_linear(number_dimensions, i_mapped);

    for(dim = 0; dim < 2; dim++)
    {
      i_mapped[dim] = i_remainder[dim]/SPAMM_N_KERNEL_BLOCKED;
      i_remainder[dim] = i_remainder[dim]%SPAMM_N_KERNEL_BLOCKED;
    }
#ifdef INDEX_DEBUG
    SPAMM_INFO("i_mapped = { % 3u, % 3u } -> i_remainder = { % 3u, % 3u }\n",
        i_mapped[0], i_mapped[1], i_remainder[0], i_remainder[1]);
#endif

    /* Add offset within kernel block. */
    offset += SPAMM_N_BLOCK*SPAMM_N_BLOCK*(i_mapped[0]*SPAMM_N_KERNEL_BLOCKED+i_mapped[1]);

    /* Add offset within basic block. */
    offset += i_remainder[0]*SPAMM_N_BLOCK+i_remainder[1];
  }

  else
  {
    for(dim = 0; dim < number_dimensions; dim++)
    {
      i_mapped[dim] = i[dim]-N_lower[dim];
    }
    offset = spamm_index_column_major_2(number_dimensions, N_upper[0]-N_lower[0], i_mapped);
  }

  free(i_mapped);
  free(i_remainder);

  return offset;
}

/** Get the linear offset into the norm array.
 *
 * @param tier The tier.
 * @param i The array of indices of the matrix element the norm should cover.
 * @param N_lower The array of the lower bounds of the bounding box.
 * @param chunk The chunk.
 *
 * @return The offset into the tier norm.
 */
unsigned int
spamm_chunk_norm_index (const unsigned int tier,
    const unsigned int *const i,
    const spamm_chunk_t *const chunk)
{
  unsigned int offset;

  unsigned int dim;
  unsigned int N_contiguous;
  unsigned int N_block;
  unsigned int number_dimensions;
  unsigned int number_tiers;
  unsigned int *i_temp;
  unsigned int *N_lower;

  assert(chunk != NULL);
  assert(tier < *spamm_chunk_get_number_tiers(chunk));

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  number_tiers = *spamm_chunk_get_number_tiers(chunk);
  N_contiguous = spamm_chunk_get_N_contiguous(chunk);
  N_lower = spamm_chunk_get_N_lower(chunk);

  /* The norms are Z-curve ordered except for the last two tiers (in case we
   * are using the linear tree). Within the kernel part of size SPAMM_N_KERNEL
   * x SPAMM_N_KERNEL, there are SPAMM_KERNEL_DEPTH tiers in which the norms
   * are row-major ordered.
   *
   * This should eventually get cleaned up a bit. It's kind of a mess right
   * now.
   */

  /* Calculate submatrix size at tier. */
  N_block = N_contiguous/(1 << tier);

  i_temp = calloc(number_dimensions, sizeof(unsigned int));

  for(dim = 0; dim < number_dimensions; dim++)
  {
    i_temp[dim] = (i[dim]-N_lower[dim])/N_block;
  }

  if(number_tiers > SPAMM_KERNEL_DEPTH && tier+SPAMM_KERNEL_DEPTH < number_tiers)
  {
    offset = spamm_index_linear(number_dimensions, i_temp);
  }

  else
  {
    for(dim = 0; dim < number_dimensions; dim++)
    {
      i_temp[dim] = (i[dim]-N_lower[dim])/N_block;
    }
    offset = spamm_index_row_major_2(number_dimensions, 1 << tier, i_temp);
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

/** Set the norms and the dilated matrix in a chunk.
 *
 * @param chunk The chunk.
 * @param flop The flop count.
 * @param mop The memory operation count
 *
 * @return The square of the norm of the chunk.
 */
spamm_norm_t
spamm_chunk_fix (spamm_chunk_t *const chunk,
    double *const flop,
    double *const mop)
{
  float *matrix;
  float *matrix_dilated;

  spamm_norm_t *norm;
  spamm_norm_t *next_norm;
  spamm_norm_t *norm2;
  spamm_norm_t *next_norm2;

  short terminate;

  unsigned int N_contiguous;
  unsigned int number_dimensions;
  unsigned int use_linear_tree;
  unsigned int number_tiers;

  int tier;

  unsigned int dim;
  unsigned int i_stream;
  unsigned int *i;
  unsigned int i_block, j_block;
  unsigned int norm_offset;
  unsigned int matrix_offset;
  unsigned int block_offset;

  assert(chunk != NULL);

  N_contiguous = spamm_chunk_get_N_contiguous(chunk);
  use_linear_tree = *spamm_chunk_get_use_linear_tree(chunk);
  number_tiers = *spamm_chunk_get_number_tiers(chunk);

  matrix = spamm_chunk_get_matrix(chunk);
  matrix_dilated = spamm_chunk_get_matrix_dilated(chunk);

  /* Norms at deepest tier. */
  norm = spamm_chunk_get_tier_norm(number_tiers-1, chunk);
  norm2 = spamm_chunk_get_tier_norm2(number_tiers-1, chunk);

  if(use_linear_tree)
  {
    /* Norms at kernel tier. */
    if(number_tiers > SPAMM_KERNEL_DEPTH)
    {
      next_norm = spamm_chunk_get_tier_norm(number_tiers-SPAMM_KERNEL_DEPTH-1, chunk);
      next_norm2 = spamm_chunk_get_tier_norm2(number_tiers-SPAMM_KERNEL_DEPTH-1, chunk);
    }

    /* Linear trees are only used in 2 dimensions. */
    i = calloc(2, sizeof(unsigned int));

    /* Fix norms on lowest tier, and update matrix_dilated. */
    for(i_stream = 0; i_stream < ipow(N_contiguous/SPAMM_N_KERNEL, 2); i_stream++)
    {
      if(number_tiers > SPAMM_KERNEL_DEPTH)
      {
        next_norm2[i_stream] = 0.0;
      }

      /* We are at the SPAMM_N_KERNEL tier. The submatrix blocks are Z-curve
       * ordered. Below this tier we store row-major ordered submatrices of
       * size SPAMM_N_BLOCK.
       */
      for(i[0] = 0; i[0] < SPAMM_N_KERNEL_BLOCKED; i[0]++) {
        for(i[1] = 0; i[1] < SPAMM_N_KERNEL_BLOCKED; i[1]++)
        {
          norm_offset = i_stream*SPAMM_N_KERNEL_BLOCKED*SPAMM_N_KERNEL_BLOCKED /* Z-curve ordered offset. */
            +i[0]*SPAMM_N_KERNEL_BLOCKED+i[1]; /* Row-major ordered offset. */

          matrix_offset = i_stream*SPAMM_N_KERNEL*SPAMM_N_KERNEL /* Z-curve ordered offset. */
            +(i[0]*SPAMM_N_KERNEL_BLOCKED+i[1])*SPAMM_N_BLOCK*SPAMM_N_BLOCK; /* Row-major ordered offset. */

          /* Loop over matrix elements within the submatrices of size SPAMM_N_BLOCK. */
          norm2[norm_offset] = 0.0;
          for(i_block = 0; i_block < SPAMM_N_BLOCK; i_block++) {
            for(j_block = 0; j_block < SPAMM_N_BLOCK; j_block++)
            {
              block_offset = i_block*SPAMM_N_BLOCK+j_block; /* Row-major order of matrix elements. */

              /* Fix norm. */
              norm2[norm_offset] += matrix[matrix_offset+block_offset]*matrix[matrix_offset+block_offset];

              /* Fix dilated matrix. */
              matrix_dilated[4*(matrix_offset+block_offset)+0] = matrix[matrix_offset+block_offset];
              matrix_dilated[4*(matrix_offset+block_offset)+1] = matrix[matrix_offset+block_offset];
              matrix_dilated[4*(matrix_offset+block_offset)+2] = matrix[matrix_offset+block_offset];
              matrix_dilated[4*(matrix_offset+block_offset)+3] = matrix[matrix_offset+block_offset];
            }
          }
          norm[norm_offset] = sqrt(norm2[norm_offset]);

          /* Fix norms at kernel tier. */
          if(number_tiers > SPAMM_KERNEL_DEPTH)
          {
            next_norm2[i_stream] += norm2[norm_offset];
          }
        }
      }
      if(number_tiers > SPAMM_KERNEL_DEPTH)
      {
        next_norm[i_stream] = sqrt(next_norm2[i_stream]);
      }
    }

    /* Update mop count. */
    *mop += ipow(N_contiguous/SPAMM_N_KERNEL, 2)*SPAMM_N_BLOCK*SPAMM_N_BLOCK;

    /* Update norms up to the root tier of this chunk. */
    for(tier = number_tiers-SPAMM_KERNEL_DEPTH-2; tier >= 0; tier--)
    {
      norm2 = spamm_chunk_get_tier_norm2(tier+1, chunk);

      next_norm = spamm_chunk_get_tier_norm(tier, chunk);
      next_norm2 = spamm_chunk_get_tier_norm2(tier, chunk);

      for(i[0] = 0; i[0] < ipow(4, tier); i[0]++)
      {
        for(i[1] = 4*i[0], next_norm2[i[0]] = 0; i[1] < 4*(i[0]+1); i[1]++)
        {
          next_norm2[i[0]] += norm2[i[1]];
        }
        next_norm[i[0]] = sqrt(next_norm2[i[0]]);
      }
    }

    free(i);

    return next_norm2[0];
  }

  else
  {
    /* Allocate matrix index. */
    number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
    i = calloc(number_dimensions, sizeof(unsigned int));

    for(norm2[0] = 0.0, terminate = 0; !terminate; )
    {
      matrix_offset = spamm_index_column_major_2(number_dimensions, N_contiguous, i);

      /* Update norm. */
      norm2[0] += matrix[matrix_offset]*matrix[matrix_offset];

      /* Update dilated matrix. */
      matrix_dilated[4*matrix_offset+0] = matrix[matrix_offset];
      matrix_dilated[4*matrix_offset+1] = matrix[matrix_offset];
      matrix_dilated[4*matrix_offset+2] = matrix[matrix_offset];
      matrix_dilated[4*matrix_offset+3] = matrix[matrix_offset];

      /* Increment matrix index. */
      for(dim = 0; dim < number_dimensions; dim++)
      {
        i[dim]++;
        if(i[dim] >= N_contiguous)
        {
          if(dim >= number_dimensions-1)
          {
            terminate = 1;
            break;
          }
          i[dim] = 0;
        }

        else
        {
          break;
        }
      }
    }
    norm[0] = sqrt(norm2[0]);

    free(i);

    return norm2[0];
  }
}
