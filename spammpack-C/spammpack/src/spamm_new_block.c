#include "config.h"
#include "spamm.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

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
