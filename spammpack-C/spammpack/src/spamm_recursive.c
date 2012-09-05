#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** @brief Get an element from a recursive matrix.
 *
 * If the matri is NULL, this function returns 0.
 *
 * @param i The row index.
 * @param j The column index.
 * @param A The matrix.
 *
 * @return The matrix element Aij.
 */
double
spamm_recursive_get (const unsigned int i, const unsigned int j, const struct spamm_recursive_t *A)
{
  struct spamm_recursive_node_t **node = NULL;

  if (A == NULL)
  {
    return 0;
  }

  node = (struct spamm_recursive_node_t**) &(A->root);

  while (1)
  {
    if (*node == NULL)
    {
      return 0;
    }

    if ((*node)->M_upper-(*node)->M_lower == A->N_contiguous)
    {
      /* Get the matrix element. */
      if ((*node)->data == NULL) { return 0.0; }
      else
      {
        return (*node)->data[spamm_index_column_major(i-(*node)->M_lower, j-(*node)->N_lower, A->N_contiguous, A->N_contiguous)];
      }
    }

    else
    {
      if (i < (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
          j < (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[0]);
      }

      else if (i <  (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
               j >= (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[1]);
      }

      else if (i >= (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
               j <  (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[2]);
      }

      else if (i >= (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
               j >= (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[3]);
      }
    }
  }
}

/** @brief Print a recursive SpAMM matrix.
 *
 * @param A The matrix.
 */
void
spamm_recursive_print (const struct spamm_recursive_t *A)
{
  int i, j;

  if (A == NULL)
  {
    printf("A = (null)\n");
  }

  else
  {
    for (i = 0; i < A->M; i++) {
      for (j = 0; j < A->N; j++)
      {
        printf(" % 1.2e", spamm_recursive_get(i, j, A));
      }
      printf("\n");
    }
  }
}
