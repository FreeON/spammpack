/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Add two spamm matrices. @f$ A \leftarrow \alpha A + \beta B @f$.
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param A The matrix A.
 * @param beta The factor @f$ \beta @f$.
 * @param B The matrix B.
 * @param number_dimensions The number of dimensions.
 */
void
spamm_recursive_add (const float alpha,
    struct spamm_recursive_node_t **A,
    const float beta,
    struct spamm_recursive_node_t **B,
    const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree)
{
  unsigned int i;

  unsigned int N_contiguous;

  float *A_matrix;
  float *B_matrix;

  /* There is nothing to do here. */
  if((*A) == NULL && (*B) == NULL)
  {
    return;
  }

  if(tier == chunk_tier)
  {
    if((*A) == NULL && (*B) != NULL)
    {
      SPAMM_FATAL("FIXME\n");
    }

    else if((*A) != NULL && (*B) == NULL)
    {
      SPAMM_FATAL("FIXME\n");
    }

    else
    {
      spamm_chunk_add(alpha, &(*A)->tree.chunk, beta, (*B)->tree.chunk);
    }
  }

  else
  {
    if((*A) == NULL && (*B) != NULL)
    {
      /* Copy B node to A. */
      spamm_recursive_copy(&(*A), beta, (*B), number_dimensions, tier,
          chunk_tier, use_linear_tree);
    }

    else if((*A) != NULL && (*B) == NULL)
    {
      /* Multiply A by alpha. */
      spamm_recursive_multiply_scalar(alpha, *A, number_dimensions, tier,
          chunk_tier, use_linear_tree);
    }

    else
    {
      /* Recurse. */
      for(i = 0; i < ipow(2, number_dimensions); i++)
      {
        spamm_recursive_add(alpha, &(*A)->tree.child[i], beta,
            &(*B)->tree.child[i], number_dimensions, tier+1, chunk_tier,
            use_linear_tree);
      }
    }
  }
}

/** Add two spamm matrices. @f$ A \leftarrow \alpha A + \beta B @f$.
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param A The matrix A.
 * @param beta The factor @f$ \beta @f$.
 * @param B The matrix B.
 */
void
spamm_add (const float alpha,
    struct spamm_matrix_t *A,
    const float beta,
    struct spamm_matrix_t *B)
{
  if(A->chunk_tier == 0)
  {
    spamm_chunk_add(alpha, &A->tree.chunk, beta, B->tree.chunk);
  }

  else
  {
    spamm_recursive_add(alpha, &A->tree.recursive_tree, beta,
        &B->tree.recursive_tree, A->number_dimensions, 0, A->chunk_tier,
        A->use_linear_tree);
  }
}
