/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Add two SpAMM chunks. @f$ A \leftarrow \alpha A + \beta B @f$.
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param A Chunk A.
 * @param beta The factor @f$ \beta @f$.
 * @param B Chunk B.
 */
void
spamm_chunk_add (const float alpha,
    spamm_chunk_t **A,
    const float beta,
    spamm_chunk_t *B)
{
  unsigned int *number_dimensions;
  unsigned int *number_tiers;
  unsigned int *N;
  unsigned int *N_lower;
  unsigned int *N_upper;

  float *norm_A;
  float *norm_B;
  float *norm2_A;
  float *norm2_B;

  float *A_matrix;
  float *B_matrix;

  unsigned int N_contiguous;

  unsigned int i;
  unsigned int i_norm;

  unsigned int tier;

  number_dimensions = spamm_chunk_get_number_dimensions(B);
  number_tiers = spamm_chunk_get_number_tiers(B);
  N = spamm_chunk_get_N(B);
  N_lower = spamm_chunk_get_N_lower(B);
  N_upper = spamm_chunk_get_N_upper(B);

  norm_A = spamm_chunk_get_norm(*A);
  norm_B = spamm_chunk_get_norm(B);
  norm2_A = spamm_chunk_get_norm2(*A);
  norm2_B = spamm_chunk_get_norm2(B);

  /* Update norms. */
  for(tier = 0, i_norm = 0; tier < *number_tiers; tier++)
  {
    for(i = 0; i < ipow(ipow(2, *number_dimensions), tier); i++)
    {
      norm2_A[i_norm] = alpha*alpha*norm2_A[i_norm]+beta*beta*norm2_B[i_norm]+2*alpha*beta*norm_A[i_norm]*norm_B[i_norm];
      norm_A[i_norm] = sqrt(norm2_A[i_norm]);
      i_norm++;
    }
  }

  A_matrix = spamm_chunk_get_matrix(*A);
  B_matrix = spamm_chunk_get_matrix(B);

  N_contiguous = spamm_chunk_get_N_contiguous(B);

  for(i = 0; i < ipow(N_contiguous, *number_dimensions); i++)
  {
    A_matrix[i] = alpha*A_matrix[i]+beta*B_matrix[i];
  }
}

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

      (*A)->norm2 = alpha*alpha*(*A)->norm2+beta*beta*(*B)->norm2+2*alpha*beta*(*A)->norm*(*B)->norm;
      (*A)->norm = sqrt((*A)->norm2);
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

        (*A)->norm2 = alpha*alpha*(*A)->norm2+beta*beta*(*B)->norm2+2*alpha*beta*(*A)->norm*(*B)->norm;
        (*A)->norm = sqrt((*A)->norm2);
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
    struct spamm_matrix_t *const A,
    const float beta,
    const struct spamm_matrix_t *const B)
{
  SPAMM_WARN("A");
  spamm_matlab_print(A);

  SPAMM_WARN("B");
  spamm_matlab_print(B);

  if(A->chunk_tier == 0)
  {
    spamm_chunk_add(alpha, &A->tree.chunk, beta, B->tree.chunk);
  }

  else
  {
    spamm_recursive_add(alpha, &A->tree.recursive_tree, beta,
        &((struct spamm_matrix_t*) B)->tree.recursive_tree,
        A->number_dimensions, 0, A->chunk_tier, A->use_linear_tree);
  }

  SPAMM_WARN("A");
  spamm_matlab_print(A);
}
