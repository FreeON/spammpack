/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>

/** Copy a SpAMM chunk. \f$ A \leftarrow \beta B \f$.
 *
 * @param A Chunk A.
 * @param beta The scalar \f$ \beta \f$.
 * @param B Chunk B.
 * @param flop The flop count.
 * @param mop The memory operation count
 */
void
spamm_chunk_copy (spamm_chunk_t **A,
    const float beta,
    spamm_chunk_t *B,
    const short use_linear_tree,
    double *const flop,
    double *const mop)
{
  unsigned int number_dimensions;
  unsigned int number_tiers;
  unsigned int *N;
  unsigned int *N_lower;
  unsigned int *N_upper;

  spamm_norm_t *norm_A;
  spamm_norm_t *norm2_A;
  spamm_norm_t *norm2_B;

  float *A_matrix;
  float *B_matrix;

  float *A_matrix_dilated;

  unsigned int N_contiguous;

  unsigned int i;
  unsigned int i_norm;

  unsigned int tier;

  assert(A != NULL);
  assert(B != NULL);
  assert(flop != NULL);
  assert(mop != NULL);

  spamm_delete_chunk(A);

  number_dimensions = *spamm_chunk_get_number_dimensions(B);
  number_tiers = *spamm_chunk_get_number_tiers(B);
  N = spamm_chunk_get_N(B);
  N_lower = spamm_chunk_get_N_lower(B);
  N_upper = spamm_chunk_get_N_upper(B);

  /* Allocate memory for new chunk. */
  *A = spamm_new_chunk(number_dimensions, use_linear_tree, N, N_lower, N_upper);

  norm_A = spamm_chunk_get_norm(*A);
  norm2_A = spamm_chunk_get_norm2(*A);
  norm2_B = spamm_chunk_get_norm2(B);

  /* Update norms. */
  for(tier = 0, i_norm = 0; tier < number_tiers; tier++)
  {
    for(i = 0; i < ipow(ipow(2, number_dimensions), tier); i++)
    {
      norm2_A[i_norm] = beta*beta*norm2_B[i_norm];
      norm_A[i_norm] = sqrt(norm2_A[i_norm]);
      i_norm++;
    }
  }

  N_contiguous = spamm_chunk_get_N_contiguous(B);

  A_matrix = spamm_chunk_get_matrix(*A);
  B_matrix = spamm_chunk_get_matrix(B);

  A_matrix_dilated = spamm_chunk_get_matrix_dilated(*A);

  /* Copy matrix elements. */
  for(i = 0; i < ipow(N_contiguous, number_dimensions); i++)
  {
    A_matrix[i] = beta*B_matrix[i];

    A_matrix_dilated[4*i+0] = A_matrix[i];
    A_matrix_dilated[4*i+1] = A_matrix[i];
    A_matrix_dilated[4*i+2] = A_matrix[i];
    A_matrix_dilated[4*i+3] = A_matrix[i];
  }

  *flop += ipow(N_contiguous, number_dimensions);
  *mop += 4*ipow(N_contiguous, number_dimensions);
}

/** Copy a matrix. \f$ A \leftarrow \beta B \f$.
 *
 * @param A The matrix to copy to.
 * @param beta The scalar beta.
 * @param B The matrix to copy from.
 * @param flop The flop count.
 * @param mop The memory operation count
 */
void
spamm_recursive_copy (struct spamm_recursive_node_t *const A,
    const float beta,
    const struct spamm_recursive_node_t *const B,
    const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree,
    double *const flop,
    double *const mop)
{
  short i;

  /* We need access to a lock on A, which is why A can not be NULL. We need to
   * allocate a new node if we need to a tier above. */
  assert(A != NULL);

  if(B == NULL) { return; }

  if(tier == chunk_tier)
  {
#ifdef _OPENMP
    omp_set_lock(&A->lock);
#endif

    spamm_chunk_copy(&A->tree.chunk, beta, B->tree.chunk, use_linear_tree, flop, mop);

#ifdef _OPENMP
    omp_unset_lock(&A->lock);
#endif
  }

  else
  {
    if(B->tree.child == NULL)
    {
      return;
    }

#ifdef _OPENMP
    omp_set_lock(&A->lock);
#endif

    if(A->tree.child == NULL)
    {
      A->tree.child = calloc(ipow(2, number_dimensions), sizeof(struct spamm_recursive_node_t*));
    }

    for(i = 0; i < ipow(2, number_dimensions); i++)
    {
      if(B->tree.child[i] != NULL && A->tree.child[i] == NULL)
      {
        A->tree.child[i] = spamm_recursive_new_node();
      }
    }

#ifdef _OPENMP
    omp_unset_lock(&A->lock);
#endif

    for(i = 0; i < ipow(2, number_dimensions); i++)
    {
      if(B->tree.child[i] != NULL)
      {
#pragma omp task untied
        spamm_recursive_copy(A->tree.child[i], beta, B->tree.child[i],
            number_dimensions, tier+1, chunk_tier, use_linear_tree, flop,
            mop);
      }
    }
#pragma omp taskwait
  }

#ifdef _OPENMP
    omp_set_lock(&A->lock);
#endif

  A->norm2 = beta*beta*B->norm2;
  A->norm = sqrt(A->norm2);

#ifdef _OPENMP
    omp_unset_lock(&A->lock);
#endif
}

/** Copy a matrix. \f$ A \leftarrow \beta B \f$.
 *
 * @param A The matrix to copy to.
 * @param beta The scalar beta.
 * @param B The matrix to copy from.
 * @param flop The flop count.
 * @param mop The memory operation count
 */
void
spamm_copy (struct spamm_matrix_t **A,
    const float beta,
    const struct spamm_matrix_t *const B,
    double *const flop,
    double *const mop)
{
  assert(A != NULL);

  if(*A == B) { return; }

  spamm_delete(A);

  if(B == NULL) { return; }

  *A = spamm_new(B->number_dimensions, B->N, B->chunk_tier, B->use_linear_tree);
  (*A)->recursive_tree = spamm_recursive_new_node();

#pragma omp parallel
  {
#pragma omp single
    {
#pragma omp task untied
      spamm_recursive_copy((*A)->recursive_tree, beta, B->recursive_tree,
          B->number_dimensions, 0, B->chunk_tier, B->use_linear_tree, flop,
          mop);
    }
  }
}
