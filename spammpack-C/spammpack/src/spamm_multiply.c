#include "spamm.h"
#include "spamm_types_private.h"

#include <stdio.h>
#include <stdlib.h>

/** Multiply two matrices, i.e. \f$ C = \alpha A \times B + \beta C\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 * @param timer The timer to use.
 * @param sgemm The external sgemm function to use.
 */
void
spamm_multiply (const float tolerance,
    const float alpha, struct spamm_recursive_t *A, struct spamm_recursive_t *B,
    const float beta, struct spamm_recursive_t *C,
    struct spamm_timer_t *timer,
    void (*sgemm) (char *, char *, int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *),
    unsigned int *number_products)
{
  if (A == NULL)
  {
    printf("[%s:%i] A is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (B == NULL)
  {
    printf("[%s:%i] B is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (C == NULL)
  {
    printf("[%s:%i] C is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (A->blocksize != B->blocksize)
  {
    printf("[%s:%i] A->blocksize != B->blocksize\n", __FILE__, __LINE__);
    exit(1);
  }

  if (A->blocksize != C->blocksize)
  {
    printf("[%s:%i] A->blocksize != C->blocksize\n", __FILE__, __LINE__);
    exit(1);
  }

  /* Multiply C by beta. */
  spamm_recursive_multiply_scalar(beta, C->root);

  /* Multiply A and B. */
  if (A->root != NULL && B->root != NULL && C->root == NULL)
  {
    C->root = spamm_recursive_new_node(0, C->blocksize);
    C->root->M_lower = 0;
    C->root->M_upper = A->root->M_upper;
    C->root->N_lower = 0;
    C->root->N_upper = A->root->N_upper;
  }

  spamm_recursive_multiply_matrix(tolerance, alpha, A->root, B->root, &(C->root), timer, sgemm, number_products);
}
