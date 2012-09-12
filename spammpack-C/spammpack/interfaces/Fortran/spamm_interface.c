/** @file */

#include "config.h"
#include "spamm.h"

#include <stdio.h>
#include <stdlib.h>

/** A linked list type. */
struct spamm_interface_object_t
{
  /** The link to the next element. */
  struct spamm_interface_object_t *previous;

  /** The link to the previous element. */
  struct spamm_interface_object_t *next;

  /** The index of this object. */
  unsigned int ID;

  /** The pointer to the spamm_hashed_t object. */
  struct spamm_hashed_t *object;
};

/** The spamm object buffer for the interface. We store the real objects here
 * so we can translate the object IDs we give out through the interface to the
 * actual object.
 */
struct spamm_interface_object_t *spamm_matrix_object = NULL;

/** Get a spamm_hashed_t object from the interface list.
 *
 * @param ID The ID of the object.
 *
 * @return A pointer to the spamm_hashed_t object or NULL if the ID does not
 * exist in the store.
 */
struct spamm_hashed_t *
spamm_interface_get_spamm_object (const int ID)
{
  struct spamm_interface_object_t *element;

  if (spamm_matrix_object == NULL)
  {
    return NULL;
  }

  for (element = spamm_matrix_object; element != NULL; element = element->next)
  {
    if (element->ID == ID)
    {
      return element->object;
    }
  }

  return NULL;
}

/** Pop a spamm_hashed_t object from the interface list. If the object exists,
 * it will be removed from the store.
 *
 * @param ID The ID of the object.
 *
 * @return A pointer to the spamm_hashed_t object or NULL if the ID does not
 * exist in the store.
 */
struct spamm_hashed_t *
spamm_interface_pop_spamm_object (const int ID)
{
  struct spamm_hashed_t *object;
  struct spamm_interface_object_t *element;

  if (ID == 0)
  {
    return NULL;
  }

  for (element = spamm_matrix_object; element != NULL; element = element->next)
  {
    if (element->ID == ID)
    {
      if (element == spamm_matrix_object)
      {
        /* Move the beginning of linked list. */
        spamm_matrix_object = element->next;
      }

      /* Relink list. */
      if (element->previous != NULL)
      {
        element->previous->next = element->next;
      }
      if (element->next != NULL)
      {
        element->next->previous = element->previous;
      }

      /* Get object... */
      object = element->object;
      free(element);

      /* and return it. */
      return object;
    }
  }

  return NULL;
}

/** Allocates a new linked list element for the matrix object store.
 *
 * @return A pointer to a newly allocated element.
 */
struct spamm_interface_object_t *
spamm_interface_new_object ()
{
  return calloc(1, sizeof(struct spamm_interface_object_t));
}

/** Add a matrix object to the interface list.
 *
 * @param A The pointer to the spamm_hashed_t object to add to the store.
 *
 * @return The newly assigned ID.
 */
unsigned int
spamm_interface_add_spamm_object (struct spamm_hashed_t *A)
{
  struct spamm_interface_object_t *element;

  if (spamm_matrix_object == NULL)
  {
    spamm_matrix_object = spamm_interface_new_object();
    spamm_matrix_object->object = A;

    /* We'll start the ID with 1, so that an ID of zero can not be found in
     * the store. */
    spamm_matrix_object->ID = 1;
    return spamm_matrix_object->ID;
  }

  else
  {
    for (element = spamm_matrix_object; element->next != NULL; element = element->next) {}
    element->next = spamm_interface_new_object();
    element->next->ID = element->ID+1;
    element->next->previous = element;
    element->next->object = A;
    return element->next->ID;
  }
}

/** Convert dense matrix to spamm matrix.
 *
 * @param A_dense The dense matrix.
 * @param A The SpAMM_Matrix object identifier.
 */
void
FC_FUNC(spamm_convert_dense_to_spamm_interface, SPAMM_CONVERT_DENSE_TO_SPAMM_INTERFACE) (const int *const M,
    const int *const N,
    float *A_dense,
    int *A)
{
  struct spamm_hashed_t *A_spamm;

  A_spamm = spamm_interface_get_spamm_object(*A);
  if (A_spamm != NULL)
  {
    printf("[%s:%i] object already exists\n", __FILE__, __LINE__);
    exit(1);
  }
  A_spamm = spamm_convert_dense_to_spamm(*M, *N, column_major, A_dense, row_major);

  *A = spamm_interface_add_spamm_object(A_spamm);
}

/** Create a new spamm matrix.
 *
 * @param A The matrix.
 */
void
FC_FUNC(spamm_new_interface, SPAMM_NEW_INTERFACE) (const int *const M, const int *const N, int *A)
{
  struct spamm_hashed_t *A_spamm;

  A_spamm = spamm_interface_pop_spamm_object(*A);
  if (A_spamm != NULL)
  {
    printf("[%s:%i] deleting already existing matrix\n", __FILE__, __LINE__);
    spamm_hashed_delete(&A_spamm);
  }
  //A_spamm = spamm_hashed_new(*M, *N, row_major);

  *A = spamm_interface_add_spamm_object(A_spamm);
}

/** Delete a spamm matrix.
 *
 * @param A The matrix.
 */
void
FC_FUNC(spamm_delete_interface, SPAMM_DELETE_INTERFACE) (int *A)
{
  struct spamm_hashed_t *A_spamm;

  A_spamm = spamm_interface_pop_spamm_object(*A);
  if (A_spamm == NULL)
  {
    return;
  }
  spamm_hashed_delete(&A_spamm);
}

/** Multiply two spamm matrices.
 *
 * @param A The matrix A.
 * @param B The matrix A.
 * @param C The matrix A.
 * @param tolerance The SpAMM tolerance.
 */
void
FC_FUNC(spamm_multiply_spamm_spamm_interface, SPAMM_MULTIPLY_SPAMM_SPAMM_INTERFACE) (int *A,
    int *B, int *C, const float *const tolerance)
{
  struct spamm_hashed_t *A_spamm;
  struct spamm_hashed_t *B_spamm;
  struct spamm_hashed_t *C_spamm;

  struct spamm_timer_t *timer;

  A_spamm = spamm_interface_get_spamm_object(*A);
  B_spamm = spamm_interface_get_spamm_object(*B);
  C_spamm = spamm_interface_get_spamm_object(*C);

  if (C_spamm == NULL)
  {
    //printf("[%s:%i] creating new C matrix\n", __FILE__, __LINE__);
    //C_spamm = spamm_hashed_new(spamm_get_number_of_rows(A_spamm), spamm_get_number_of_columns(B_spamm), row_major);
    *C = spamm_interface_add_spamm_object(C_spamm);
  }

  timer = spamm_timer_new();
  spamm_hashed_multiply(*tolerance, 1.0, A_spamm, B_spamm, 1.0, C_spamm, timer, kernel_standard_SSE);
  spamm_timer_delete(&timer);
}

/** Add two spamm matrices. @f$ A \leftarrow \alpha A + \beta B @f$.
 *
 * @param A The matrix A.
 * @param B The matrix B.
 * @param alpha The factor @f$ \alpha @f$.
 * @param beta The factor @f$ \beta @f$.
 */
void
FC_FUNC(spamm_add_spamm_spamm_interface, SPAMM_ADD_SPAMM_SPAMM_INTERFACE) (int *A,
    int *B, float *alpha, float *beta)
{
  struct spamm_hashed_t *A_spamm;
  struct spamm_hashed_t *B_spamm;

  A_spamm = spamm_interface_get_spamm_object(*A);
  B_spamm = spamm_interface_get_spamm_object(*B);

  spamm_hashed_add(*alpha, A_spamm, *beta, B_spamm);
}

/** Get Frobenius norm of matrix.
 *
 * @param A The matrix.
 *
 * @return The Frobenius norm.
 */
float
FC_FUNC(spamm_get_norm_interface, SPAMM_GET_NORM_INTERFACE) (int *A)
{
  struct spamm_hashed_t *A_spamm;

  A_spamm = spamm_interface_get_spamm_object(*A);
  return spamm_get_norm(A_spamm);
}

/** Print a matrix.
 *
 * @param A The matrix.
 */
void
FC_FUNC(spamm_print_interface, SPAMM_PRINT_INTERFACE) (int *A)
{
  struct spamm_hashed_t *A_spamm;

  A_spamm = spamm_interface_get_spamm_object(*A);
  spamm_print(A_spamm);
}

/** Set the number of threads (within OpenMP).
 *
 * @param num_threads The number of threads to use for subsequent
 * calculations.
 */
void
FC_FUNC(spamm_set_num_threads, SPAMM_SET_NUM_THREADS) (int *num_threads)
{
}
