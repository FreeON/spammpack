/** @file */

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
 * @param B The SpAMM_Matrix object identifier.
 */
void
spamm_convert_dense_to_spamm_interface (const int *const M,
    const int *const N,
    const float *const A_dense,
    int *B)
{
  struct spamm_hashed_t *B_spamm;

  B_spamm = spamm_interface_get_spamm_object(*B);
  if (B_spamm != NULL)
  {
    printf("[%s:%i] object already exists\n", __FILE__, __LINE__);
    exit(1);
  }
  B_spamm = spamm_convert_dense_to_spamm(*M, *N, column_major, A_dense, row_major);

  *B = spamm_interface_add_spamm_object(B_spamm);
}

void
spamm_convert_dense_to_spamm_interface_ (const int *const M,
    const int *const N,
    const float *const A_dense,
    int *B)
{
  spamm_convert_dense_to_spamm_interface(M, N, A_dense, B);
}

/** Create a new spamm matrix.
 *
 * @param A The matrix.
 */
void
spamm_new_interface (int *A)
{
}

void
spamm_new_interface_ (int *A)
{
  spamm_new_interface(A);
}
