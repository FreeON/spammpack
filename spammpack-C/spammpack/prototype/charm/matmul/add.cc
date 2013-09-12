/** @file
 *
 * The implementation of the Add class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "add.h"

/** The constructor. The Add object is set up so that a subsequent call to
 * Add::add() will result in
 *
 * @f[ A \leftarrow A + \alpha B @f]
 *
 * @param A Matrix A.
 * @param B Matrix B.
 */
Add::Add (CProxy_Matrix A, CProxy_Matrix B)
{
  this->A = A;
  this->B = B;
}

/** Add the two matrices.
 *
 * @param alpha The factor @f$ \alpha @f$.
 */
void Add::add (double alpha)
{
}

#include "add.def.h"
