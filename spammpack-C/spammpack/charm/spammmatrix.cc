#include "spammmatrix.h"

/** Initialize matrix with random elements.
 *
 * @param N The size of the matrix.
 */
SpAMMMatrix::SpAMMMatrix (const int N)
{
  CkPrintf("Initializing %dx%d matrix\n", N, N);
  this->N = N;
}

SpAMMMatrix::SpAMMMatrix (CkMigrateMessage *msg) {}

#include "spammmatrix.def.h"
