/** @file
 *
 * The implementation of the Reduction class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "reduction.h"
#include "logger.h"
#include "messages.h"

Reduction::Reduction (int N, CProxy_ReductionData a)
{
  this->N = N;
  this->index = thisIndex.x*N*N + thisIndex.y*N + thisIndex.z;
  this->callCount = 0;
  this->result = NULL;
  this->a = a;
  INFO("R(%d,%d,%d) constructor\n", thisIndex.x, thisIndex.y, thisIndex.z);
}

Reduction::Reduction (CkMigrateMessage *msg)
{
}

Reduction::~Reduction ()
{
  INFO("R(%d,%d,%d) destructor\n", thisIndex.x, thisIndex.y, thisIndex.z);
  delete[] result;
}

void Reduction::reduce (CkCallback &cb)
{
  INFO("R(%d,%d,%d) reduce\n", thisIndex.x, thisIndex.y, thisIndex.z);
  if(callCount > 0)
  {
    ABORT("R(%d,%d,%d) this chare has been called before\n", thisIndex.x, thisIndex.y, thisIndex.z);
  }
  callCount++;

  result = new double[N*N];
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      for(int k = 0; k < N; k++)
      {
        DoubleMsg * aik = a(i, k).get();
        DoubleMsg * akj = a(k, j).get();
        result[i*N+j] += aik->x * akj->x;
        delete aik;
        delete akj;
      }
    }
  }

  contribute(cb);
  migrateMe(0);
}

void Reduction::pup (PUP::er &p)
{
  CBase_Reduction::pup(p);
  p|N;
  p|index;
  p|callCount;
  p|a;

  if(p.isUnpacking())
  {
    INFO("R(%d,%d,%d) pup: unpacking\n", thisIndex.x, thisIndex.y, thisIndex.z);
  }
  else
  {
    if(p.isSizing())
    {
      INFO("R(%d,%d,%d) pup: sizing\n", thisIndex.x, thisIndex.y, thisIndex.z);
    }
    else
    {
      INFO("R(%d,%d,%d) pup: packing\n", thisIndex.x, thisIndex.y, thisIndex.z);
    }
  }

  int numberElements = (result == NULL ? 0 : N*N);
  p|numberElements;

  if(numberElements == 0)
  {
    if(p.isUnpacking()) { result = NULL; }
  }
  else
  {
    if(p.isUnpacking()) { result = new double[numberElements]; }
    PUParray(p, result, numberElements);
  }
}

#include "reduction.def.h"
