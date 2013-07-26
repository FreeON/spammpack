/** @file
 *
 * The implementation of the ReductionData class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "reductiondata.h"
#include "logger.h"
#include "messages.h"

ReductionData::ReductionData (int N)
{
  this->N = N;
  this->index = thisIndex.x*N + thisIndex.y;
  this->x = new double[1];
  this->x[0] = rand()/(double) RAND_MAX;
  INFO("RD(%d,%d) constructor\n", thisIndex.x, thisIndex.y);
}

ReductionData::ReductionData (CkMigrateMessage *msg)
{
}

ReductionData::~ReductionData ()
{
  INFO("RD(%d,%d) destructor\n", thisIndex.x, thisIndex.y);
  delete[] x;
}

DoubleMsg * ReductionData::get ()
{
  return new DoubleMsg(x[0]);
}

void ReductionData::pup (PUP::er &p)
{
  CBase_ReductionData::pup(p);
  p|N;
  p|index;

  if(p.isUnpacking())
  {
    INFO("RD(%d,%d) pup: unpacking\n", thisIndex.x, thisIndex.y);
  }
  else
  {
    if(p.isSizing())
    {
      INFO("RD(%d,%d) pup: sizing\n", thisIndex.x, thisIndex.y);
    }
    else
    {
      INFO("RD(%d,%d) pup: packing\n", thisIndex.x, thisIndex.y);
    }
  }

  int numberElements = (x == NULL ? 0 : N*N);
  p|numberElements;

  if(numberElements == 0)
  {
    if(p.isUnpacking()) { x = NULL; }
  }
  else
  {
    if(p.isUnpacking()) { x = new double[numberElements]; }
    PUParray(p, x, numberElements);
  }
}

#include "reductiondata.def.h"
