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
  this->x = rand()/(double) RAND_MAX;
  INFO("(%d) constructor\n", index);
}

ReductionData::ReductionData (CkMigrateMessage *msg)
{
}

ReductionData::~ReductionData ()
{
  INFO("(%d) destructor\n", index);
}

DoubleMsg * ReductionData::get ()
{
  return new DoubleMsg(x);
}

void ReductionData::pup (PUP::er &p)
{
  CBase_ReductionData::pup(p);
  p|N;
  p|index;
  p|x;

  if(p.isUnpacking())
  {
    INFO("(%d) pup unpacking\n", index);
  }
  else
  {
    if(p.isSizing()) { INFO("(%d) pup sizing\n", index); }
    else { INFO("(%d) pup packing\n", index); }
  }
}

#include "reductiondata.def.h"
