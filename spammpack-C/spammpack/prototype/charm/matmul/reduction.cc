#include "reduction.h"
#include "logger.h"
#include "messages.h"

ReductionData::ReductionData ()
{
  x = rand()/(double) RAND_MAX;
}

ReductionData::ReductionData (CkMigrateMessage *msg)
{
}

ReductionData::~ReductionData ()
{
}

DoubleMsg * ReductionData::get ()
{
  return new DoubleMsg(x);
}

void ReductionData::pup (PUP::er &p)
{
  CBase_ReductionData::pup(p);
  p|x;
}

Reduction::Reduction (int N, CProxy_ReductionData a)
{
  this->N = N;
  this->index = thisIndex.x*N*N + thisIndex.y*N + thisIndex.z;
  this->callCount = 0;
  this->result = NULL;
  this->a = a;
  INFO("(%d) constructor\n", index);
}

Reduction::Reduction (CkMigrateMessage *msg)
{
}

Reduction::~Reduction ()
{
  INFO("(%d) destructor\n", index);
  delete[] result;
}

void Reduction::reduce (CkCallback &cb)
{
  INFO("(%d) reduce\n", index);
  if(callCount > 0)
  {
    ABORT("(%d) this chare has been called before\n");
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
    INFO("(%d) pup unpacking\n", index);
  }
  else
  {
    if(p.isSizing()) { INFO("(%d) pup sizing\n", index); }
    else { INFO("(%d) pup packing\n", index); }
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
