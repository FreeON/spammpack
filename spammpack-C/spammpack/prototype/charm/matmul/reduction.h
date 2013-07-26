#ifndef __REDUCTION_H
#define __REDUCTION_H

#include "reduction.decl.h"

class ReductionData : public CBase_ReductionData
{
  private:

    double x;

  public:

    ReductionData ();
    ReductionData (CkMigrateMessage *msg);
    ~ReductionData ();
    DoubleMsg * get ();
    void pup (PUP::er &p);
};

class Reduction : public CBase_Reduction
{
  private:

    int N;
    int index;
    int callCount;
    double *result;
    CProxy_ReductionData a;

  public:

    Reduction (int N, CProxy_ReductionData a);
    Reduction (CkMigrateMessage *msg);
    ~Reduction ();
    void reduce (CkCallback &cb);
    void pup (PUP::er &p);
};

#endif
