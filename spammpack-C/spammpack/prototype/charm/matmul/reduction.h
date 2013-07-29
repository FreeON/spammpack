/** @file
 *
 * The header file for the Reduction class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __REDUCTION_H
#define __REDUCTION_H

#include "reductiondata.h"

#include "reduction.decl.h"

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
