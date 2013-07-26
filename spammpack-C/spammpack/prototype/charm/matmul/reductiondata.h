/** @file
 *
 * The header file for the ReductionData class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __REDUCTIONDATA_H
#define __REDUCTIONDATA_H

#include "reductiondata.decl.h"

class ReductionData : public CBase_ReductionData
{
  private:

    int N;
    int index;
    double x;

  public:

    ReductionData (int N);
    ReductionData (CkMigrateMessage *msg);
    ~ReductionData ();
    DoubleMsg * get ();
    void pup (PUP::er &p);
};

#endif
