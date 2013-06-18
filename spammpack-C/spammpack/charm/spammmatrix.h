#ifndef __SPAMMMATRIX_H
#define __SPAMMMATRIX_H

#include "spammmatrix.decl.h"

class SpAMMMatrix : public CBase_SpAMMMatrix
{
  private:

    int N;

  public:

    SpAMMMatrix (const int N);
    SpAMMMatrix (CkMigrateMessage *msg);
};

#endif
