#include "spamm.decl.h"

class SpAMMMatrix : public CBase_SpAMMMatrix
{
  public:

    SpAMMMatrix (const int N)
    {
    }

    SpAMMMatrix (CkMigrateMessage *msg) {}
};

#include "spamm.def.h"
