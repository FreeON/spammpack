#include "migration.decl.h"

#include <math.h>
#include <stdlib.h>

#define NUMBER_ELEMENTS 100000

class Work : public CBase_Work
{
  private:

    double A[NUMBER_ELEMENTS];

  public:

    Work (void)
    {
      memset(A, 0, NUMBER_ELEMENTS*sizeof(double));
    }

    Work (CkMigrateMessage *msg) {}

    void doSomething (CkCallback &cb)
    {
      for(int i = 0; i < NUMBER_ELEMENTS; i++)
      {
        A[i] += rand()/(double) RAND_MAX;
        if(A[i] < 0) { CkExit(); }
        A[i] = sqrt(A[i]);
        for(int i = 0; i < 10; i++)
        {
          if(sin(i) < -2) { CkExit(); }
        }
      }
      contribute(cb);
    }

    void pup (PUP::er &p)
    {
      CBase_Work::pup(p);
      PUParray(p, A, NUMBER_ELEMENTS);
    }
};

class Main : public CBase_Main
{
  private:

    CProxy_Work work;

  public:

    Main (CkArgMsg *msg)
    {
      thisProxy.run();
    }

    void run (void)
    {
      const int N = 1000;

      work = CProxy_Work::ckNew(N);

      for(int iteration = 0; iteration < 50; iteration++)
      {
        CkPrintf("iteration %d\n", iteration+1);
        work.doSomething(CkCallbackResumeThread());
      }
      CkPrintf("done\n");
      CkExit();
    }
};

#include "migration.def.h"
