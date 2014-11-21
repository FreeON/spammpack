#include "openmp.decl.h"

#include <omp.h>

class Work : public CBase_Work
{
  public:

    Work (void) {}
    Work (CkMigrateMessage *msg) {}

    void doSomething (CkCallback &cb)
    {
#pragma omp parallel for
      for(int i = 0; i < 4; i++)
      {
        CkPrintf("array index %d, PE %d, thread %d, i = %d\n", thisIndex,
            CkMyPe(), omp_get_thread_num(), i);
      }
      contribute(cb);
    }
};

class Main : public CBase_Main
{
  public:

    Main (CkArgMsg *arg)
    {
#pragma omp parallel
      {
#pragma omp master
        {
          CkPrintf("there are %d threads\n", omp_get_num_threads());
        }
      }
      thisProxy.run();
    }

    void run (void)
    {
      CProxy_Work work = CProxy_Work::ckNew(4);
      work.doSomething(CkCallbackResumeThread());
      CkWaitQD();
      CkExit();
    }
};

#include "openmp.def.h"
