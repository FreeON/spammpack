#include "reductions.decl.h"

#include <time.h>

class Driver : public CBase_Driver
{
  private:

    struct timespec start_timer;
    struct timespec end_timer;

  public:

    Driver (CkArgMsg *msg)
    {
      int array_size = 60;
      CProxy_Worker w = CProxy_Worker::ckNew(array_size, array_size);
      CkCallback *cb = new CkCallback(CkReductionTarget(Driver, done), thisProxy);
      w.ckSetReductionClient(cb);
      clock_gettime(CLOCK_MONOTONIC_RAW, &start_timer);
      w.reduce();
    }

    void done (int result)
    {
      struct timespec timer;
      clock_gettime(CLOCK_MONOTONIC_RAW, &end_timer);
      CkPrintf("% 2d %f %d\n", CkNumPes(),
          end_timer.tv_sec+end_timer.tv_nsec/1.0e9
          -(start_timer.tv_sec+start_timer.tv_nsec/1.0e9),
          result);
      CkExit();
    }
};

class Worker : public CBase_Worker
{
  public:

    Worker () {}

    Worker (CkMigrateMessage *msg) {}

    void reduce (void)
    {
      int partial = 1;
      for(int i = 0; i < 200000; i++)
      {
        if(rand() < -0.5) { CkExit(); }
      }
      contribute(1*sizeof(int), &partial, CkReduction::sum_int);
    }
};

#include "reductions.def.h"
