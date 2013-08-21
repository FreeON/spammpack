#include "migration.decl.h"

#include <math.h>
#include <stdlib.h>

#define NUMBER_ELEMENTS 100000

class Work : public CBase_Work
{
  private:

    unsigned int seed;
    double A[NUMBER_ELEMENTS];

  public:

    Work (void)
    {
      seed = 1;
      memset(A, 0, NUMBER_ELEMENTS*sizeof(double));
      usesAtSync = true;
    }

    Work (CkMigrateMessage *msg) {}

    void doSomething ()
    {
      /* Do some work. */
      for(int i = 0; i < NUMBER_ELEMENTS; i++)
      {
        A[i] += rand_r(&seed)/(double) RAND_MAX;
        if(A[i] < 0) { CkExit(); }
        A[i] = sqrt(A[i]);
        for(int i = 0; i < 10; i++)
        {
          if(sin(i) < -2) { CkExit(); }
        }
      }
      /* Done. */
      AtSync();
    }

    void pup (PUP::er &p)
    {
      CBase_Work::pup(p);
      p|seed;
      PUParray(p, A, NUMBER_ELEMENTS);
    }

    void ResumeFromSync (void)
    {
      contribute();
    }
};

class Main : public CBase_Main
{
  private:

    int iteration;
    int maxIteration;
    CProxy_Work work;

  public:

    Main (CkArgMsg *msg)
    {
      const int N = 1000;

      work = CProxy_Work::ckNew(N);
      work.ckSetReductionClient(new CkCallback(CkReductionTarget(Main, done), thisProxy));

      iteration = 0;
      maxIteration = 50;

      /* First iteration. */
      CkPrintf("iteration %d\n", iteration+1);
      work.doSomething();
    }

    void done (void)
    {
      CkPrintf("iteration %d done\n", iteration+1);
      iteration++;
      if(iteration < maxIteration)
      {
        CkPrintf("iteration %d\n", iteration+1);
        work.doSomething();
      }

      else
      {
        CkPrintf("done\n");
        CkExit();
      }
    }
};

#include "migration.def.h"
