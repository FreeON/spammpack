#include "migration.decl.h"

#include <math.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NUMBER_ELEMENTS 100000

void setManualLB (void)
{
  CkPrintf("(%d) turning manual LB on\n", CkMyPe());
  TurnManualLBOn();
}

class DataMsg : public CMessage_DataMsg
{
  public:

    int i;
    int PE;

    DataMsg (int i, int PE)
    {
      this->i = i;
      this->PE = PE;
    }
};

class Data : public CBase_Data
{
  private:

    int numberElements;

  public:

    Data (void)
    {
      numberElements = (int) (floor(rand()/(double) RAND_MAX * 100 + 100000));
    }

    Data (CkMigrateMessage *msg) {}

    DataMsg * info (void)
    {
      return new DataMsg(numberElements, CkMyPe());
    }

    void pup (PUP::er &p)
    {
      CBase_Data::pup(p);
      p|numberElements;
    }
};

class Work : public CBase_Work
{
  private:

    int mismatchedPE;
    unsigned int seed;
    double A[NUMBER_ELEMENTS];
    CProxy_Data data;

  public:

    Work (CProxy_Data data)
    {
      seed = 1;
      memset(A, 0, NUMBER_ELEMENTS*sizeof(double));
      this->data = data;
      usesAtSync = true;
    }

    Work (CkMigrateMessage *msg) {}

    void doSomething (CkCallback &cb)
    {
      /* Get some data. */
      DataMsg *msg = data(thisIndex).info();

      /* Check PE. */
      mismatchedPE = (CkMyPe() == msg->PE ? 0 : 1);

      /* Do some work. */
#pragma omp parallel for
      for(int i = 0; i < NUMBER_ELEMENTS; i++)
      {
        CkPrintf("(%d) i = %d, running on thread %d\n", thisIndex, i, omp_get_thread_num());
        A[i] += rand_r(&seed)/(double) RAND_MAX;
        if(A[i] < 0) { CkExit(); }
        A[i] = sqrt(A[i]);
        for(int i = 0; i < 10; i++)
        {
          if(sin(i+msg->i) < -2) { CkExit(); }
        }
      }

      /* Cleanup. */
      delete msg;

      /* Done. */
      contribute(cb);
    }

    void pup (PUP::er &p)
    {
      CBase_Work::pup(p);
      p|mismatchedPE;
      p|seed;
      PUParray(p, A, NUMBER_ELEMENTS);
      p|data;
    }
};

class Main : public CBase_Main
{
  private:

    CProxy_Data data;
    CProxy_Work work;

  public:

    Main (CkArgMsg *msg)
    {
      const int N = 1000;

#ifdef _OPENMP
      CkPrintf("OpenMP hybrid on %d OpenMP nodes\n", omp_get_num_threads());
#endif

      data = CProxy_Data::ckNew(N);
      work = CProxy_Work::ckNew(data, N);
      thisProxy.iterate();
    }

    void iterate (void)
    {
      LBDatabase *db = LBDatabaseObj();

      for(int iteration = 0; iteration < 5; iteration++)
      {
        CkPrintf("iteration %d\n", iteration+1);
        work.doSomething(CkCallbackResumeThread());
        //db->StartLB();
        //CkWaitQD();
      }

      CkPrintf("done\n");
      CkExit();
    }
};

#include "migration.def.h"
