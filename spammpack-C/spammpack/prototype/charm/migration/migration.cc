#include "migration.decl.h"

#include <math.h>
#include <stdlib.h>

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

    unsigned int seed;
    double A[NUMBER_ELEMENTS];
    CProxy_Main mainProxy;
    CProxy_Data data;

  public:

    Work (CProxy_Main mainProxy, CProxy_Data data)
    {
      seed = 1;
      memset(A, 0, NUMBER_ELEMENTS*sizeof(double));
      this->mainProxy = mainProxy;
      this->data = data;
      usesAtSync = true;
    }

    Work (CkMigrateMessage *msg) {}

    void doSomething (CkCallback &cb)
    {
      /* Get some data. */
      DataMsg *msg = data(thisIndex).info();

      /* Check PE. */
      if(CkMyPe() != msg->PE)
      {
        mainProxy.addMismatched();
      }

      /* Do some work. */
      for(int i = 0; i < NUMBER_ELEMENTS; i++)
      {
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
      p|seed;
      PUParray(p, A, NUMBER_ELEMENTS);
      p|mainProxy;
      p|data;
    }
};

class Main : public CBase_Main
{
  private:

    int numberMismatched;
    CProxy_Data data;
    CProxy_Work work;

  public:

    Main (CkArgMsg *msg)
    {
      const int N = 1000;

      data = CProxy_Data::ckNew(N);
      work = CProxy_Work::ckNew(thisProxy, data, N);
      thisProxy.iterate();
    }

    void iterate (void)
    {
      LBDatabase *db = LBDatabaseObj();

      for(int iteration = 0; iteration < 5; iteration++)
      {
        CkPrintf("iteration %d, ", iteration+1);
        numberMismatched = 0;
        work.doSomething(CkCallbackResumeThread());
        CkPrintf("%d mismatched elements\n", numberMismatched);
        db->StartLB();
        CkWaitQD();
      }

      CkPrintf("done\n");
      CkExit();
    }

    void addMismatched (void)
    {
      numberMismatched++;
    }
};

#include "migration.def.h"
