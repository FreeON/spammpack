#include "migration.decl.h"

#include <math.h>
#include <stdlib.h>

#ifdef WITH_DATA
#warning Compiling with DataMsg.
#else
#warning Compiling _without_ DataMsg.
#endif

#ifdef WITH_DATA
class DataMsg : public CMessage_DataMsg
{
  public:

    int i;

    DataMsg (int i)
    {
      this->i = i;
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
      return new DataMsg(numberElements);
    }

    void pup (PUP::er &p)
    {
      CBase_Data::pup(p);
      p|numberElements;
    }
};
#endif

class Work : public CBase_Work
{
  private:

    int numberElements;
    double *A;
#ifdef WITH_DATA
    CProxy_Data data;
#endif

  public:

#ifdef WITH_DATA
    Work (CProxy_Data data)
    {
      this->data = data;
#else
    Work (void)
    {
#endif
      numberElements = 0;
      A = NULL;
    }

    Work (CkMigrateMessage *msg) {}

    ~Work (void)
    {
      delete[] A;
    }

    void doSomething (CkCallback &cb)
    {
#ifdef WITH_DATA
      DataMsg *msg = data(thisIndex).info();

      if(numberElements == 0)
      {
        numberElements = msg->i;
        A = new double[numberElements];
        memset(A, 0, numberElements*sizeof(double));
      }

      else if(numberElements != msg->i)
      {
        CkPrintf("numberElement mismatch\n");
        CkExit();
      }

      for(int i = 0; i < numberElements; i++)
      {
        A[i] += rand()/(double) RAND_MAX;
        if(A[i] < 0) { CkExit(); }
        A[i] = sqrt(A[i]);
        for(int i = 0; i < 10; i++)
        {
          if(sin(i+msg->i) < -2) { CkExit(); }
        }
      }

      delete msg;
#else
#warning Compiling _without_ DataMsg.
      if(numberElements == 0)
      {
        numberElements = 100000;
        A = new double[numberElements];
        memset(A, 0, numberElements*sizeof(double));
      }

      for(int i = 0; i < numberElements; i++)
      {
        A[i] += rand()/(double) RAND_MAX;
        if(A[i] < 0) { CkExit(); }
        A[i] = sqrt(A[i]);
        for(int i = 0; i < 10; i++)
        {
          if(sin(i) < -2) { CkExit(); }
        }
      }
#endif
      contribute(cb);
    }

    void pup (PUP::er &p)
    {
      CBase_Work::pup(p);
      p|numberElements;
      if(p.isUnpacking()) { A = new double[numberElements]; }
      PUParray(p, A, numberElements);
#ifdef WITH_DATA
      p|data;
#endif
    }
};

class Main : public CBase_Main
{
  private:

#ifdef WITH_DATA
    CProxy_Data data;
#endif
    CProxy_Work work;

  public:

    Main (CkArgMsg *msg)
    {
      thisProxy.run();
    }

    void run (void)
    {
      const int N = 1000;

#ifdef WITH_DATA
      data = CProxy_Data::ckNew(N);
      work = CProxy_Work::ckNew(data, N);
#else
      work = CProxy_Work::ckNew(N);
#endif

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
