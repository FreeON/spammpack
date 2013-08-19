#include "migration.decl.h"

class DataMsg : public CMessage_DataMsg
{
  public:

    int numberElements;
    double *A;
};

class Data : public CBase_Data
{
  private:

    int numberElements;
    double *A;

  public:

    Data (void)
    {
      numberElements = (int) floor(rand()/(double) RAND_MAX * 3000 + 10000);
      A = new double[numberElements];
      for(int i = 0; i < numberElements; i++)
      {
        A[i] = rand()/(double) RAND_MAX;
      }
    }

    Data (CkMigrateMessage *msg) {}

    ~Data (void)
    {
      delete[] A;
    }

    DataMsg * info (void)
    {
      DataMsg *m = new (numberElements) DataMsg();
      memcpy(m->A, A, numberElements*sizeof(double));
      m->numberElements = numberElements;
      return m;
    }

    void pup (PUP::er &p)
    {
      CBase_Data::pup(p);
      p|numberElements;
      if(p.isUnpacking()) { A = new double[numberElements]; }
      /* Following the example in section 6.3, we pup the A array even if
       * numberElements == 0. */
      PUParray(p, A, numberElements);
    }
};

class Work : public CBase_Work
{
  private:

    int numberElements;
    double *A;
    CProxy_Data data;

  public:

    Work (CProxy_Data data)
    {
      numberElements = 0;
      A = NULL;
      this->data = data;
    }

    Work (CkMigrateMessage *msg) {}

    ~Work (void)
    {
      delete[] A;
    }

    void doSomething (CkCallback &cb)
    {
      DataMsg *msg = data(thisIndex).info();

      if(numberElements == 0)
      {
        numberElements = msg->numberElements;
        A = new double[numberElements];
      }

      else
      {
        /* Sanity check. */
        if(numberElements != msg->numberElements)
        {
          CkPrintf("mismatching numberElements\n");
          CkExit();
        }
      }

      for(int i = 0; i < numberElements; i++)
      {
        A[i] = msg->A[i]+rand()/(double) RAND_MAX;
        /* Do some extra work. */
        if(A[i] < 0) { CkExit(); }
        for(int j = 0; j < 100; j++)
        {
          if(sin(j) < -2) { CkExit(); }
        }
      }

      delete msg;

      contribute(cb);
    }

    void pup (PUP::er &p)
    {
      CBase_Work::pup(p);
      p|numberElements;
      if(p.isUnpacking()) { A = new double[numberElements]; }
      /* Following the example in section 6.3, we pup the A array even if
       * numberElements == 0. */
      PUParray(p, A, numberElements);
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
      thisProxy.run();
    }

    void run (void)
    {
      const int N = 1000;

      data = CProxy_Data::ckNew(N);
      work = CProxy_Work::ckNew(data, N);

      for(int iteration = 0; iteration < 10; iteration++)
      {
        CkPrintf("iteration %d\n", iteration+1);
        work.doSomething(CkCallbackResumeThread());
      }
      CkPrintf("done\n");
      CkExit();
    }
};

#include "migration.def.h"
