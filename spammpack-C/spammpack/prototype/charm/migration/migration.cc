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
      numberElements = (int) floor(rand()/(double) RAND_MAX * 3000 + 100000);
      A = new double[numberElements];
      for(int i = 0; i < numberElements; i++)
      {
        A[i] = rand()/(double) RAND_MAX;
      }
    }

    Data (CkMigrateMessage *msg) {}

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
      if(numberElements > 0)
      {
        if(p.isUnpacking())
        {
          A = new double[numberElements];
        }
        PUParray(p, A, numberElements);
      }
      else
      {
        if(p.isUnpacking()) { A = NULL; }
      }
    }
};

class Work : public CBase_Work
{
  private:

    int numberElements;
    double *data;
    CProxy_Data stuff;

  public:

    Work (CProxy_Data stuff)
    {
      numberElements = 0;
      data = NULL;
      this->stuff = stuff;
    }

    Work (CkMigrateMessage *msg) {}

    ~Work (void)
    {
      delete[] data;
    }

    void doSomething (CkCallback &cb)
    {
      DataMsg *msg = stuff(thisIndex).info();

      if(numberElements == 0)
      {
        numberElements = msg->numberElements;
        data = new double[numberElements];
      }

      for(int i = 0; i < numberElements; i++)
      {
        data[i] = msg->A[i]+rand()/(double) RAND_MAX;
      }

      delete msg;

      contribute(cb);
    }

    void pup (PUP::er &p)
    {
      CBase_Work::pup(p);
      p|numberElements;
      if(numberElements > 0)
      {
        if(p.isUnpacking())
        {
          data = new double[numberElements];
        }
        PUParray(p, data, numberElements);
      }

      else
      {
        if(p.isUnpacking()) { data = NULL; }
      }
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
      work.doSomething(CkCallbackResumeThread());
      CkExit();
    }
};

#include "migration.def.h"
