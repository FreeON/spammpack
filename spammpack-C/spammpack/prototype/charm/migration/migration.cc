#include "migration.decl.h"

class IntMsg : public CMessage_IntMsg
{
  public:

    int i;

    IntMsg (int i)
    {
      this->i = i;
    }
};

class Data : public CBase_Data
{
  private:

    int i;

  public:

    Data (void)
    {
      i = thisIndex;
    }

    Data (CkMigrateMessage *msg) {}

    IntMsg * info (void)
    {
      return new IntMsg(i);
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
#define SUSPEND
#ifdef SUSPEND
      /* When this code block is used, the code segfaults during object
       * migration when the load balancer kicks in. My best guess here is that
       * since the call to Data::info() suspends this chare, it can sometimes
       * accidentally happen that the load balancer migrates the chare before
       * it continues.
       */
      IntMsg *msg = stuff(thisIndex).info();
      delete msg;
#endif

      if(numberElements == 0)
      {
        numberElements = (int) (floor(rand()/(double) RAND_MAX * 100 + 100000));
        data = new double[numberElements];
      }

      for(int i = 0; i < numberElements; i++)
      {
        data[i] = rand()/(double) RAND_MAX;
        if(data[i] < 0) { CkExit(); }
        data[i] = sqrt(data[i]);
      }

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
