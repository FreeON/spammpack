#include "synctest.decl.h"

class ValueMsg : public CMessage_ValueMsg
{
  public:
    int i;
};

class Main : public CBase_Main
{
  public:

    Main (CkArgMsg *msg)
    {
      CkPrintf("starting...\n");
      thisProxy.run();
    }

    void run()
    {
      CProxy_Worker r = CProxy_Worker::ckNew();
      CkPrintf("calling get method...\n");
      ValueMsg *m = r.get();
      CkPrintf("received %d from chare\n", m->i);
      delete m;
    }
};

class Worker : public CBase_Worker
{
  public:

    Worker ()
    {
      CkPrintf("initializing chare...\n");
    }

    ValueMsg * get ()
    {
      ValueMsg *m = new ValueMsg();
      m->i = 1;
      return m;
    }
};

#include "synctest.def.h"
