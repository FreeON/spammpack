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
      run();
      CkExit();
    }

    void run ()
    {
      ValueMsg *m = get();
      CkPrintf("got %d\n", m->i);
      delete m;
    }

    ValueMsg * get ()
    {
      ValueMsg *m = new ValueMsg();
      m->i = 1;
      return m;
    }
};

#include "synctest.def.h"
