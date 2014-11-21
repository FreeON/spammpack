#include "charmtest_01.decl.h"

int numElements;
CProxy_Main mainProxy;

class Main : public CBase_Main
{
  public:

    Main (CkArgMsg *msg)
    {
      numElements = 5;

      CkPrintf("Hello World starting with %d chares "
          "on %d PEs\n", numElements, CkNumPes());

      mainProxy = thisProxy;
      CProxy_Hello helloArray = CProxy_Hello::ckNew(numElements);
      helloArray[0].sayHi(-1);
    }

    void done ()
    {
      CkExit();
    }

    Main (CkMigrateMessage *msg) {}
};

class Hello : public CBase_Hello
{
  public:

    Hello ()
    {
      CkPrintf("(%d) initializing, numElements = %d\n", thisIndex,
          numElements);
    }

    void sayHi (int from)
    {
      CkPrintf("(%d) Hello on processor %d (called by % d)\n",
          thisIndex, CkMyPe(), from);

      if(thisIndex < (numElements-1))
      {
        thisProxy[thisIndex+1].sayHi(thisIndex);
      }

      else
      {
        CkPrintf("(%d) I am calling main to let it know we are done\n",
            thisIndex);
        mainProxy.done();
      }
    }

    Hello (CkMigrateMessage *msg) {}
};

#include "charmtest_01.def.h"
