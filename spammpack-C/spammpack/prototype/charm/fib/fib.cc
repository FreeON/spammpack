#include "fib.decl.h"

class ValueMsg : public CMessage_ValueMsg
{
  public:

    int value;
};

class Main : public CBase_Main
{
  public:

    Main (CkArgMsg *msg)
    {
      if(msg->argc != 2)
      {
        CkPrintf("missing argument\n");
        CkExit();
      }

      CProxy_Fib::ckNew(true, strtol(msg->argv[1], NULL, 10), CkFuture());
    }
};

class Fib : public CBase_Fib
{
  public:

    Fib (bool amIRoot, int n, CkFuture f)
    {
      run(amIRoot, n, f);
    }

    void run (bool amIRoot, int n, CkFuture f)
    {
      ValueMsg *m = new ValueMsg();

      if(n < 0)
      {
        CkPrintf("n can not be negative\n");
        CkExit();
      }

      else if(n <= 1)
      {
        /* Direct evaluation. */
        m->value = n;
      }

      else
      {
        /* Recursive evaluation. */
        CkFuture f1 = CkCreateFuture();
        CkFuture f2 = CkCreateFuture();
        CProxy_Fib::ckNew(false, n-1, f1);
        CProxy_Fib::ckNew(false, n-2, f2);
        ValueMsg *m1 = (ValueMsg*) CkWaitFuture(f1);
        ValueMsg *m2 = (ValueMsg*) CkWaitFuture(f2);
        m->value = m1->value+m2->value;
        delete m1;
        delete m2;
      }

      if(amIRoot)
      {
        /* Done. */
        CkPrintf("F[%d] = %d\n", n, m->value);
        delete m;
        CkExit();
      }

      else
      {
        CkSendToFuture(f, m);
      }
    }
};

#include "fib.def.h"
