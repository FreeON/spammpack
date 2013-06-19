#include "fib.decl.h"

class ValueMsg
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

      run(strtol(msg->argv[1], NULL, 10));
    }

    void run (int n)
    {
      printf("calculating F(%d)... ", n); fflush(stdout);
      CkFuture f = CkCreateFuture();
      printf("created future... "); fflush(stdout);
      CProxy_Fib::ckNew(n, f);
      printf("created instance of Fib proxy... "); fflush(stdout);
      ValueMsg *m = (ValueMsg*) CkWaitFuture(f);
      printf("F(%d) = %d\n", n, m->value); fflush(stdout);
    }
};

class Fib : public CBase_Fib
{
  public:

    Fib (int n, CkFuture f)
    {
      printf("Fib: constructor with n = %d\n", n); fflush(stdout);
      run(n, f);
    }

    void run (int n, CkFuture f)
    {
      ValueMsg *m = new ValueMsg();

      printf("run: calculating F(%d)\n", n); fflush(stdout);

      if(n < 0)
      {
        CkPrintf("n can not be negative\n");
        CkExit();
      }

      else if(n <= 1)
      {
        printf("Fib: direct evaluation\n"); fflush(stdout);
        m->value = n;
      }

      else
      {
        CkFuture f1 = CkCreateFuture();
        CkFuture f2 = CkCreateFuture();
        CProxy_Fib::ckNew(n-1, f1);
        CProxy_Fib::ckNew(n-2, f2);
        ValueMsg *m1 = (ValueMsg*) CkWaitFuture(f1);
        ValueMsg *m2 = (ValueMsg*) CkWaitFuture(f2);
        m->value = m1->value+m2->value;
        delete m1;
        delete m2;
      }

      CkSendToFuture(f, m);
    }
};

#include "fib.def.h"
