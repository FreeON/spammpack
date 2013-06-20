#include "fib.decl.h"

class ValueMsg
{
  public:

    int value;
};

class Main : public CBase_Main
{
  public:

    Main (CkMigrateMessage *msg)
    {
      CkPrintf("migrating...\n");
    }

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

    Fib (CkMigrateMessage *msg)
    {
      CkPrintf("[Fib] migrating...\n");
    }

    Fib (bool amIRoot, int n, CkFuture f)
    {
      printf("Fib%s: constructor with n = %d\n", (amIRoot ? " (root)" : ""), n); fflush(stdout);
      run(amIRoot, n, f);
    }

    void run (bool amIRoot, int n, CkFuture f)
    {
      ValueMsg *m = new ValueMsg();

      /* Set queueing strategy. */
      CkSetQueueing(m, CK_QUEUEING_FIFO);

      printf("[Fib::run] calculating F(%d)... ", n); fflush(stdout);

      if(n < 0)
      {
        CkPrintf("n can not be negative\n");
        CkExit();
      }

      else if(n <= 1)
      {
        printf("direct evaluation... "); fflush(stdout);
        m->value = n;
      }

      else
      {
        printf("recursive evaluation... "); fflush(stdout);
        CkFuture f1 = CkCreateFuture();
        CkFuture f2 = CkCreateFuture();
        printf("created futures... "); fflush(stdout);
        CProxy_Fib::ckNew(false, n-1, f1);
        CProxy_Fib::ckNew(false, n-2, f2);
        printf("waiting for results...\n"); fflush(stdout);
        ValueMsg *m1 = (ValueMsg*) CkWaitFuture(f1);
        printf("received m1: %d... ", m1->value); fflush(stdout);
        ValueMsg *m2 = (ValueMsg*) CkWaitFuture(f2);
        printf("received m2: %d... ", m2->value); fflush(stdout);
        m->value = m1->value+m2->value;
        printf("done adding\n"); fflush(stdout);
        //delete m1;
        //delete m2;
        //printf("done deleting\n"); fflush(stdout);
      }

      if(amIRoot)
      {
        printf("\n"); fflush(stdout);
        CkPrintf("F[%d] = %d\n", n, m->value);
        delete m;
        CkExit();
      }

      else
      {
        printf("sending to future\n"); fflush(stdout);
        CkSendToFuture(f, m);
      }
    }
};

#include "fib.def.h"
