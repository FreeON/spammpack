#include "Futures.decl.h"

class IntMsg : public CMessage_IntMsg
{
  public:

    unsigned int counter;

    IntMsg ()
    {
      counter = 0;
    }
};

class Node : public CBase_Node
{
  private:

    CProxy_Node *child[8];

  public:

    Node ()
    {
      for(int i = 0; i < 8; i++)
      {
        child[i] = NULL;
      }
    }

    void compute (int tier, int depth, CkFuture f)
    {
      IntMsg *result = new IntMsg();
      if(tier == depth)
      {
        result->counter = 1;
      }

      else
      {
        CkFuture f_child[8];
        for(int i = 0; i < 8; i++)
        {
          if(child[i] == NULL)
          {
            child[i] = new CProxy_Node;
            *(child[i]) = CProxy_Node::ckNew();
          }
          f_child[i] = CkCreateFuture();
          child[i]->compute(tier+1, depth, f_child[i]);
        }

        for(int i = 0; i < 8; i++)
        {
          IntMsg *m = (IntMsg*) CkWaitFuture(f_child[i]);
          result->counter += m->counter;
          delete m;
          CkReleaseFuture(f_child[i]);
        }
      }

      CkSendToFuture(f, result);
    }
};

class Futures : public CBase_Futures
{
  public:

    Futures (CkArgMsg *msg)
    {
      int depth = 1;
      if(msg->argc > 1)
      {
        depth = strtol(msg->argv[1], NULL, 10);
      }
      CkPrintf("running on %d PEs\n", CkNumPes());
      thisProxy.run(depth);
    }

    void run (int depth)
    {
      CkPrintf("depth = %d\n", depth);
      CProxy_Node root = CProxy_Node::ckNew();

      CkFuture f = CkCreateFuture();
      root.compute(0, depth, f);
      IntMsg *m = (IntMsg*) CkWaitFuture(f);
      CkPrintf("done, counter = %d (%d)\n", m->counter, 1 << (3*depth));
      delete m;
      CkExit();
    }
};

#include "Futures.def.h"
