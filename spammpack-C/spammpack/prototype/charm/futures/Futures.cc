#include "Futures.decl.h"

#include <bitset>

unsigned int counter;

class EmptyMsg : public CMessage_EmptyMsg
{
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

    void compute (int tier, unsigned int index, int depth, CkFuture f)
    {
      if(tier == depth)
      {
        /* Do some work. */
        for(int i = 0; i < 1000; i++)
        {
          double x = rand()/(double) RAND_MAX;
          if(x == -1)
          {
            CkExit();
          }
        }
        counter++;
#ifdef DEBUG
        std::bitset<32> bits = index;
        CkPrintf("returning from %s\n", bits.to_string().c_str());
#endif
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
          child[i]->compute(tier+1, (index << 3) | i, depth, f_child[i]);
        }

        for(int i = 0; i < 8; i++)
        {
          EmptyMsg *m = (EmptyMsg*) CkWaitFuture(f_child[i]);
          delete m;
          CkReleaseFuture(f_child[i]);
        }
      }

      CkSendToFuture(f, new EmptyMsg());
    }
};

class Futures : public CBase_Futures
{
  public:

    Futures (CkArgMsg *msg)
    {
      int depth = 1;
      counter = 0;
      if(msg->argc > 1)
      {
        depth = strtol(msg->argv[1], NULL, 10);
      }
      thisProxy.run(depth);
    }

    void run (int depth)
    {
      CkPrintf("Testing depth = %d\n", depth);
      CProxy_Node root = CProxy_Node::ckNew();

      CkFuture f = CkCreateFuture();
      root.compute(0, 1, depth, f);
      EmptyMsg *m = (EmptyMsg*) CkWaitFuture(f);
      delete m;
      CkPrintf("done, counted %d\n", counter);
      CkExit();
    }
};

#include "Futures.def.h"
