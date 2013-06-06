#include "linkedList.decl.h"

class LinkedListNode : public CBase_LinkedListNode
{
  private:

    int ID;
    CProxy_LinkedListNode *next;

  public:

    LinkedListNode ()
    {
      ID = -1;
      next = NULL;
      CkPrintf("initializing node\n");
    }

    void set (const int ID)
    {
      CkPrintf("setting ID to %i\n", ID);
      this->ID = ID;
      CkPrintf("done setting ID to %i\n", ID);
    }

    void append (const int ID)
    {
      if(next != NULL) { next->append(ID); }
      else
      {
        *next = CProxy_LinkedListNode::ckNew();
        next->set(ID);
      }
      CkPrintf("done appending\n");
    }
};

class Main : public CBase_Main
{
  private:

    CProxy_LinkedListNode first;

  public:

    Main (CkArgMsg *msg)
    {
      first = CProxy_LinkedListNode::ckNew();
      first.set(0);
      for(int i = 1; i < 10; i++)
      {
        first.append(i);
      }
    }
};

#include "linkedList.def.h"
