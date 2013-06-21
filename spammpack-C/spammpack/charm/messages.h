#ifndef __MESSAGES_
#define __MESSAGES_

#include "messages.decl.h"

class EmptyMsg : public CMessage_EmptyMsg
{
  public:

    EmptyMsg ();
};

class GetMsg : public CMessage_GetMsg
{
  public:

    float a;
};

#endif
