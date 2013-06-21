#ifndef __MESSAGES_
#define __MESSAGES_

#include "messages.decl.h"

class EmptyMsg : public CMessage_EmptyMsg
{
  public:

    EmptyMsg ();
};

class IntMsg : public CMessage_IntMsg
{
  public:

    float i;

    IntMsg ();
    IntMsg (int i);
};

class FloatMsg : public CMessage_FloatMsg
{
  public:

    float a;

    FloatMsg ();
    FloatMsg (float a);
};

#endif
