#ifndef __MESSAGES_H
#define __MESSAGES_H

#include "messages.decl.h"

class DoubleMsg : public CMessage_DoubleMsg
{
  public:

    double x;
    DoubleMsg (double x);
};

class EmptyMsg : public CMessage_EmptyMsg
{
};

class IntMsg : public CMessage_IntMsg
{
  public:

    int i;
    IntMsg (int i);
};

#endif
