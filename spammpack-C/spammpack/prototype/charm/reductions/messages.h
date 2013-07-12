#ifndef __MESSAGES_H
#define __MESSAGES_H

#include "messages.decl.h"

class DoubleMsg : public CMessage_DoubleMsg
{
  public:

    double x;
    DoubleMsg ();
};

#endif
