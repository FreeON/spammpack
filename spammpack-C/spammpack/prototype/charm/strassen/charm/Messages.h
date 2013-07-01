#ifndef __MESSAGES_H
#define __MESSAGES_H

#include "Messages.decl.h"

class DoubleMsg : public CMessage_DoubleMsg
{
  public :

    double x;
    DoubleMsg ();
    DoubleMsg (double x);
};

#endif
