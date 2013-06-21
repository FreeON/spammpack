#include "messages.h"

EmptyMsg::EmptyMsg () {}

GetMsg::GetMsg ()
{
  this->a = 0;
}

GetMsg::GetMsg (float a)
{
  this->a = a;
}

#include "messages.def.h"
