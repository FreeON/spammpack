#include "messages.h"

EmptyMsg::EmptyMsg () {}

IntMsg::IntMsg ()
{
  this->i = 0;
}

IntMsg::IntMsg (int i)
{
  this->i = i;
}

FloatMsg::FloatMsg ()
{
  this->a = 0;
}

FloatMsg::FloatMsg (float a)
{
  this->a = a;
}

#include "messages.def.h"
