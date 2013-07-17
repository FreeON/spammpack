#include "messages.h"

DoubleMsg::DoubleMsg (double x)
{
  this->x = x;
}

IntMsg::IntMsg (int i)
{
  this->i = i;
}

MatrixInfoMsg::MatrixInfoMsg ()
{
  rootNull = true;
}

MatrixInfoMsg::MatrixInfoMsg (CProxy_Node root)
{
  this->rootNull = false;
  this->root = root;
}

#include "messages.def.h"
