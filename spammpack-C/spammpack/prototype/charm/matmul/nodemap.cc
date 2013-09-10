/** @file
 *
 * The implementation of the NodeMap class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "nodemap.h"
#include "logger.h"

NodeMap::NodeMap (int intialPE)
{
  INFO("I am %p\n", this);
  this->initialPE = intialPE;
}

int NodeMap::procNum (int, const CkArrayIndex &element)
{
  INFO("returning current PE\n");
  return CkMyPe();
}

void NodeMap::pupulateInitial (int, int numInitial, void *msg, CkArrMgr *mgr)
{
  INFO("inserting %d initial elements\n", numInitial);

  if(numInitial == 0) { return; }

  for(int i = 0; i < numInitial; i++)
  {
  }
}

#include "nodemap.def.h"
