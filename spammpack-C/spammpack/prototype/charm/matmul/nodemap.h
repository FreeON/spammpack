/** @file
 *
 * The header file for the NodeMap class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */


#ifndef __NODEMAP_H
#define __NODEMAP_H

#include "nodemap.decl.h"

class NodeMap : public CkArrayMap
{
  private:

    int initialPE;

  public:

    NodeMap (int intialPE);
    int procNum (int, const CkArrayIndex &element);
    void pupulateInitial (int, int numInitial, void *msg, CkArrMgr *mgr);
};

#endif
