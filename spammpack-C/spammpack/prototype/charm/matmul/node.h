/** @file
 *
 * The header file for the Node class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __NODE_H
#define __NODE_H

#include "node.decl.h"
#include "types.h"
#include <map>
#include <list>

/** The Node class. */
class Node : public CBase_Node
{
  private:

    /** The matrix size. */
    int N;

    /** The submatrix size at the lowest tier. */
    int blocksize;

    /** The tree depth of the matrix. */
    int depth;

    /** The tier of this Node. */
    int tier;

    /** The lower row index. */
    int iLower;

    /** The upper row index. */
    int iUpper;

    /** The lower column index. */
    int jLower;

    /** The upper column index. */
    int jUpper;

    /** The linear index of this node. */
    int index;

    std::map<int, bool> callbackSet;
    std::map<int, CkCallback> cb;

    std::map<int, std::list<int> > childWorking;

    double *block;

  public:

    Node (int N, int depth, int blocksize, int tier, int iLower, int iUpper,
        int jLower, int jUpper);
    Node (CkMigrateMessage *msg);
    NodeBlockMsg * getBlock ();
    DoubleMsg * get (int i, int j);
    NodeInfoMsg * info ();
    EmptyMsg * printPE ();
    void initialize (int initType, CkCallback &cb);
    void printLeafPes (int index, CkCallback &cb);
};

#endif
