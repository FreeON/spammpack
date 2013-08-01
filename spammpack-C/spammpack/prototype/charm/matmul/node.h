/** @file
 *
 * The header file for the Node class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __NODE_H
#define __NODE_H

#include "types.h"

#include "node.decl.h"

#include <string>

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
    unsigned int index;

    /** The norm of this matrix block. */
    double norm;

    /** The square of the norm of this matrix block. */
    double norm_2;

    /** A flag to indicate whether tierNode is set. */
    bool tierNodeSet;

    /** The nodes for the next tier. */
    CProxy_Node tierNode;

    /** The submatrix of size blocksize x blocksize. */
    double *block;

  public:

    Node ();
    Node (int N, int depth, int blocksize, int tier);
    ~Node ();
    Node (CkMigrateMessage *msg);
    virtual void pup (PUP::er &p);
    NodeBlockMsg * getBlock ();
    DoubleMsg * get (int i, int j);
    NodeInfoMsg * info ();
    void printPE (CkCallback &cb);
    void random (CkCallback &cb);
    void zero (CkCallback &cb);
    void set (int numberElements, double *A, CkCallback &cb);
    void initialize (int initType, CkCallback &cb);
    void printLeafPes (CkCallback &cb);
    void add (int blocksize, double *A);
    void updateNorms (CkCallback &cb);
    void setTierNode (CProxy_Node tierNode, CkCallback &cb);
    void print (std::string tag);
};

#endif
