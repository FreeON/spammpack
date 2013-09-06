/** @file
 *
 * The header file for the messages.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __MESSAGES_H
#define __MESSAGES_H

#include "messages.decl.h"
#include "node.h"

/** A message containing some matrix information. */
class MatrixInfoMsg : public CMessage_MatrixInfoMsg
{
  public:

    /** The matrix size. */
    int N;

    /** The submatrix size at the lowest tier. */
    int blocksize;

    /** The tree depth of the matrix. */
    int depth;

    /** The padded size of the matrix. */
    int NPadded;

    MatrixInfoMsg (int N, int blocksize, int depth, int NPadded);
    bool equal (MatrixInfoMsg *b);
};

/** A message containing the Node chare array on a particular tier. */
class MatrixNodeMsg : public CMessage_MatrixNodeMsg
{
  public:

    /** The Node chare array. */
    CProxy_Node nodes;

    MatrixNodeMsg (CProxy_Node nodes);
};

/** A message containing information on matrix nodes. */
class NodeInfoMsg : public CMessage_NodeInfoMsg
{
  public:

    /** The linear index of this node. */
    int index;

    /** The norm of this matrix block. */
    double norm;

    /** The square of the norm of this matrix block. */
    double norm_2;

    NodeInfoMsg (int index, double norm, double norm_2);
};

/** A message containing a dense matrix. */
class DenseMatrixMsg : public CMessage_DenseMatrixMsg
{
  public:

    /** The dense matrix. */
    double *A;
};

/** A message for PEMaps. */
class PEMapMsg : public CMessage_PEMapMsg
{
  public:

    /** The PEMap, updated by calling updatePEMap(). */
    int *PEMap;

    /** The norms of the PEMap. */
    double *PEMap_norm;
};

#endif
