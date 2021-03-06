/** @file
 *
 * Some data types.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __TYPES_H
#define __TYPES_H

/** The matrix types. */
enum matrix_t
{
  /** A full matrix. */
  full,

  /** A matrix with decay. */
  decay,

  /** A diagonal matrix. */
  diagonal
};

/** Data type for Matrix PEMap. */
struct PEMap_Node_t
{
  /** The array index of the Node. */
  int index[2];

  /** The PE of the Node. */
  int PE;

  /** The norm of the Node. */
  double norm;
};

/** Data type for convolution PEMap. */
struct PEMap_MultiplyElement_t
{
  /** The array index of the convolution element. */
  int index[3];

  /** The PE. */
  int PE;

  /** The norm product. */
  double norm_product;
};

/** The operation to test. */
enum operation_t
{
  /** A multiplication. */
  multiply,

  /** An addition. */
  add,

  /** A trace. */
  trace,

  /** Add a scaled identity matrix. */
  addIdentity,

  /** Scale a matrix by a number. */
  scale,

  /** An SP2 iteration. */
  SP2
};

#endif
