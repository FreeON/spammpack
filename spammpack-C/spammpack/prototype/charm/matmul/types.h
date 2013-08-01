/** @file
 *
 * Some data types.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __TYPES_H
#define __TYPES_H

/** The initialization type. */
enum init_t
{
  /** Random matrix. */
  initRandom,
  
  /** Zero matrix. */
  initZero,

  /** Matrices with decay. */
  initDecay
};

/** The matrix types. */
enum matrix_t
{
  /** A full matrix. */
  full,

  /** A matrix with decay. */
  decay
};

#endif
