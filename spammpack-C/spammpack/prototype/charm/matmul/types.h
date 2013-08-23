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

#endif
