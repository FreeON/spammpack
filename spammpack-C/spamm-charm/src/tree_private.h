/** @file
 *
 * The header file for the data structure of trees.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef TREE_PRIVATE_H
#define TREE_PRIVATE_H

/** A tree node. */
struct tree_t
{
  /** The norm^2 of this node. */
  double norm_2;

  /** The tier. */
  int tier;

  /** Depending on where this tree_t object is, it can have pointers to
   * children nodes, or store actual data. */
  union
  {
    /** Links to children nodes. */
    struct tree_t *child[2][2];

    /** The matrix chunk. */
    void *chunk;
  } data;
};

#endif
