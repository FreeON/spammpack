/** @file
 *
 * The header file for the SpAMM_Node class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __BLOCK_H
#define __BLOCK_H

#include <charm++.h>

class Block
{
  private:

    int blocksize;
    double *block;

  public:

    Block (void);
    void pup (PUP::er &p);
    void set (const int blocksize, const double *const A);
    void multiply (Block A, Block B);
};

#endif
