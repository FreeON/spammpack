/** @file
 *
 * The header file for the SpAMM_Node class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __SPAMM_NODE_H
#define __SPAMM_NODE_H

#include <charm++.h>

class SpAMM_Node
{
  public:

    void pup (PUP::er &p);
    void set (const int blocksize, const double *const A);
    void multiply (SpAMM_Node A, SpAMM_Node B);
};

#endif
