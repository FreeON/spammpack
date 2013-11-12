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

/** The Block class. */
class Block
{
  private:

    /** The blocksize. We are using square matrices. */
    int blocksize;

    /** The dense matrix. */
    double *block;

    /** The norm. */
    double norm;

    /** The square of the norm. */
    double norm_2;

    void updateNorm (void);

  public:

    Block (void);
    Block & operator= (const Block &rhs);
    void pup (PUP::er &p);
    void set (const int blocksize, const double *const A);
    int getBlocksize (void);
    double getNorm (void);
    double * toDense (void);
    void scale (const double alpha);
    void multiply (Block A, Block B);
    void add (const double alpha, const double beta, const Block A);
    double trace (void);
    void addIdentity (const int numberRows, const int blocksize,
        const double alpha);
    void print (const char *const format, ...);
};

#endif
