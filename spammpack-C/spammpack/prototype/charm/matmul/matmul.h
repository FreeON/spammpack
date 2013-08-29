/** @file
 *
 * The header file for the matmul program.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __MATMUL_H
#define __MATMUL_H

#include "matmul.decl.h"

/** The main entry method. */
class Main : public CBase_Main
{
  public:

    Main (CkArgMsg *msg);
    void run (int N, int blocksize, int numberIterations, double tolerance,
        int matrixType, double decayConstant, bool verify, double verifyTolerance,
        bool loadBalance, bool printPEMap);
};

#endif
