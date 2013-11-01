/** @file
 *
 * The header file for the spamm program.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __SPAMM_H
#define __SPAMM_H

#include "spamm.decl.h"

/** The main class. */
class SpAMM : public CBase_SpAMM
{
  public:

    SpAMM (CkArgMsg *msg);
    void run (int N, int blocksize, int numberIterations, double tolerance,
        int matrixType, double decayConstant, int operation, bool verify,
        double verifyTolerance, bool loadBalance, int initialPE,
        bool alignPEs, bool printPEMap);
    void runSP2 (int lengthFilename, char *filename, int Ne, int blocksize,
        int maxIterations, double tolerance, bool loadBalance, int initialPE,
        bool alignPEs, bool printPEMap, int lengthPEMap, char *filenamePEMap);
};

#endif
