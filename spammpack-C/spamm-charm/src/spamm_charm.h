/** @file
 *
 * The header file for the spamm program.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __SPAMM_CHARM_H
#define __SPAMM_CHARM_H

#include "spamm_charm.decl.h"

/** The main class. */
class SpAMM_Charm : public CBase_SpAMM_Charm
{
  public:

    SpAMM_Charm (CkArgMsg *msg);
    void run (int N, int blocksize, int N_basic, int numberIterations,
        double tolerance, int matrixType, double decayConstant, int operation,
        bool verify, double verifyTolerance, bool loadBalance, int initialPE,
        bool alignPEs, bool printPEMap);
    void runSP2 (int lengthFockianFilename, char *fockianFilename,
        int lengthDensityFilename, char *densityFilename, int Ne, int N,
        int blocksize, int N_basic, int maxIterations, double tolerance,
        bool loadBalance, int initialPE, bool alignPEs, bool printPEMap,
        int lengthPEMap, char *filenamePEMap, double F_min, double F_max);
};

#endif
