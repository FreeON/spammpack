/** @file
 *
 * The header file for the BCSRInfo program.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __BCSRINFO_H
#define __BCSRINFO_H

#include "bcsrinfo.decl.h"

/** The main class. */
class BCSRInfo : public CBase_BCSRInfo
{
  public:

    BCSRInfo (CkArgMsg *msg);
};

#endif
