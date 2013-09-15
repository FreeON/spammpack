/** @file
 *
 * A program to print out some information on a BCSR file.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "bcsr.h"
#include "bcsrinfo.h"

BCSRInfo::BCSRInfo (CkArgMsg *args)
{
  BCSR A((char*) "A.OrthoF");
  A.toStr();
  CkExit();
}

#include "bcsrinfo.def.h"
