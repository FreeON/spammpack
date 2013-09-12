/** @file
 *
 * The implementation of the BCSR class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "bcsr.h"
#include "logger.h"

#include <charm++.h>
#include <stdio.h>
#include <stdlib.h>

BCSR::BCSR (char *filename)
{
  FILE *fd = NULL;

  if((fd = fopen(filename, "r")) == NULL)
  {
    printf("error opening BCSR file\n");
    exit(0);
  }

  int result;

  if((result = fread(&NSMat,         sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&NAtoms,        sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&numberNonZero, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&numberBlocks,  sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  if((result = fclose(fd)) != 0)
  {
    printf("error closing BCSR file\n");
    exit(0);
  }

  /* For debugging. */
  toStr();
  CkExit();
}

void BCSR::toDense (int *N, double **PDense)
{
}

void BCSR::toStr (void)
{
  printf("BCSR: NSMat = %d\n", NSMat);
  printf("BCSR: NAtoms = %d\n", NAtoms);
  printf("BCSR: numberNonZero = %d\n", numberNonZero);
  printf("BCSR: numberBlocks = %d\n", numberBlocks);
}
