/** @file
 *
 * The implementation of the BCSR class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "bcsr.h"
#include "logger.h"

#include <stdio.h>
#include <stdlib.h>

BCSR::BCSR (char *filename)
{
  FILE *fd = NULL;

  if((fd = fopen(filename, "r")) == NULL)
  {
    ABORT("error opening BCSR file \"%s\"\n", filename);
  }

  int result;
  int record_length;

  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  if((result = fread(&NSMat, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&NAtoms, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&numberNonZero, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&numberBlocks, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  rowPointer = new int[NAtoms+1];
  colPointer = new int[numberBlocks];
  blockPointer = new int[numberBlocks];
  matrix = new double[numberNonZero];

  double *dummy = new double[numberBlocks];

  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(rowPointer, sizeof(int), NAtoms+1, fd)) != NAtoms+1) { ABORT("read error\n"); }
  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(colPointer, sizeof(int), numberBlocks, fd)) != numberBlocks) { ABORT("read error\n"); }
  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(blockPointer, sizeof(int), numberBlocks, fd)) != numberBlocks) { ABORT("read error\n"); }
  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(dummy, sizeof(double), numberBlocks, fd)) != numberBlocks) { ABORT("read error\n"); }
  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(matrix, sizeof(double), numberNonZero, fd)) != numberNonZero) { ABORT("read error\n"); }
  if((result = fread(&record_length, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  delete[] dummy;

  if((result = fclose(fd)) != 0)
  {
    ABORT("error closing BCSR file\n");
  }
}

void BCSR::toDense (int *N, double **PDense)
{
}

void BCSR::toStr (void)
{
  printf("BCSR: NSMat         = %d\n", NSMat);
  printf("BCSR: NAtoms        = %d\n", NAtoms);
  printf("BCSR: numberNonZero = %d\n", numberNonZero);
  printf("BCSR: numberBlocks  = %d\n", numberBlocks);
}
