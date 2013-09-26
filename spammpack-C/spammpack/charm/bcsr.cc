/** @file
 *
 * The implementation of the BCSR class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "bcsr.h"
#include "logger.h"
#include "index.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** The constructor.
 *
 * @param filename The filename of the BCSR file.
 */
BCSR::BCSR (char *filename)
{
  FILE *fd = NULL;

  if((fd = fopen(filename, "r")) == NULL)
  {
    ABORT("error opening BCSR file \"%s\"\n", filename);
  }

  int result;

  if((result = fread(&NSMat, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&NAtoms, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&numberNonZero, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&numberBlocks, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&numberBasisFunctions, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&M, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&N, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  blockSize = new int[NAtoms];
  offset = new int[NAtoms];

  if((result = fread(blockSize, sizeof(int), NAtoms, fd)) != NAtoms) { ABORT("read error\n"); }
  if((result = fread(offset, sizeof(int), NAtoms, fd)) != NAtoms) { ABORT("read error\n"); }

  rowPointer = new int[NAtoms+1];
  columnPointer = new int[numberBlocks];
  blockPointer = new int[numberBlocks];
  matrix = new double[numberNonZero];

  double *dummy = new double[numberBlocks];

  if((result = fread(rowPointer, sizeof(int), NAtoms+1, fd)) != NAtoms+1) { ABORT("read error\n"); }
  if((result = fread(columnPointer, sizeof(int), numberBlocks, fd)) != numberBlocks) { ABORT("read error\n"); }
  if((result = fread(blockPointer, sizeof(int), numberBlocks, fd)) != numberBlocks) { ABORT("read error\n"); }
  if((result = fread(dummy, sizeof(double), numberBlocks, fd)) != numberBlocks) { ABORT("read error\n"); }
  if((result = fread(matrix, sizeof(double), numberNonZero, fd)) != numberNonZero) { ABORT("read error\n"); }

  delete[] dummy;

  /* Fortran counts starting at 1, C starts from 0. */
  for(int i = 0; i < NAtoms; i++)
  {
    offset[i]--;
  }
  for(int i = 0; i < NAtoms+1; i++)
  {
    rowPointer[i]--;
  }
  for(int i = 0; i < numberBlocks; i++)
  {
    columnPointer[i]--;
    blockPointer[i]--;
  }

  if((result = fclose(fd)) != 0)
  {
    ABORT("error closing BCSR file\n");
  }
}

/** Convert a BCSR matrix into a dense matrix.
 *
 * @param M [out] The number of rows.
 * @param N [out] The number of columns.
 * @param ADense [out] The dense matrix.
 */
void BCSR::toDense (int *M, int *N, double **ADense)
{
  *M = this->M;
  *N = this->N;

  *ADense = new double[this->M*this->N];
  memset(*ADense, 0, sizeof(double)*this->M*this->N);

  for(int atom = 0; atom < NAtoms; atom++)
  {
    int MAtom = blockSize[atom];
    int rowOffset = offset[atom];

    for(int iRow = rowPointer[atom]; iRow < rowPointer[atom+1]; iRow++)
    {
      int NAtom = blockSize[columnPointer[iRow]];
      int columnOffset = offset[columnPointer[iRow]];

      for(int iSMat = 0; iSMat < NSMat; iSMat++)
      {
        int pointer = blockPointer[iRow]+iSMat*MAtom*NAtom;
        int i_block;
        int j_block;

        switch(NSMat)
        {
          case 1:
            i_block = rowOffset;
            j_block = columnOffset;
            break;

          case 2:
            i_block = rowOffset;
            j_block = columnOffset+numberBasisFunctions;
            break;

          case 3:
            i_block = rowOffset+numberBasisFunctions;
            j_block = columnOffset;
            break;

          case 4:
            i_block = rowOffset+numberBasisFunctions;
            j_block = columnOffset+numberBasisFunctions;
            break;

          default:
            ABORT("error\n");
            break;
        }

        for(int i = 0; i < MAtom; i++) {
          for(int j = 0; j < NAtom; j++)
          {
            (*ADense)[BLOCK_INDEX_NONSQUARE(i+i_block, j+j_block, 0, 0, this->M, this-N)] = matrix[pointer+i+j*MAtom];
          }
        }
      }
    }
  }
}

/** Print some information on the BCSR matrix.
 */
void BCSR::toStr (void)
{
  printf("BCSR: M             = %d\n", M);
  printf("BCSR: N             = %d\n", N);
  printf("BCSR: NBasF         = %d\n", numberBasisFunctions);
  printf("BCSR: NSMat         = %d\n", NSMat);
  printf("BCSR: NAtoms        = %d\n", NAtoms);
  printf("BCSR: block sizes   = {");
  for(int i = 0; i < NAtoms; i++)
  {
    printf(" %d", blockSize[i]);
  }
  printf(" }\n");
  printf("BCSR: offset        = {");
  for(int i = 0; i < NAtoms; i++)
  {
    printf(" %d", offset[i]);
  }
  printf(" }\n");
  printf("BCSR: numberNonZero = %d\n", numberNonZero);
  printf("BCSR: numberBlocks  = %d\n", numberBlocks);
}
