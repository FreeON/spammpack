#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#define ABORT(format, ...) printf(format, ##__VA_ARGS__); exit(-1)

struct bcsr_t
{
  /** The number of rows. */
  int M;

  /** The number of columns. */
  int N;

  /** The number of spin matrices. */
  int NSMat;

  /** The number of atoms. */
  int NAtoms;

  /** The number of non-zero elements. */
  int numberNonZero;

  /** The number of blocks. */
  int numberBlocks;

  /** The number of basis functions. */
  int numberBasisFunctions;

  /** The block sizes per atom. */
  int *blockSize;

  /** The offsets. */
  int *offset;

  /** The array of row indices. */
  int *rowPointer;

  /** The array of column indices. */
  int *columnPointer;

  /** The array of block indices. */
  int *blockPointer;

  /** The non-zero matrix elements. */
  double *matrix;
};

/** Load a BCSR style matrix from disk. Not the format is the newer
 * stand-alone BCSR.
 *
 * @param filename The filename of the matrix file.
 *
 * @return The newly allocated BCSR matrix.
 */
struct bcsr_t *
load_BCSR (const char *const filename)
{
  FILE *fd = NULL;

  if((fd = fopen(filename, "r")) == NULL)
  {
    ABORT("error opening BCSR file \"%s\"\n", filename);
  }

  int result;
  struct bcsr_t *A = calloc(1, sizeof(struct bcsr_t));

  if((result = fread(&A->NSMat, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&A->NAtoms, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&A->numberNonZero, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&A->numberBlocks, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&A->numberBasisFunctions, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&A->M, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&A->N, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  A->blockSize = calloc(A->NAtoms, sizeof(int));
  A->offset = calloc(A->NAtoms, sizeof(int));

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

  INFO("read BCSR matrix, size %dx%d, %d nonzeros, %1.4f percent nonzero elements\n",
      M, N, numberNonZero, 100*numberNonZero/(double) (N*N));
}

int
main (int argc, char **argv)
{
  char *F_filename = NULL;

  int c;
  const char short_options[] = "h";
  const struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
  };

  while((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
  {
    switch(c)
    {
      default:
        printf("illegal command line argument\n");
        exit(-1);
        break;
    };
  }

  if(F_filename == NULL)
  {
    printf("missing Fockian matrix file\n");
    exit(-1);
  }

  return 0;
}
