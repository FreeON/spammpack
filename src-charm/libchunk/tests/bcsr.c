#include "bcsr.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INFO(format, ...) print_log("INFO", __FILE__, __LINE__, __func__, format, ##__VA_ARGS__)
#define ABORT(format, ...) print_log("ERROR", __FILE__, __LINE__, __func__, format, ##__VA_ARGS__); exit(-1)

/** The linear offset inside a dense non-square matrix block. Column-major
 * order. */
#define BLOCK_INDEX_NONSQUARE(i, j, iLower, jLower, M, N) ((i-iLower)+(j-jLower)*M)

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

/** Print something.
 */
static void
print_log (
    const char *const tag,
    const char *const filename,
    const int line_number,
    const char *const function_name,
    const char *const format,
    ...
    )
{
#define FORMAT_LENGTH 1000

  va_list ap;
  char new_format[FORMAT_LENGTH];

  snprintf(new_format, FORMAT_LENGTH, "[%s - %s:%d (%s)] %s",
      tag, filename, line_number, function_name, format);

  va_start(ap, format);
  vprintf(new_format, ap);
  va_end(ap);
}

/** Get the number of rows of the matrix.
 *
 * @param A The BCSR matrix.
 *
 * @return The number of rows.
 */
int bcsr_get_M (const struct bcsr_t *A)
{
  return A->M;
}

/** Get the number of columns of the matrix.
 *
 * @param A The BCSR matrix.
 *
 * @return The number of columns.
 */
int bcsr_get_N (const struct bcsr_t *A)
{
  return A->N;
}

/** Load a BCSR style matrix from disk. Not the format is the newer
 * stand-alone BCSR.
 *
 * @param filename The filename of the matrix file.
 *
 * @return The newly allocated BCSR matrix.
 */
struct bcsr_t *
bcsr_load (const char *const filename)
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

  if((result = fread(A->blockSize, sizeof(int), A->NAtoms, fd)) != A->NAtoms) { ABORT("read error\n"); }
  if((result = fread(A->offset, sizeof(int), A->NAtoms, fd)) != A->NAtoms) { ABORT("read error\n"); }

  A->rowPointer = calloc(A->NAtoms+1, sizeof(int));
  A->columnPointer = calloc(A->numberBlocks, sizeof(int));
  A->blockPointer = calloc(A->numberBlocks, sizeof(int));
  A->matrix = calloc(A->numberNonZero, sizeof(double));

  double *dummy = calloc(A->numberBlocks, sizeof(double));

  if((result = fread(A->rowPointer, sizeof(int), A->NAtoms+1, fd)) != A->NAtoms+1) { ABORT("read error\n"); }
  if((result = fread(A->columnPointer, sizeof(int), A->numberBlocks, fd)) != A->numberBlocks) { ABORT("read error\n"); }
  if((result = fread(A->blockPointer, sizeof(int), A->numberBlocks, fd)) != A->numberBlocks) { ABORT("read error\n"); }
  if((result = fread(dummy, sizeof(double), A->numberBlocks, fd)) != A->numberBlocks) { ABORT("read error\n"); }
  if((result = fread(A->matrix, sizeof(double), A->numberNonZero, fd)) != A->numberNonZero) { ABORT("read error\n"); }

  free(dummy);

  /* Fortran counts starting at 1, C starts from 0. */
  for(int i = 0; i < A->NAtoms; i++)
  {
    A->offset[i]--;
  }
  for(int i = 0; i < A->NAtoms+1; i++)
  {
    A->rowPointer[i]--;
  }
  for(int i = 0; i < A->numberBlocks; i++)
  {
    A->columnPointer[i]--;
    A->blockPointer[i]--;
  }

  if((result = fclose(fd)) != 0)
  {
    ABORT("error closing BCSR file\n");
  }

  INFO("read BCSR matrix, size %dx%d, %d nonzeros, %1.4f percent nonzero elements\n",
      A->M, A->N, A->numberNonZero, 100*A->numberNonZero/(double) (A->N*A->N));

  return A;
}

/** Convert a BCSR matrix into a dense matrix.
 *
 * @param M [out] The number of rows.
 * @param N [out] The number of columns.
 * @param A_BCSR [in] The BCSR matrix.
 *
 * @return The dense matrix. Free with free().
 */
double *
bcsr_to_dense (
    int *M,
    int *N,
    const struct bcsr_t *A_BCSR
    )
{
  *M = A_BCSR->M;
  *N = A_BCSR->N;

  double *A_dense = calloc(A_BCSR->M*A_BCSR->N, sizeof(double));

  for(int atom = 0; atom < A_BCSR->NAtoms; atom++)
  {
    int MBlock = A_BCSR->blockSize[atom];
    int rowOffset = A_BCSR->offset[atom];

    for(int iRow = A_BCSR->rowPointer[atom]; iRow < A_BCSR->rowPointer[atom+1]; iRow++)
    {
      int NBlock = A_BCSR->blockSize[A_BCSR->columnPointer[iRow]];
      int columnOffset = A_BCSR->offset[A_BCSR->columnPointer[iRow]];

      for(int iSMat = 0; iSMat < A_BCSR->NSMat; iSMat++)
      {
        int pointer = A_BCSR->blockPointer[iRow]+iSMat*MBlock*NBlock;
        int i_block;
        int j_block;

        switch(iSMat+1)
        {
          case 1:
            i_block = rowOffset;
            j_block = columnOffset;
            break;

          case 2:
            i_block = rowOffset;
            j_block = columnOffset+A_BCSR->numberBasisFunctions;
            break;

          case 3:
            i_block = rowOffset+A_BCSR->numberBasisFunctions;
            j_block = columnOffset;
            break;

          case 4:
            i_block = rowOffset+A_BCSR->numberBasisFunctions;
            j_block = columnOffset+A_BCSR->numberBasisFunctions;
            break;

          default:
            ABORT("error\n");
            break;
        }

        for(int i = 0; i < MBlock; i++) {
          for(int j = 0; j < NBlock; j++)
          {
            A_dense[BLOCK_INDEX_NONSQUARE(i+i_block, j+j_block, 0, 0, A_BCSR->M, A_BCSR-N)] =
              A_BCSR->matrix[pointer+i+j*MBlock];
          }
        }
      }
    }
  }

  return A_dense;
}
