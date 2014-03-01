#include <getopt.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#define ABORT(format, ...) print_log("ERROR", __FILE__, __LINE__, __func__, format, ##__VA_ARGS__); exit(-1)
#define INFO(format, ...) print_log("INFO", __FILE__, __LINE__, __func__, format, ##__VA_ARGS__)

#define FORMAT_LENGTH 1000

void
print_log (const char *const tag,
    const char *const filename,
    const int line_number,
    const char *const function_name,
    const char *const format,
    ...)
{
  va_list ap;
  char new_format[FORMAT_LENGTH];

  snprintf(new_format, FORMAT_LENGTH, "[%s - %s:%d (%s)] %s",
      tag, filename, line_number, function_name, format);

  va_start(ap, format);
  vprintf(new_format, ap);
  va_end(ap);
}

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
}

void
print_help_and_exit (void)
{
  printf("Usage: chunk_SP2 [options]\n");
  printf("\n");
  printf("{ -h | --help }     This help\n");
  exit(0);
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
      case 'h':
        print_help_and_exit();
        break;

      default:
        printf("illegal command line argument\n");
        exit(-1);
        break;
    };
  }

  if(F_filename == NULL)
  {
    ABORT("missing Fockian matrix file\n");
  }

  return 0;
}
