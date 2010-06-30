#include "spamm.h"
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** The maximum number of characters per input line.
 */
#define LINE_MAX 2000

/** Read a matrix from file in MatrixMarket format.
 *
 * @param filename The filename of the matrix to read.
 * @param M_block Number of rows of matrix blocks in output matrix.
 * @param N_block Number of columns of matrix blocks in output matrix.
 * @param M_child Number of rows of children per node in matrix.
 * @param N_child Number of columns of children per node in matrix.
 * @param threshold The matrix threshold below which elements are considered
 *        zero.
 * @param A The output matrix.
 */
void
spamm_read_MM (const char *filename,
    const unsigned int M_block, const unsigned int N_block,
    const unsigned int M_child, const unsigned int N_child,
    const floating_point_t threshold, struct spamm_t *A)
{
  FILE *fd;
  int linenumber;
  int number_nonzero;
  int number_dropped;
  int i, j, M, N;
  floating_point_t Aij;
  char line[LINE_MAX];
  char *token;

  assert(A != NULL);
  assert(filename != NULL);

  fd = fopen(filename, "r");
  if (fd == NULL)
  {
    LOG_FATAL("can not open file %s for reading: %s\n", filename, strerror(errno));
    exit(1);
  }

  /* Reset matrix size. */
  M = -1;
  N = -1;

  linenumber = 0;
  number_nonzero = 0;
  number_dropped = 0;
  while (fgets(line, LINE_MAX, fd) != NULL)
  {
    linenumber++;

    if (linenumber == 1)
    {
      /* Load header line. */
      LOG_DEBUG("header: %s", line);
      continue;
    }

    /* Remove comments. */
    for (i = 0; i < strlen(line); ++i)
    {
      if (line[i] == '%')
      {
        line[i] = '\0';
        break;
      }
    }

    /* Skip empty lines. */
    if (strlen(line) == 0) { continue; }

    if (M < 0 && N < 0)
    {
      /* Figure out matrix size. */
      token = strtok(line, " \t");
      if (token == NULL) { LOG_FATAL("syntax error, line %i\n", linenumber); }
      M = strtol(token, NULL, 10);
      token = strtok(NULL, " \t");
      if (token == NULL) { LOG_FATAL("syntax error, line %i\n", linenumber); }
      N = strtol(token, NULL, 10);
      spamm_new(M, N, M_block, N_block, M_child, N_child, threshold, A);
      continue;
    }

    /* Load elements. */
    token = strtok(line, " \t");
    if (token == NULL) { LOG_FATAL("syntax error, line %i\n", linenumber); }
    i = strtol(token, NULL, 10)-1;
    token = strtok(NULL, " \t");
    if (token == NULL) { LOG_FATAL("syntax error, line %i\n", linenumber); }
    j = strtol(token, NULL, 10)-1;
    token = strtok(NULL, " \t");
    if (token == NULL) { LOG_FATAL("syntax error, line %i\n", linenumber); }
    Aij = strtod(token, NULL);
    if (spamm_set(i, j, Aij, A) != SPAMM_RESULT_OK) { number_dropped++; }
    number_nonzero++;
  }

  LOG_DEBUG("loaded %i nonzero elements, %i dropped below threshold of %e\n", number_nonzero, number_dropped, threshold);
  fclose(fd);
}
