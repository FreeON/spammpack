#include "spamm.h"
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_MAX 2000

void
spamm_read_MM (const char *filename, const int M_block, const int N_block,
    const int M_child, const int N_child, const float_t threshold,
    struct spamm_t *A)
{
  FILE *fd;
  int linenumber;
  int i, j, M, N;
  float_t Aij;
  char line[LINE_MAX];
  char *token;

  assert(A != NULL);
  assert(filename != NULL);

  fd = fopen(filename, "r");
  if (fd == NULL)
  {
    spamm_log("can not open file %s for reading: %s\n", __FILE__, __LINE__, filename, strerror(errno));
    exit(1);
  }

  /* Reset matrix size. */
  M = -1;
  N = -1;

  linenumber = 0;
  while (fgets(line, LINE_MAX, fd) != NULL)
  {
    linenumber++;

    if (linenumber == 1)
    {
      /* Load header line. */
      spamm_log("header: %s", __FILE__, __LINE__, line);
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
      if (token == NULL) { spamm_log("syntax error, line %i\n", __FILE__, __LINE__, linenumber); }
      M = strtol(token, NULL, 10);
      token = strtok(NULL, " \t");
      if (token == NULL) { spamm_log("syntax error, line %i\n", __FILE__, __LINE__, linenumber); }
      N = strtol(token, NULL, 10);
      spamm_new(M, N, M_block, N_block, M_child, N_child, threshold, A);
      continue;
    }

    /* Load elements. */
    token = strtok(line, " \t");
    if (token == NULL) { spamm_log("syntax error, line %i\n", __FILE__, __LINE__, linenumber); }
    i = strtol(token, NULL, 10)-1;
    token = strtok(NULL, " \t");
    if (token == NULL) { spamm_log("syntax error, line %i\n", __FILE__, __LINE__, linenumber); }
    j = strtol(token, NULL, 10)-1;
    token = strtok(NULL, " \t");
    if (token == NULL) { spamm_log("syntax error, line %i\n", __FILE__, __LINE__, linenumber); }
    Aij = strtod(token, NULL);
    spamm_set(i, j, Aij, A);
  }

  fclose(fd);
}
