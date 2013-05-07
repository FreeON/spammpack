/** @file
 *
 * Generate matrices for the tests.
 */

#include "test.h"

#include <assert.h>
#include <math.h>
#include <spamm.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

const char *matrix_type_name[] = {
  "full",
  "diagonally-banded",
  "sparse-random"
};

/** Return a matrix type name. */
const char *const
get_matrix_type_name (const enum matrix_t matrix_type)
{
  assert(matrix_type >= 0 && matrix_type < NUMBER_MATRIX_TYPES);

  return matrix_type_name[matrix_type];
}

/** Print the available matrix types.
 *
 * @return An newly allocated string that contains a list of available matrix
 * types. This string has to be free'ed in the caller.
 */
char *
print_matrix_types ()
{
  char *matrix_type_string;
  int string_length;
  int i;

  for(i = 0, string_length = strlen("{ "); i < NUMBER_MATRIX_TYPES; i++)
  {
    string_length += strlen(matrix_type_name[i]) + strlen(", ");
  }
  string_length += strlen(" }") + 1;

  matrix_type_string = calloc(string_length, sizeof(char));
  strncat(matrix_type_string, "{ ", string_length-1);
  for(i = 0; i < NUMBER_MATRIX_TYPES; i++)
  {
    strncat(matrix_type_string, matrix_type_name[i], string_length-strlen(matrix_type_string));
    if(i+1 < NUMBER_MATRIX_TYPES)
    {
      strncat(matrix_type_string, ",", string_length-strlen(matrix_type_string));
    }
    strncat(matrix_type_string, " ", string_length-strlen(matrix_type_string));
  }
  strncat(matrix_type_string, "}", string_length-strlen(matrix_type_string));

  return matrix_type_string;
}

/** Parse a string to a matrix type.
 *
 * @param type_name The name of the matrix type.
 *
 * @return The matrix type.
 */
enum matrix_t
parse_matrix_type (const char *const type_name)
{
  int i;
  enum matrix_t matrix_type;
  short result = SPAMM_ERROR;

  for(i = 0; i < NUMBER_MATRIX_TYPES; i++)
  {
    if(strcasecmp(type_name, matrix_type_name[i]) == 0)
    {
      matrix_type = i;
      result = SPAMM_OK;
      break;
    }
  }
  if(result != SPAMM_OK)
  {
    SPAMM_FATAL("parse error, unknown matrix type \"%s\"\n", type_name);
  }

  return matrix_type;
}

/** Generate a matrix shape.
 *
 * @param number_dimensions The number of dimensions.
 * @param is_square Whether the matrix is square shaped.
 *
 * @return The shape of the matrix.
 */
unsigned int *
generate_shape (const unsigned int number_dimensions,
    const short is_square)
{
  unsigned int *N;
  unsigned int dim;

  N = calloc(number_dimensions, sizeof(unsigned int));

  for(dim = 0; dim < number_dimensions; dim++)
  {
    N[dim] = 150+(int) ((0.5-(float) rand()/(float) RAND_MAX)*30);
  }

  return N;
}
