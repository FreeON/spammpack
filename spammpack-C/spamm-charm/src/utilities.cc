/** @file
 *
 * The implementation of the utility functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "index.h"
#include "logger.h"
#include "utilities.h"

#include <bitset>
#include <charm++.h>
#include <errno.h>
#include <sstream>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

/** Print a dense matrix in matrix market format.
 *
 * @param N The matrix size.
 * @param A The matrix.
 * @param format The format string for a message that is printed in front of
 * the matrix. See printf() for details.
 */
void printDense (int N, double *A, const char *const format, ...)
{
  std::ostringstream o;
  o.setf(std::ios::scientific);

  va_list ap;
  const int message_length = 2000;
  char message[message_length];
  char numberBuffer[40];

  va_start(ap, format);
  vsnprintf(message, message_length, format, ap);
  va_end(ap);

  if(N <= 32)
  {
    o << message << " = [" << std::endl;
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++)
      {
        snprintf(numberBuffer, 40, "% 1.3e", A[BLOCK_INDEX(i, j, 0, 0, N)]);
        o << " " << numberBuffer;
      }
      o << std::endl;
    }
    o << "]" << std::endl;
  }

  else
  {
    o << "%%%%MatrixMarket matrix coordinate double general" << std::endl;
    o << "%% " << message << " = " << N << " x " << N << " --> " << N*N << " elements" << std::endl;
    o << N << " " << N << " " << N*N << std::endl;
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++)
      {
        snprintf(numberBuffer, 40, "% e", A[BLOCK_INDEX(i, j, 0, 0, N)]);
        o << i+1 << " " << j+1 << " " << numberBuffer << std::endl;
      }
    }
  }

  CkPrintf(o.str().c_str());
}

/** Convert an integer to a binary as a string.
 *
 * @param i The integer.
 *
 * @return The string represenation.
 */
std::string toBinary (unsigned int i)
{
  std::string bitString = std::bitset<8*sizeof(unsigned int)>(i).to_string();
  while(bitString[0] == '0')
  {
    bitString.erase(0, 1);
  }
  return bitString;
}

/** Convert some number of bytes into a human readable form. For instance, n =
 * 1024 would result in "1 kiB".
 *
 * @param n The number of bytes to convert.
 *
 * @return The string representation of the human readable form.
 */
std::string humanReadableSize (size_t n)
{
  std::ostringstream o;
  double size = 0;
  const int numberUnits = 5;
  std::string unit[numberUnits] = { "B", "kiB", "MiB", "GiB", "TiB" };

  for(int i = 0; i < numberUnits; i++)
  {
    size /= 1024.;

    if(n/1024 > 0)
    {
      size += (n%1024);
      n /= 1024;
    }

    else
    {
      size += n;
      o << size << " " << unit[i];
      return o.str();
    }
  }
  ABORT("can not convert this number\n");
  return NULL;
}

/** Load a matrix in coordinate format, i.e.
 *
 * i j A_{ij}
 *
 * from file and allocate a dense array.
 *
 * @param filename The file name.
 * @param N [out] The matrix size.
 * @param ADense [out] The dense matrix.
 */
void loadCoordinateFile (char *filename, int *N, double **ADense)
{
  FILE *fd;

  if((fd = fopen(filename, "r")) == NULL)
  {
    ABORT("error opening density file \"%s\"\n", filename);
  }

  INFO("reading density matrix from \"%s\"\n", filename);

  *N = -1;
  char linebuffer[2000];
  int i, j;
  double Aij;
  int linenumber = 0;
  int result;
  while(fgets(linebuffer, 2000, fd) == linebuffer)
  {
    linenumber++;
    if((result = sscanf(linebuffer, "%d %d %le\n", &i, &j, &Aij)) == 3)
    {
      DEBUG("read %d %d %e\n", i, j, Aij);
      if(i > *N) { *N = i; }
      if(j > *N) { *N = j; }
    }

    else
    {
      break;
    }
  }

  if(result == EOF)
  {
    if(ferror(fd) != 0)
    {
      ABORT("error reading file: %s\n", strerror(errno));
    }
  }

  if(result == 0)
  {
    while(linebuffer[strlen(linebuffer)-1] == '\n')
    {
      linebuffer[strlen(linebuffer)-1] = '\0';
    }

    ABORT("syntax error, line %d: \"%s\"\n", linenumber, linebuffer);
  }

  if(*N < 1)
  {
    ABORT("could not read coordinate file\n");
  }

  INFO("reading %dx%d matrix\n", *N, *N);

  rewind(fd);

  *ADense = new double[(*N)*(*N)];

  while(fgets(linebuffer, 2000, fd) == linebuffer)
  {
    sscanf(linebuffer, "%d %d %le\n", &i, &j, &Aij);
    (*ADense)[BLOCK_INDEX(i-1, j-1, 0, 0, *N)] = Aij;
  }

  fclose(fd);
}

/** Get the spectral bounds of the matrix by using the Gershgorin circle
 * theorem.
 *
 * Estimate spectral bounds via Gersgorin approximation, @f$ \left[
 * F_{min}-F_{max} \right] @f$.
 *
 * In detail:
 * @f[
 *   R_{i} = \sum_{j \neq i} \left| F_{ij} \right|
 * @f]
 * @f[
 *   F_{\mathrm{max}} = \max_{i} \left\{ F_{ii} + R_{i} \right\}
 * @f]
 * @f[
 *   F_{\mathrm{min}} = \min_{i} \left\{ F_{ii} - R_{i} \right\}
 * @f]
 *
 * @param method The method to use. method = 0 is Gershgorin; method = 1 is
 * full eigensolve.
 * @param minBound [out] The lower bound.
 * @param maxBound [out] The upper bound.
 * @param N The matrix size.
 * @param A The dense matrix.
 */
void getSpectralBounds (int method, double *minBound, double *maxBound, int N, double *A)
{
  for(int i = 0; i < N; i++)
  {
    double R = 0;
    for(int j = 0; j < N; j++)
    {
      if(i != j)
      {
        R += fabs(A[BLOCK_INDEX(i, j, 0, 0, N)]);
      }
    }
    if(i == 0)
    {
      *minBound = A[BLOCK_INDEX(i, i, 0, 0, N)] - R;
      *maxBound = A[BLOCK_INDEX(i, i, 0, 0, N)] + R;
    }

    else
    {
      if(*minBound > A[BLOCK_INDEX(i, i, 0, 0, N)] - R)
      {
        *minBound = A[BLOCK_INDEX(i, i, 0, 0, N)] - R;
      }

      if(*maxBound < A[BLOCK_INDEX(i, i, 0, 0, N)] + R)
      {
        *maxBound = A[BLOCK_INDEX(i, i, 0, 0, N)] + R;
      }
    }
  }
}
