/** @file
 *
 * The headers of the utility functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __UTILITIES_H
#define __UTILITIES_H

#include <string>

/** The minimum function.
 *
 * @param a The first value.
 * @param b The second value.
 *
 * @return The smaller value of a and b.
 */
#define MIN(a, b) ((a) < (b) ? (a) : (b))

void printDense (int N, double *A, const char *const format, ...);

std::string toBinary (unsigned int i);

std::string humanReadableSize (size_t n);

void loadCoordinateFile (char *filename, int *N, double **ADense);

void getSpectralBounds (int method, double *minBound, double *maxBound, int N, double *A);

#endif
