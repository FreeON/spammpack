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

void printDense (int N, double *A, const char *const format, ...);

std::string toBinary (unsigned int i);

std::string humanReadableSize (unsigned long n);

#endif