/** @file
 *
 * The implementation of the chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"

#define CONCAT2(chunk_type, name) chunk ## _ ## chunk_type ## _ ## name
#define CONCAT(chunk_type, name) CONCAT2(chunk_type, name)
#define FUNC(name) CONCAT(CHUNK_IMPLEMENTATION, name)

#include "chunk.h"
#include "chunk_tiled.h"
#include "chunk_tree.h"

#include <stdio.h>

size_t
chunk_sizeof (const int N_chunk, const int N_basic)
{
  return FUNC(sizeof)(N_chunk, N_basic);
}

void *
chunk_alloc (const int N_chunk,
    const int N_basic,
    const int N,
    const unsigned int i_lower,
    const unsigned int j_lower)
{
  return FUNC(alloc)(N_chunk, N_basic, N, i_lower, j_lower);
}

void
chunk_set (void *const chunk, const double *const A)
{
  FUNC(set)(chunk, A);
}

int
chunk_get_N_chunk (void *const chunk)
{
  return FUNC(get_N_chunk)(chunk);
}

int
chunk_get_N_basic (void *const chunk)
{
  return FUNC(get_N_basic)(chunk);
}

double
chunk_get_norm_2 (const void *const chunk)
{
  return FUNC(get_norm_2)(chunk);
}

void
chunk_print (const void *const chunk,
    const char *const format, ...)
{
  FUNC(print)(chunk, "can not pass the format...\n");
}

void
chunk_add (const double alpha, void *const A,
    const double beta, const void *const B)
{
  FUNC(add)(alpha, A, beta, B);
}

void
chunk_multiply (const double tolerance,
    const void *const A,
    const void *const B,
    void *const C,
    const short symbolic_only)
{
  FUNC(multiply)(tolerance, A, B, C, symbolic_only);
}

double
chunk_trace (const void *const chunk)
{
  return FUNC(trace)(chunk);
}

void
chunk_scale (const double alpha, void *const chunk)
{
  FUNC(scale)(alpha, chunk);
}

void
chunk_add_identity (const double alpha, void *const chunk)
{
  FUNC(add_identity)(alpha, chunk);
}

double *
chunk_to_dense (const void *const chunk)
{
  return FUNC(to_dense)(chunk);
}

void
chunk_delete (void **const chunk)
{
  FUNC(delete)(chunk);
}

size_t
chunk_get_complexity (const void *const chunk)
{
  return FUNC(get_complexity)(chunk);
}
