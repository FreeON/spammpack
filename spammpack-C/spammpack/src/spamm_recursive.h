/** @file */

#ifndef __SPAMM_NAIVE_H
#define __SPAMM_NAIVE_H

#include "spamm_timer.h"
#include "spamm_types.h"

struct spamm_recursive_t *
spamm_recursive_new (const unsigned int M, const unsigned int N, const unsigned int blocksize);

struct spamm_recursive_node_t *
spamm_recursive_new_node (const unsigned int tier, const unsigned int blocksize);

double
spamm_recursive_get (const unsigned int i, const unsigned int j, const struct spamm_recursive_t *A);

void
spamm_recursive_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_recursive_t *A);

void
spamm_recursive_multiply (const float tolerance,
    const float alpha, struct spamm_recursive_t *A, struct spamm_recursive_t *B,
    const float beta, struct spamm_recursive_t *C,
    struct spamm_timer_t *timer,
    void (*sgemm) ());

void
spamm_recursive_print (const struct spamm_recursive_t *A);

void spamm_recursive_sgemm (char * transA, char * transB,
    int *M, int *N, int *K,
    float *alpha, float *A, int *LDA, float *B, int *LDB,
    float *beta, float *C, int *LDC);

#endif
