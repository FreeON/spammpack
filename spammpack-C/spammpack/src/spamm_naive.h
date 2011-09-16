/** @file */

#ifndef __SPAMM_NAIVE_H
#define __SPAMM_NAIVE_H

#include "spamm_timer.h"
#include "spamm_types.h"

struct spamm_naive_t *
spamm_naive_new (const unsigned int M, const unsigned int N, const unsigned int blocksize);

struct spamm_naive_node_t *
spamm_naive_new_node (const unsigned int tier, const unsigned int blocksize);

double
spamm_naive_get (const unsigned int i, const unsigned int j, const struct spamm_naive_t *A);

void
spamm_naive_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_naive_t *A);

void
spamm_naive_multiply (const float tolerance,
    const float alpha, struct spamm_naive_t *A, struct spamm_naive_t *B,
    const float beta, struct spamm_naive_t *C,
    struct spamm_timer_t *timer,
    void (*sgemm) ());

void
spamm_naive_multiply_scalar (const float beta, struct spamm_naive_node_t *node);

void
spamm_naive_multiply_matrix (const float tolerance,
    const float alpha,
    struct spamm_naive_node_t *node_A,
    struct spamm_naive_node_t *node_B,
    struct spamm_naive_node_t **node_C,
    struct spamm_timer_t *timer,
    void (*sgemm) (char *, char *, int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *));

void
spamm_naive_print (const struct spamm_naive_t *A);

#endif
