#ifndef __BCSR_H
#define __BCSR_H

struct bcsr_t;

int bcsr_get_M (const struct bcsr_t *A);
int bcsr_get_N (const struct bcsr_t *A);

struct bcsr_t *
bcsr_load (const char *const filename);

double *
bcsr_to_dense (
    int *M,
    int *N,
    const struct bcsr_t *A_BCSR
    );

#endif
