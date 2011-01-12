#ifndef __SPAMM_CONVERT_H
#define __SPAMM_CONVERT_H

enum spamm_dense_type_t
{
  row_major, column_major
};

struct spamm_t *
spamm_convert_dense_to_spamm (const unsigned int M, const unsigned int N,
    const enum spamm_dense_type_t type, float *A_dense);

#endif
