/** @file */

#ifndef __SPAMM_SORT_H
#define __SPAMM_SORT_H

void
spamm_sort_masked_unsigned_int (const unsigned int length,
    unsigned int *list,
    const unsigned int mask);

void
spamm_sort_float (const unsigned int length,
    float *list);

#endif
