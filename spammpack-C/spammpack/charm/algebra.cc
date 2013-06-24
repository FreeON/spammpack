#include "utilities.h"

void matmul (const int chunksize, const float *const A, const float *const B, float *const C)
{
  LOG_DEBUG("multiplying blocks\n");
  for(int i = 0; i < chunksize; i++) {
    for(int j = 0; j < chunksize; j++)
    {
      for(int k = 0; k < chunksize; k++)
      {
        C[i+j*chunksize] += A[i+k*chunksize]*B[k+j*chunksize];
      }
    }
  }
}
