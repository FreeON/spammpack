#include "config.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>

int
work (const int max)
{
  int i, j, k;

  k = 0;

  for(i = 0; i < max; i++)
  {
    for(j = 0; j < 200000; j++)
    {
      k++;
      if(k > 1000) k = 0;
    }
  }

  return k;
}

int
tree (const int max, const int tier)
{
  int i ;

  if(tier < 5)
  {
    for(i = 0; i < 2; i++)
    {
#pragma omp task untied if(tier < CHUNK_TREE_MAX_TIER)
      {
        tree(max/2, tier+1);
      }
    }
#pragma omp taskwait
  }
  else
  {
    return work(max);
  }
}

int
main (int argc, char **argv)
{
  int k;
  int num_threads;
  const int max = 1 << 16;

  printf("work = %d\n", max);
#pragma omp parallel
  {
#pragma omp master
    {
#ifdef _OPENMP
      num_threads = omp_get_num_threads();
      printf("starting task (using %d threads)...\n", num_threads);
#else
      num_threads = 1;
      printf("starting task (serial)...\n");
#endif
#pragma omp task untied
      {
        k = tree(max, 0);
      }
#pragma omp taskwait
    }
  }
  printf("k = %d\n", k);

  return 0;
}
