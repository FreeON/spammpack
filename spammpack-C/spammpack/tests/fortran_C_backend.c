#include <stdlib.h>
#include <stdio.h>

int
c_func_1_ ()
{
  return 1;
}

void
c_func_2_ (int *n)
{
  *n = 1;
}

int64_t
c_func_3_ ()
{
  return 1;
}

void
c_func_4_ (int *n, int64_t *i)
{
  printf("[%s:%i] n = %p = %li\n", __FILE__, __LINE__, n, (int64_t) n);
  *i = (int64_t) n;
}

void
c_func_5_ (int64_t *i)
{
  int *pointer = (int*) *i;
  *pointer = 4;
}
