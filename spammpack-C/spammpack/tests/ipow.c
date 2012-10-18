#include <spamm.h>

#include <stdio.h>

unsigned int
get_reference (unsigned int b, unsigned int n)
{
  unsigned int i;
  unsigned int reference = b;

  for (i = 0; i < (n-1); i++)
  {
    reference *= b;
  }

  return reference;
}

int
main ()
{
  int result = 0;

  unsigned int b1 = 2;
  unsigned int b2 = 3;

  unsigned int n = 8;

  if (get_reference(b1, n) != ipow(b1, n))
  {
    printf("wrong result: found %u, should have been %u\n", ipow(b1, n), get_reference(b1, n));
    result = SPAMM_ERROR;
  }

  if (get_reference(b2, n) != ipow(b2, n))
  {
    printf("wrong result: found %u, should have been %u\n", ipow(b2, n), get_reference(b2, n));
    result = SPAMM_ERROR;
  }

  return result;
}
