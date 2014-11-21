#include <stdio.h>

void
print_bits (int i)
{
  int index;

  printf("% i = 0b", i);
  for (index = 8; index >= 0; index--)
  {
    if (i & (1 << index)) { printf("1"); }
    else                  { printf("0"); }
  }
  printf("\n");
}

int
main ()
{
  print_bits(16);
  print_bits(-16);
  print_bits(30);
  print_bits((30 & (-16)) + 16);
}
