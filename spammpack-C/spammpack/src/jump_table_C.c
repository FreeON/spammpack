#include <stdio.h>
#include <stdlib.h>

void
jump_table ()
{
  int i = rand()/(float) RAND_MAX*100;

  switch (i)
  {
    case 0: printf("at case 0\n"); break;
    case 1: printf("at case 1\n"); break;
    case 2: printf("at case 2\n"); break;
    case 3: printf("at case 3\n"); break;
    case 4: printf("at case 4\n"); break;
    case 5: printf("at case 5\n"); break;
    case 6: printf("at case 6\n"); break;
    case 7: printf("at case 7\n"); break;
    case 8: printf("at case 8\n"); break;
    case 9: printf("at case 9\n"); break;
    case 10: printf("at case 10\n"); break;
    case 11: printf("at case 11\n"); break;
    case 12: printf("at case 12\n"); break;
    case 13: printf("at case 13\n"); break;
    case 14: printf("at case 14\n"); break;
    case 15: printf("at case 15\n"); break;
    case 16: printf("at case 16\n"); break;
    case 17: printf("at case 17\n"); break;
    case 18: printf("at case 18\n"); break;
    case 19: printf("at case 19\n"); break;

    default:
      printf("unknown\n");
      break;
  }
}
