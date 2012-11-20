#include "config.h"
#include <spamm.h>

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define VERIFY_RESULT

#define TEST_ABS_TOLERANCE 2e-8
#define TEST_REL_TOLERANCE 2e-6

int
main (int argc, char **argv)
{
  int result = 0;

  unsigned int dim;
  unsigned int i[2];
  unsigned int k;

  unsigned int N[] = { 513, 513 };

  unsigned int chunk_tier = 5;

  short use_linear_tree = 0;
  short use_sgemm = 0;
  short use_diagonal = 0;

  float gamma = 1.0;

  double alpha = 1.2;
  double beta = 0.5;

  short random_matrix = 1;
  short print_debug = 0;

  float tolerance = 0.0;

  double *A_dense;
  double *B_dense;
  double *C_dense;

  struct spamm_matrix_t *A;
  struct spamm_matrix_t *B;
  struct spamm_matrix_t *C;

  unsigned int max_i[] = { 0, 0 };
  unsigned int max_rel_i[] = { 0, 0 };

  double max_diff;
  double max_rel_diff;

  enum spamm_kernel_t kernel = kernel_standard_SSE;
  struct spamm_timer_t *timer;

  int option_index;
  int parse_result;
  char *short_options = "hk:N:lc:rds1g:";
  static struct option long_options[] = {
    { "help",       no_argument,        NULL, 'h' },
    { "kernel",     required_argument,  NULL, 'k' },
    { "N",          required_argument,  NULL, 'N' },
    { "linear",     no_argument,        NULL, 'l' },
    { "contiguous", required_argument,  NULL, 'c' },
    { "no-random",  no_argument,        NULL, 'r' },
    { "debug",      no_argument,        NULL, 'd' },
    { "sgemm",      no_argument,        NULL, 's' },
    { "diagonal",   no_argument,        NULL, '1' },
    { "gamma",      required_argument,  NULL, 'g' },
    { NULL,         0,                  NULL,  0  }
  };

  while(1)
  {
    parse_result = getopt_long(argc, argv, short_options, long_options, &option_index);

    if(parse_result == -1) { break; }

    switch(parse_result)
    {
      case 'h':
        printf("Usage:\n");
        printf("\n");
        printf("{ -k | --kernel } kernel      Use the kernel\n");
        printf("{ -N | --N } N                Set N\n");
        printf("{ -l | --linear }             Use a linear tier\n");
        printf("{ -c | --contiguous } c       Set contiguous tier to c\n");
        printf("{ -r | --no-random }          Do not create random matrix\n");
        printf("{ -d | --debug }              Print matrices\n");
        printf("{ -s | --sgemm }              Use sgemm\n");
        printf("{ -1 | --diagonal }           Create diagonally dominant matrices\n");
        printf("{ -g | --gamma } g            Set decay for diagonal to g\n");
        exit(0);
        break;

      case 'k':
        kernel = spamm_kernel_get_kernel(optarg);
        break;

      case 'N':
        for(dim = 0; dim < 2; dim++)
        {
          N[dim] = strtol(optarg, NULL, 10);
        }
        break;

      case 'l':
        use_linear_tree = 1;
        break;

      case 'c':
        chunk_tier = strtol(optarg, NULL, 10);
        break;

      case 'r':
        random_matrix = 0;
        break;

      case 'd':
        print_debug = 1;
        break;

      case 's':
        use_sgemm = 1;
        break;

      case '1':
        use_diagonal = 1;
        break;

      case 'g':
        gamma = strtof(optarg, NULL);
        break;

      default:
        printf("unknown option\n");
        break;
    }
  }

  A_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);
  B_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);
  C_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);

  A = spamm_new(2, N, chunk_tier, use_linear_tree);
  B = spamm_new(2, N, chunk_tier, use_linear_tree);
  C = spamm_new(2, N, chunk_tier, use_linear_tree);

  for(i[0] = 0; i[0] < N[0]; i[0]++)
  {
    if(use_diagonal)
    {
      if(random_matrix)
      {
        A_dense[i[0]*N[1]+i[0]] = (float) rand()/(float) RAND_MAX;
        B_dense[i[0]*N[1]+i[0]] = (float) rand()/(float) RAND_MAX;
        C_dense[i[0]*N[1]+i[0]] = (float) rand()/(float) RAND_MAX;
      }

      else
      {
        A_dense[i[0]*N[1]+i[0]] = i[0]*N[1]+i[0]+1;
        B_dense[i[0]*N[1]+i[0]] = i[0]*N[1]+i[0]+1;
        C_dense[i[0]*N[1]+i[0]] = i[0]*N[1]+i[0]+1;
      }
    }

    for(i[1] = 0; i[1] < N[1]; i[1]++)
    {
      if(use_diagonal)
      {
        if(i[0] != i[1])
        {
          A_dense[i[0]*N[1]+i[1]] = A_dense[i[0]*N[1]+i[0]]*(fabs((float) i[0]-(float) i[1]) > gamma ? expf(-fabsf((float) i[0]-(float) i[1])/gamma) : 1);
          B_dense[i[0]*N[1]+i[1]] = B_dense[i[0]*N[1]+i[0]]*(fabs((float) i[0]-(float) i[1]) > gamma ? expf(-fabsf((float) i[0]-(float) i[1])/gamma) : 1);
          C_dense[i[0]*N[1]+i[1]] = C_dense[i[0]*N[1]+i[0]]*(fabs((float) i[0]-(float) i[1]) > gamma ? expf(-fabsf((float) i[0]-(float) i[1])/gamma) : 1);
        }
      }

      else
      {
        if(random_matrix)
        {
          A_dense[i[0]*N[1]+i[1]] = (float) rand()/(float) RAND_MAX;
          B_dense[i[0]*N[1]+i[1]] = (float) rand()/(float) RAND_MAX;
          C_dense[i[0]*N[1]+i[1]] = (float) rand()/(float) RAND_MAX;
        }

        else
        {
          A_dense[i[0]*N[1]+i[1]] = i[0]*N[1]+i[1];
          B_dense[i[0]*N[1]+i[1]] = i[0]*N[1]+i[1];
          C_dense[i[0]*N[1]+i[1]] = i[0]*N[1]+i[1];
        }
      }

      spamm_set(i, A_dense[i[0]*N[1]+i[1]], A);
      spamm_set(i, B_dense[i[0]*N[1]+i[1]], B);
      spamm_set(i, C_dense[i[0]*N[1]+i[1]], C);
    }
  }

  spamm_print_info(A);
  spamm_print_info(B);
  spamm_print_info(C);

  printf("checking SpAMM matrices... ");
  fflush(stdout);
  spamm_check(A, 1e-7);
  spamm_check(B, 1e-7);
  spamm_check(C, 1e-7);
  printf("done\n");

  if(print_debug)
  {
    printf("A_dense =\n");
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf(" %5.1f", A_dense[i[0]*N[1]+i[1]]);
      }
      printf("\n");
    }

    printf("A =\n");
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf(" %5.1f", spamm_get(i, A));
      }
      printf("\n");
    }

    printf("B_dense =\n");
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf(" %5.1f", B_dense[i[0]*N[1]+i[1]]);
      }
      printf("\n");
    }

    printf("C_dense =\n");
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf(" %5.1f", C_dense[i[0]*N[1]+i[1]]);
      }
      printf("\n");
    }

    printf("A =\n");
    spamm_print_tree(A);
  }

  printf("multiplying reference... ");
  fflush(stdout);
  for(i[0] = 0; i[0] < N[0]; i[0]++) {
    for(i[1] = 0; i[1] < N[1]; i[1]++)
    {
      C_dense[i[0]*N[1]+i[1]] *= beta;
    }
  }

  for(i[0] = 0; i[0] < N[0]; i[0]++) {
    for(i[1] = 0; i[1] < N[1]; i[1]++) {
      for(k = 0; k < N[0]; k++)
      {
        C_dense[i[0]*N[1]+i[1]] += alpha*A_dense[i[0]*N[1]+k]*B_dense[k*N[1]+i[1]];
      }
    }
  }
  printf("done\n");

  if(print_debug)
  {
    printf("C_dense =\n");
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf(" %5.1f", C_dense[i[0]*N[1]+i[1]]);
      }
      printf("\n");
    }
  }

  timer = spamm_timer_new();
  spamm_timer_add_event(0x8000003b, timer);
  spamm_multiply(tolerance, alpha, A, B, beta, C, timer, (use_sgemm ? SGEMM : NULL), kernel, NULL);
  spamm_timer_delete(&timer);

  if(print_debug)
  {
    printf("C =\n");
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf(" %5.1f", spamm_get(i, C));
      }
      printf("\n");
    }
  }

#ifdef VERIFY_RESULT
  max_diff = 0;
  max_rel_diff = 0;
  printf("verifying result... ");
  fflush(stdout);
  for(i[0] = 0; i[0] < N[0]; i[0]++) {
    for(i[1] = 0; i[1] < N[1]; i[1]++) {
      for(k = 0; k < N[0]; k++)
      {
        if(fabs(C_dense[i[0]*N[1]+i[1]]-spamm_get(i, C)) > max_diff)
        {
          max_diff = fabs(C_dense[i[0]*N[1]+i[1]]-spamm_get(i, C));
          max_i[0] = i[0];
          max_i[1] = i[1];
        }

        if(C_dense[i[0]*N[1]+i[1]] != 0)
        {
          if(fabs((C_dense[i[0]*N[1]+i[1]]-spamm_get(i, C))/C_dense[i[0]*N[1]+i[1]]) > max_rel_diff)
          {
            max_rel_diff = fabs((C_dense[i[0]*N[1]+i[1]]-spamm_get(i, C))/C_dense[i[0]*N[1]+i[1]]);
            max_rel_i[0] = i[0];
            max_rel_i[1] = i[1];
          }
        }
      }
    }
  }
  printf("done\n");

  printf("max diff =      %e, rel. diff = %e, A[%u][%u] = %e, A_reference[%u][%u] = %e\n",
      max_diff,
      (C_dense[max_i[0]*N[1]+max_i[1]] != 0.0 ? max_diff/C_dense[max_i[0]*N[1]+max_i[1]] : 0.0),
      max_i[0], max_i[1], spamm_get(max_i, C),
      max_i[0], max_i[1], C_dense[max_i[0]*N[1]+max_i[1]]);
  printf("max rel. diff = %e, diff =      %e, A[%u][%u] = %e, A_reference[%u][%u] = %e\n",
      max_rel_diff,
      fabs(C_dense[max_rel_i[0]*N[1]+max_rel_i[1]]-spamm_get(max_rel_i, C)),
      max_rel_i[0], max_rel_i[1], spamm_get(max_rel_i, C),
      max_rel_i[0], max_rel_i[1], C_dense[max_rel_i[0]*N[1]+max_rel_i[1]]);

  if(max_diff > TEST_ABS_TOLERANCE && max_rel_diff > TEST_REL_TOLERANCE)
  {
    printf("test failed (abs. tolerance = %e, rel. tolerance = %e)\n", TEST_ABS_TOLERANCE, TEST_REL_TOLERANCE);
    result = -1;
  }
#endif

  free(A_dense);
  free(B_dense);
  free(C_dense);

  spamm_delete(&A);
  spamm_delete(&B);
  spamm_delete(&C);

  return result;
}
