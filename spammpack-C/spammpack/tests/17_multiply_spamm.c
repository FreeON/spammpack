#include "config.h"
#include <spamm.h>

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TEST_ABS_TOLERANCE 2e-8
#define TEST_REL_TOLERANCE 2e-6

inline unsigned int
matrix_index (const unsigned int i,
    const unsigned int j,
    const unsigned int M,
    const unsigned int N)
{
  return i+j*M;
}

int
main (int argc, char **argv)
{
  int result = 0;

  unsigned int dim;
  unsigned int i[2];
  unsigned int k;

  unsigned int N[] = { 513, 513 };

  unsigned int chunk_tier = 5;

  /* Boolean program parameters. */
  short use_linear_tree = 0;
  short use_sgemm = 0;
  short use_diagonal = 0;
  short verify_result = 0;
  short check_matrices = 0;
  short random_matrix = 1;
  short print_debug = 0;

  float gamma = 1.0;

  double alpha = 1.2;
  double beta = 0.5;

  float alpha_float = alpha;
  float beta_float = beta;

  float tolerance = 0.0;

  double *A_dense;
  double *B_dense;
  double *C_dense;

  float *A_float;
  float *B_float;
  float *C_float;

  struct spamm_matrix_t *A;
  struct spamm_matrix_t *B;
  struct spamm_matrix_t *C;

  unsigned int max_i[] = { 0, 0 };
  unsigned int max_rel_i[] = { 0, 0 };

  double max_diff;
  double max_rel_diff;
  double max_diff_float;

  enum spamm_kernel_t kernel = kernel_standard_SSE;
  struct spamm_timer_t *timer;
  char *timer_string;

  int option_index;
  int parse_result;
  char *short_options = "hk:N:la:b:t:c:rds1g:vx";
  static struct option long_options[] = {
    { "help",       no_argument,        NULL, 'h' },
    { "kernel",     required_argument,  NULL, 'k' },
    { "N",          required_argument,  NULL, 'N' },
    { "linear",     no_argument,        NULL, 'l' },
    { "alpha",      required_argument,  NULL, 'a' },
    { "beta",       required_argument,  NULL, 'b' },
    { "tolerance",  required_argument,  NULL, 't' },
    { "contiguous", required_argument,  NULL, 'c' },
    { "no-random",  no_argument,        NULL, 'r' },
    { "debug",      no_argument,        NULL, 'd' },
    { "sgemm",      no_argument,        NULL, 's' },
    { "diagonal",   no_argument,        NULL, '1' },
    { "gamma",      required_argument,  NULL, 'g' },
    { "verify",     no_argument,        NULL, 'v' },
    { "check",      no_argument,        NULL, 'x' },
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
        printf("{ -a | --alpha } alpha        Set alpha\n");
        printf("{ -b | --beta } beta          Set beta\n");
        printf("{ -t | --tolerance } tau      Multiply with tolerance tau\n");
        printf("{ -c | --contiguous } c       Set contiguous tier to c\n");
        printf("{ -r | --no-random }          Do not create random matrix\n");
        printf("{ -d | --debug }              Print matrices\n");
        printf("{ -s | --sgemm }              Use sgemm\n");
        printf("{ -1 | --diagonal }           Create diagonally dominant matrices\n");
        printf("{ -g | --gamma } g            Set decay for diagonal to g\n");
        printf("{ -v | --verify }             Verify result\n");
        printf("{ -x | --check }              Check matrices\n");
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

      case 'a':
        alpha = strtod(optarg, NULL);
        break;

      case 'b':
        beta = strtod(optarg, NULL);
        break;

      case 't':
        tolerance = strtod(optarg, NULL);
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

      case 'v':
        verify_result = 1;
        break;

      case 'x':
        check_matrices = 1;
        break;

      default:
        printf("unknown option\n");
        break;
    }
  }

  A_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);
  B_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);
  C_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);

  A_float = (float*) malloc(sizeof(float)*N[0]*N[1]);
  B_float = (float*) malloc(sizeof(float)*N[0]*N[1]);
  C_float = (float*) malloc(sizeof(float)*N[0]*N[1]);

  printf("creating random matrices... ");
  fflush(stdout);
  for(i[0] = 0; i[0] < N[0]; i[0]++)
  {
    if(use_diagonal)
    {
      if(random_matrix)
      {
        A_dense[matrix_index(i[0], i[0], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
        B_dense[matrix_index(i[0], i[0], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
        C_dense[matrix_index(i[0], i[0], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
      }

      else
      {
        A_dense[matrix_index(i[0], i[0], N[0], N[1])] = i[0]*N[1]+i[0]+1;
        B_dense[matrix_index(i[0], i[0], N[0], N[1])] = i[0]*N[1]+i[0]+1;
        C_dense[matrix_index(i[0], i[0], N[0], N[1])] = i[0]*N[1]+i[0]+1;
      }
    }

    for(i[1] = 0; i[1] < N[1]; i[1]++)
    {
      if(use_diagonal)
      {
        if(i[0] != i[1])
        {
          A_dense[matrix_index(i[0], i[1], N[0], N[1])] = A_dense[matrix_index(i[0], i[0], N[0], N[1])]*(fabs((float) i[0]-(float) i[1]) > gamma ? expf(-fabsf((float) i[0]-(float) i[1])/gamma) : 1);
          B_dense[matrix_index(i[0], i[1], N[0], N[1])] = B_dense[matrix_index(i[0], i[0], N[0], N[1])]*(fabs((float) i[0]-(float) i[1]) > gamma ? expf(-fabsf((float) i[0]-(float) i[1])/gamma) : 1);
          C_dense[matrix_index(i[0], i[1], N[0], N[1])] = C_dense[matrix_index(i[0], i[0], N[0], N[1])]*(fabs((float) i[0]-(float) i[1]) > gamma ? expf(-fabsf((float) i[0]-(float) i[1])/gamma) : 1);
        }
      }

      else
      {
        if(random_matrix)
        {
          A_dense[matrix_index(i[0], i[1], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
          B_dense[matrix_index(i[0], i[1], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
          C_dense[matrix_index(i[0], i[1], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
        }

        else
        {
          A_dense[matrix_index(i[0], i[1], N[0], N[1])] = i[0]*N[1]+i[1];
          B_dense[matrix_index(i[0], i[1], N[0], N[1])] = i[0]*N[1]+i[1];
          C_dense[matrix_index(i[0], i[1], N[0], N[1])] = i[0]*N[1]+i[1];
        }
      }
    }
  }

  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      A_float[matrix_index(i[0], i[1], N[0], N[1])] = A_dense[i[0]*N[1]+i[1]];
      B_float[matrix_index(i[0], i[1], N[0], N[1])] = B_dense[i[0]*N[1]+i[1]];
      C_float[matrix_index(i[0], i[1], N[0], N[1])] = C_dense[i[0]*N[1]+i[1]];
    }
  }
  printf("done\n");

  printf("creating SpAMM matrices... ");
  fflush(stdout);

  A = spamm_convert_dense_to_spamm(2, N, chunk_tier, use_linear_tree, row_major, A_float);
  B = spamm_convert_dense_to_spamm(2, N, chunk_tier, use_linear_tree, row_major, B_float);
  C = spamm_convert_dense_to_spamm(2, N, chunk_tier, use_linear_tree, row_major, C_float);

  printf("done\n");

  printf("multiply: alpha = %f, beta = %f, tolerance = %f\n", alpha, beta, tolerance);

  //spamm_print_info(A);
  //spamm_print_info(B);
  //spamm_print_info(C);

  if(check_matrices)
  {
    printf("checking SpAMM matrices...\n");
    printf("checking A\n");
    spamm_check(A, 1e-7);
    printf("checking B\n");
    spamm_check(B, 1e-7);
    printf("checking C\n");
    spamm_check(C, 1e-7);
    printf("done\n");
  }

  if(print_debug)
  {
    printf("A_dense = zeros(%u,%u);\n", N[0], N[1]);
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf("A_dense(%u,%u) = %e;\n", i[0]+1, i[1]+1, A_dense[matrix_index(i[0], i[1], N[0], N[1])]);
      }
    }

    printf("A =\n");
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf(" %5.2e", spamm_get(i, A));
      }
      printf("\n");
    }

    printf("B_dense = zeros(%u,%u);\n", N[0], N[1]);
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf("B_dense(%u,%u) = %e;\n", i[0]+1, i[1]+1, B_dense[matrix_index(i[0], i[1], N[0], N[1])]);
      }
    }

    printf("C_dense = zeros(%u,%u);\n", N[0], N[1]);
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf("C_dense(%u,%u) = %e;\n", i[0]+1, i[1]+1, C_dense[matrix_index(i[0], i[1], N[0], N[1])]);
      }
    }

    printf("A =\n");
    spamm_print_tree(A);

    printf("C =\n");
    spamm_print_tree(C);
  }

  if(verify_result)
  {
    printf("multiplying reference... ");
    fflush(stdout);
#ifdef DGEMM
    DGEMM("N", "N", &N[0], &N[0], &N[0], &alpha, A_dense, &N[0], B_dense, &N[0], &beta, C_dense, &N[0]);
#else
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        C_dense[matrix_index(i[0], i[1], N[0], N[1])] *= beta;
        for(k = 0; k < N[0]; k++)
        {
          C_dense[matrix_index(i[0], i[1], N[0], N[1])] += alpha*A_dense[matrix_index(i[0], k, N[0], N[1])]*B_dense[matrix_index(k, i[1], N[0], N[1])];
        }
      }
    }
#endif
    printf("done\n");

    printf("multiplying sgemm... ");
    fflush(stdout);
#ifdef SGEMM
    SGEMM("N", "N", &N[0], &N[0], &N[0], &alpha_float, A_float, &N[0], B_float, &N[0], &beta_float, C_float, &N[0]);
#else
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        C_float[matrix_index(i[0], i[1], N[0], N[1])] *= beta;
        for(k = 0; k < N[0]; k++)
        {
          C_float[matrix_index(i[0], i[1], N[0], N[1])] += alpha*A_float[matrix_index(i[0], k, N[0], N[1])]*B_float[matrix_index(k, i[1], N[0], N[1])];
        }
      }
    }
#endif
    printf("done\n");
  }

  if(print_debug)
  {
    printf("C_ref_dense = zeros(%u,%u);\n", N[0], N[1]);
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf("C_ref_dense(%u,%u) = %e;\n", i[0]+1, i[1]+1, C_dense[matrix_index(i[0], i[1], N[0], N[1])]);
      }
    }
  }

  timer = spamm_timer_new();
  spamm_timer_add_event(0x8000003b, timer);
  spamm_timer_start(timer);
  spamm_multiply(tolerance, alpha, A, B, beta, C, (use_sgemm ? SGEMM : NULL), kernel, NULL);
  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
  printf("timer: %s\n", timer_string);
  free(timer_string);
  spamm_timer_delete(&timer);

  if(print_debug)
  {
    printf("C =\n");
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf(" %5.2e", spamm_get(i, C));
      }
      printf("\n");
    }
  }

  if(check_matrices)
  {
    printf("checking C\n");
    spamm_check(C, 1e-7);
  }

  if(verify_result)
  {
    max_diff = 0;
    max_rel_diff = 0;
    max_diff_float = 0;
    printf("verifying result... ");
    fflush(stdout);
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++) {
        for(k = 0; k < N[0]; k++)
        {
          if(fabs(C_dense[matrix_index(i[0], i[1], N[0], N[1])]-spamm_get(i, C)) > max_diff)
          {
            max_diff = fabs(C_dense[matrix_index(i[0], i[1], N[0], N[1])]-spamm_get(i, C));
            max_i[0] = i[0];
            max_i[1] = i[1];
          }

          if(fabs(C_dense[matrix_index(i[0], i[1], N[0], N[1])]-C_float[matrix_index(i[0], i[1], N[0], N[1])]) > max_diff_float)
          {
            max_diff_float = fabs(C_dense[matrix_index(i[0], i[1], N[0], N[1])]-C_float[matrix_index(i[0], i[1], N[0], N[1])]);
          }

          if(C_dense[matrix_index(i[0], i[1], N[0], N[1])] != 0)
          {
            if(fabs((C_dense[matrix_index(i[0], i[1], N[0], N[1])]-spamm_get(i, C))/C_dense[matrix_index(i[0], i[1], N[0], N[1])]) > max_rel_diff)
            {
              max_rel_diff = fabs((C_dense[matrix_index(i[0], i[1], N[0], N[1])]-spamm_get(i, C))/C_dense[matrix_index(i[0], i[1], N[0], N[1])]);
              max_rel_i[0] = i[0];
              max_rel_i[1] = i[1];
            }
          }
        }
      }
    }
    printf("done\n");

    printf("max diff =       %e, rel. diff = %e, A[%u][%u] = %e, A_reference[%u][%u] = %e\n",
        max_diff,
        (C_dense[matrix_index(max_i[0], max_i[1], N[0], N[1])] != 0.0 ? max_diff/C_dense[matrix_index(max_i[0], max_i[1], N[0], N[1])] : 0.0),
        max_i[0], max_i[1], spamm_get(max_i, C),
        max_i[0], max_i[1], C_dense[matrix_index(max_i[0], max_i[1], N[0], N[1])]);
    printf("max float diff = %e\n", max_diff_float);
    printf("max rel. diff =  %e, diff =      %e, A[%u][%u] = %e, A_reference[%u][%u] = %e\n",
        max_rel_diff,
        fabs(C_dense[matrix_index(max_rel_i[0], max_rel_i[1], N[0], N[1])]-spamm_get(max_rel_i, C)),
        max_rel_i[0], max_rel_i[1], spamm_get(max_rel_i, C),
        max_rel_i[0], max_rel_i[1], C_dense[matrix_index(max_rel_i[0], max_rel_i[1], N[0], N[1])]);

    if(max_diff > TEST_ABS_TOLERANCE && max_rel_diff > TEST_REL_TOLERANCE)
    {
      printf("test failed (abs. tolerance = %e, rel. tolerance = %e)\n", TEST_ABS_TOLERANCE, TEST_REL_TOLERANCE);
      result = -1;
    }
  }

  free(A_dense);
  free(B_dense);
  free(C_dense);

  free(A_float);
  free(B_float);
  free(C_float);

  spamm_delete(&A);
  spamm_delete(&B);
  spamm_delete(&C);

  return result;
}
