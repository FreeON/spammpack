#include "config.h"

#include "test.h"

#include <spamm.h>

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define REL_TOLERANCE 1e-8

#define TEST_ABS_TOLERANCE 2e-8
#define TEST_REL_TOLERANCE 2e-6

unsigned int
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

#if !defined(DGEMM) && !defined(SGEMM)
  unsigned int k;
#endif

  unsigned int N[] = { 513, 513 };

  unsigned int chunk_tier = 5;

  /* Boolean program parameters. */
  short use_linear_tree = 0;
  short use_sgemm = 0;
  short use_diagonal = 0;
  short verify_result = 1;
  short check_matrices = 1;
  short random_matrix = 1;
  short print_debug = 0;
  short null_C = 0;

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
  unsigned int max_float_i[] = { 0, 0 };

  double max_diff;
  double max_rel_diff;
  double max_diff_float;

  struct spamm_timer_t *timer;
  char *timer_string;

  time_t random_seed = time(NULL);

  enum matrix_t matrix_type = full;

  double flops;

  int option_index;
  int parse_result;
  char *short_options = "hN:la:b:t:c:rds1g:vxnm:2:";
  static struct option long_options[] = {
    { "help",         no_argument,        NULL, 'h' },
    { "N",            required_argument,  NULL, 'N' },
    { "linear",       no_argument,        NULL, 'l' },
    { "alpha",        required_argument,  NULL, 'a' },
    { "beta",         required_argument,  NULL, 'b' },
    { "tolerance",    required_argument,  NULL, 't' },
    { "chunk",        required_argument,  NULL, 'c' },
    { "no-random",    no_argument,        NULL, 'r' },
    { "debug",        no_argument,        NULL, 'd' },
    { "sgemm",        no_argument,        NULL, 's' },
    { "diagonal",     no_argument,        NULL, '1' },
    { "gamma",        required_argument,  NULL, 'g' },
    { "verify",       no_argument,        NULL, 'v' },
    { "check",        no_argument,        NULL, 'x' },
    { "nullC",        no_argument,        NULL, 'n' },
    { "matrix-type",  required_argument,  NULL, 'm' },
    { "seed",         required_argument,  NULL, '2' },
    { NULL,           0,                  NULL,  0  }
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
        printf("{ -N | --N } N                Set N\n");
        printf("{ -l | --linear }             Use a linear tier\n");
        printf("{ -a | --alpha } alpha        Set alpha\n");
        printf("{ -b | --beta } beta          Set beta\n");
        printf("{ -t | --tolerance } tau      Multiply with tolerance tau\n");
        printf("{ -c | --chunk } c            Set chunk tier to c\n");
        printf("{ -r | --no-random }          Do not create random matrix\n");
        printf("{ -d | --debug }              Print matrices\n");
        printf("{ -s | --sgemm }              Use sgemm\n");
        printf("{ -1 | --diagonal }           Create diagonally dominant matrices\n");
        printf("{ -g | --gamma } g            Set decay for diagonal to g\n");
        printf("{ -v | --verify }             Verify result\n");
        printf("{ -x | --check }              Check matrices\n");
        printf("{ -n | --nullC }              Initialize C as empty matrix\n");
        printf("{ -m | --matrix-type } TYPE   Matrix type: %s\n", print_matrix_types());
        printf("{ --seed } SEED               Seed random number generator with SEED\n");
        exit(0);
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
        alpha_float = alpha;
        break;

      case 'b':
        beta = strtod(optarg, NULL);
        beta_float = beta;
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
        verify_result = (verify_result+1)%2;
        break;

      case 'x':
        check_matrices = (check_matrices+1)%2;
        break;

      case 'n':
        null_C = 1;
        break;

      case 'm':
        matrix_type = parse_matrix_type(optarg);
        break;

      case '2':
        random_seed = strtol(optarg, NULL, 10);
        break;

      default:
        printf("unknown option\n");
        break;
    }
  }

  /* Initialize random number generator. */
  srand(random_seed);

  A_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);
  B_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);
  C_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);

  A_float = (float*) malloc(sizeof(float)*N[0]*N[1]);
  B_float = (float*) malloc(sizeof(float)*N[0]*N[1]);
  C_float = (float*) malloc(sizeof(float)*N[0]*N[1]);

  printf("creating random matrices (%ux%u)... ", N[0], N[1]);
  fflush(stdout);
  for(i[0] = 0; i[0] < N[0]; i[0]++)
  {
    if(use_diagonal)
    {
      if(random_matrix)
      {
        A_dense[matrix_index(i[0], i[0], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
        B_dense[matrix_index(i[0], i[0], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
        if(!null_C)
        {
          C_dense[matrix_index(i[0], i[0], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
        }
      }

      else
      {
        A_dense[matrix_index(i[0], i[0], N[0], N[1])] = i[0]*N[1]+i[0]+1;
        B_dense[matrix_index(i[0], i[0], N[0], N[1])] = i[0]*N[1]+i[0]+1;
        if(!null_C)
        {
          C_dense[matrix_index(i[0], i[0], N[0], N[1])] = i[0]*N[1]+i[0]+1;
        }
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
          if(!null_C)
          {
            C_dense[matrix_index(i[0], i[1], N[0], N[1])] = C_dense[matrix_index(i[0], i[0], N[0], N[1])]*(fabs((float) i[0]-(float) i[1]) > gamma ? expf(-fabsf((float) i[0]-(float) i[1])/gamma) : 1);
          }
        }
      }

      else
      {
        if(random_matrix)
        {
          A_dense[matrix_index(i[0], i[1], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
          B_dense[matrix_index(i[0], i[1], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
          if(!null_C)
          {
            C_dense[matrix_index(i[0], i[1], N[0], N[1])] = (float) rand()/(float) RAND_MAX;
          }
        }

        else
        {
          A_dense[matrix_index(i[0], i[1], N[0], N[1])] = i[0]*N[1]+i[1];
          B_dense[matrix_index(i[0], i[1], N[0], N[1])] = i[0]*N[1]+i[1];
          if(!null_C)
          {
            C_dense[matrix_index(i[0], i[1], N[0], N[1])] = i[0]*N[1]+i[1];
          }
        }
      }
    }
  }

  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      A_float[matrix_index(i[0], i[1], N[0], N[1])] = A_dense[matrix_index(i[0], i[1], N[0], N[1])];
      B_float[matrix_index(i[0], i[1], N[0], N[1])] = B_dense[matrix_index(i[0], i[1], N[0], N[1])];
      C_float[matrix_index(i[0], i[1], N[0], N[1])] = C_dense[matrix_index(i[0], i[1], N[0], N[1])];
    }
  }
  printf("done\n");

  printf("creating SpAMM matrices... ");
  fflush(stdout);

  A = spamm_convert_dense_to_spamm(2, N, chunk_tier, use_linear_tree, column_major, A_float);
  B = spamm_convert_dense_to_spamm(2, N, chunk_tier, use_linear_tree, column_major, B_float);
  if(null_C)
  {
    C = spamm_new(2, N, chunk_tier, use_linear_tree);
  }

  else
  {
    C = spamm_convert_dense_to_spamm(2, N, chunk_tier, use_linear_tree, column_major, C_float);
  }

  printf("done\n");

  printf("multiply: alpha = %f, beta = %f, tolerance = %f, chunk_tier = %u, use_linear_tree = %u\n",
      alpha, beta, tolerance, chunk_tier, use_linear_tree);

  if(check_matrices)
  {
    printf("checking A... ");
    if(spamm_check(A, REL_TOLERANCE) != SPAMM_OK)
    {
      SPAMM_FATAL("failed\n");
    }

    else
    {
      printf("ok\n");
    }

    printf("checking B... ");
    if(spamm_check(B, REL_TOLERANCE) != SPAMM_OK)
    {
      SPAMM_FATAL("failed\n");
    }

    else
    {
      printf("ok\n");
    }

    printf("checking C... ");
    if(spamm_check(C, REL_TOLERANCE) != SPAMM_OK)
    {
      SPAMM_FATAL("failed\n");
    }

    else
    {
      printf("ok\n");
    }
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

    printf("B =\n");
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf(" %5.2e", spamm_get(i, B));
      }
      printf("\n");
    }

    printf("C_dense = zeros(%u,%u);\n", N[0], N[1]);
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf("C_dense(%u,%u) = %e;\n", i[0]+1, i[1]+1, C_dense[matrix_index(i[0], i[1], N[0], N[1])]);
      }
    }

    printf("C =\n");
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf(" %5.2e", spamm_get(i, C));
      }
      printf("\n");
    }

    printf("A (tree)\n");
    spamm_print_tree(A);

    printf("B (tree)\n");
    spamm_print_tree(B);

    printf("C (tree)\n");
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
    printf("C_ref_float = zeros(%u,%u);\n", N[0], N[1]);
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf("C_ref_float(%u,%u) = %e;\n", i[0]+1, i[1]+1, C_float[matrix_index(i[0], i[1], N[0], N[1])]);
      }
    }

    printf("C_ref_dense = zeros(%u,%u);\n", N[0], N[1]);
    for(i[0] = 0; i[0] < N[0]; i[0]++) {
      for(i[1] = 0; i[1] < N[1]; i[1]++)
      {
        printf("C_ref_dense(%u,%u) = %e;\n", i[0]+1, i[1]+1, C_dense[matrix_index(i[0], i[1], N[0], N[1])]);
      }
    }
  }

  printf("multiplying SpAMM... ");
  fflush(stdout);

  timer = spamm_timer_new();
  spamm_timer_add_event(0x8000003b, timer);
  spamm_timer_start(timer);
  spamm_multiply(tolerance, alpha, A, B, beta, C, (use_sgemm ? SGEMM : NULL), &flops);
  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
  printf("%s\n", timer_string);
  free(timer_string);
  spamm_timer_delete(&timer);

  if(print_debug)
  {
    printf("C (tree)\n");
    spamm_print_tree(C);

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
    printf("checking C... ");
    if(spamm_check(C, REL_TOLERANCE) != SPAMM_OK)
    {
      SPAMM_FATAL("failed\n");
    }

    else
    {
      printf("ok\n");
    }
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
        if(fabs(C_dense[matrix_index(i[0], i[1], N[0], N[1])]-spamm_get(i, C)) > max_diff)
        {
          max_diff = fabs(C_dense[matrix_index(i[0], i[1], N[0], N[1])]-spamm_get(i, C));
          max_i[0] = i[0];
          max_i[1] = i[1];
        }

        if(fabs(C_dense[matrix_index(i[0], i[1], N[0], N[1])]-C_float[matrix_index(i[0], i[1], N[0], N[1])]) > max_diff_float)
        {
          max_diff_float = fabs(C_dense[matrix_index(i[0], i[1], N[0], N[1])]-C_float[matrix_index(i[0], i[1], N[0], N[1])]);
          max_float_i[0] = i[0];
          max_float_i[1] = i[1];
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
    printf("done\n");

    printf("max float diff =       %e, C_dense[%u][%u] = %e, C_float[%u][%u] = %e\n",
        max_diff_float,
        max_float_i[0], max_float_i[1], C_dense[matrix_index(max_float_i[0], max_float_i[1], N[0], N[1])],
        max_float_i[0], max_float_i[1], C_float[matrix_index(max_float_i[0], max_float_i[1], N[0], N[1])]);
    printf("max SpAMM diff =       %e, rel. diff = %e, A[%u][%u] = %e, A_reference[%u][%u] = %e\n",
        max_diff,
        (C_dense[matrix_index(max_i[0], max_i[1], N[0], N[1])] != 0.0 ? max_diff/C_dense[matrix_index(max_i[0], max_i[1], N[0], N[1])] : 0.0),
        max_i[0], max_i[1], spamm_get(max_i, C),
        max_i[0], max_i[1], C_dense[matrix_index(max_i[0], max_i[1], N[0], N[1])]);
    printf("max SpAMM rel. diff =  %e, diff =      %e, A[%u][%u] = %e, A_reference[%u][%u] = %e\n",
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
