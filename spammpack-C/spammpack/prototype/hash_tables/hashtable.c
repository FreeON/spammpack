#include <spamm_timer.h>
#include <getopt.h>
#include <glib.h>
#include <stdio.h>
#include <stdlib.h>

unsigned int
get_linear_index (unsigned int i, unsigned int j)
{
  unsigned int bitmask = 1;
  unsigned int index = 0;

  while (i != 0 || j != 0)
  {
    if ((j & 1) == 1)
    {
      index |= bitmask;
    }
    bitmask <<= 1;

    if ((i & 1) == 1)
    {
      index |= bitmask;
    }
    bitmask <<= 1;

    if (i != 0) { i >>= 1; }
    if (j != 0) { j >>= 1; }
  }

  return index;
}

int
main (int argc, char **argv)
{
  unsigned int N = 1; /* The matrix size. */
  unsigned int N_read = 1; /* How many values to read from the hashtable. */
  short N_read_set = 0;

  unsigned int i, j, n;

  int parse_result;
  struct option long_options[] =
  {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "read", required_argument, NULL, 'r' },
    { 0, 0, 0, 0 }
  };
  char *short_options = "hN:r:";

  unsigned int *index = NULL;
  gpointer value;

  GHashTable *hashtable;
  GList *key_list;
  GList *key_list_item;
  unsigned int key_list_length = 0;
  unsigned int *keys = NULL;
  unsigned int temp;

  struct spamm_timer_t *timer = NULL;

  /* Read command line options. */
  while ((parse_result = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
  {
    switch (parse_result)
    {
      case 'h':
        printf("Usage:\n");
        printf("\n");
        printf("{ --help | -h}     This help\n");
        printf("-N N               Test for NxN matrix\n");
        printf("--read N           Read N values randomly from hashtable\n");
        return 0;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'r':
        N_read = strtol(optarg, NULL, 10);
        N_read_set = 1;
        break;

      default:
        printf("unknown option\n");
        return 1;
    }
  }

  if (optind < argc) {
    printf("non-option ARGV-elements: ");
    while (optind < argc)
    {
      printf("%s ", argv[optind++]);
    }
    printf("\n");
    return 1;
  }

  /* Allocate new hashtable. */
  hashtable = g_hash_table_new(g_int_hash, g_int_equal);

  /* Intialize timer. */
  timer = spamm_timer_new(papi_total_cycles);

  /* Fill hashtable with matrix indices. */
  printf("filling hashtable with %ux%u matrix, i.e. %d keys... ", N, N, N*N);
  fflush(NULL);
  spamm_timer_start(timer);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      index = (unsigned int*) malloc(sizeof(unsigned int));
      *index = get_linear_index(i, j);
      g_hash_table_insert(hashtable, index, &N);
    }
  }
  spamm_timer_stop(timer);
  printf("%llu timer units\n", spamm_timer_get(timer));

  /* Get all keys. */
  key_list = g_hash_table_get_keys(hashtable);

  /* Copy keys to array. */
  key_list_length = g_list_length(key_list);
  keys = (unsigned int*) malloc(sizeof(unsigned int)*key_list_length);
  for (i = 0, key_list_item = g_list_first(key_list); key_list_item != NULL; key_list_item = g_list_next(key_list_item))
  {
    keys[i++] = *((unsigned int*) key_list_item->data);
  }
  if (i != key_list_length)
  {
    printf("i is not key_list_length, i = %u, key_list_length = %u\n", i, key_list_length);
    return 1;
  }

  /* Randomize key order. */
  printf("randomizing keys... ");
  fflush(NULL);
  for (i = 0; i < key_list_length-1; i++) {
    for (n = 0; n < 100; n++)
    {
      j = (unsigned int) (rand()/(double) RAND_MAX*(key_list_length-i-1)+i+1);
      if (j < i+1 || j >= key_list_length)
      {
        printf("illegal j value\n");
        return 1;
      }

      temp = keys[i];
      keys[i] = keys[j];
      keys[j] = temp;
    }
  }
  printf("done\n");

  /* Set N_read in case it wasn't on the command line. */
  if (N_read_set == 0)
  {
    N_read = key_list_length;
  }

  /* For comparison, read random elements from keys. This test is supposed to
   * benchmark the lower bound on reading random keys, i.e. each key read
   * should be O(1), but there might exist a small N dependence because of
   * cache or TLB miss issues. */
  printf("reading %u random keys from array... ", N_read);
  fflush(NULL);
  spamm_timer_start(timer);
  for (i = 0; i < N_read; i++)
  {
    j = (unsigned int) (rand()/(double) RAND_MAX*key_list_length);
    if (keys[j] == 0)
    {
      printf("found 0 key... ");
      fflush(NULL);
    }
  }
  spamm_timer_stop(timer);
  printf("%llu timer units, %e timer units per key read\n", spamm_timer_get(timer), spamm_timer_get(timer)/(double) N_read);

  /* Read randomly from hashtable. */
  printf("reading %u random keys... ", N_read);
  fflush(NULL);
  spamm_timer_start(timer);
  for (i = 0; i < N_read; i++)
  {
    j = keys[i%key_list_length];
    value = g_hash_table_lookup(hashtable, &j);
    if (value == NULL)
    {
      printf("strange, encountered NULL value for key %u\n", j);
      return 1;
    }
  }
  spamm_timer_stop(timer);
  printf("%llu timer units, %e timer units per key read\n", spamm_timer_get(timer), spamm_timer_get(timer)/(double) N_read);

  return 0;
}
