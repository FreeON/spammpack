struct spamm_t
{
  int M, N;
  int M_padded, N_padded;
  int M_block, N_block;
  struct spamm_node_t *root;
};

struct spamm_node_t
{
  int child_M, child_N;
  int block_M, block_N;
  struct spamm_node_t *child[4];
  double *block_dense;
};

void
spamm_log (const char *format, const char *filename, const int linenumber, ...);

int
spamm_dense_index (const int i, const int j, const int stride);

void
spamm_new (const int M, const int N, const int M_block, const int N_block, struct spamm_t *A);

void
spamm_delete (struct spamm_t *A);

void
spamm_dense_to_spamm (const int M, const int N, const int M_block, const int N_block, const double *A_dense, struct spamm_t *A);

double
spamm_get (const int i, const int j, const struct spamm_t *A);

void
spamm_set (const int i, const int j, const double Aij, const struct spamm_t *A);

void
spamm_print_dense (const int M, const int N, const double *A_dense);

void
spamm_print_spamm (const struct spamm_t *A);
