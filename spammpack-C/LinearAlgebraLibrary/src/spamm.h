struct spamm_t
{
  int M, N;
  struct spamm_node_t *root;
};

struct spamm_node_t
{
  int child_M, child_N;
  struct spamm_node_t *child;
};

int
spamm_dense_index (const int i, const int j, const int stride);

void
spamm_new (struct spamm_t *A);

void
spamm_delete (struct spamm_t *A);

void
spamm_dense_to_spamm (const double *A_dense, struct spamm_t *A);

void
spamm_print_dense (const int M, const int N, const double *A_dense);

void
spamm_print_spamm (const struct spamm_t *A);
