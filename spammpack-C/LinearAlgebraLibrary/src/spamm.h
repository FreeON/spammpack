struct spamm_t
{
  /* Matrix dimensions. */
  int M, N;

  /* The padded dimensions. */
  int M_padded, N_padded;

  /* The dimensions of the data blocks at the leaf level. */
  int M_block, N_block;

  /* The root node. */
  struct spamm_node_t *root;
};

struct spamm_node_t
{
  /* The dimensions covered in the padded matrix in this node. The indices are
   * meant as [M_lower, M_upper[, i.e. the upper limit is not included in the
   * interval. */
  int M_lower, M_upper;
  int N_lower, N_upper;

  /* The dimensions of the data blocks at the leaf level. */
  int M_block, N_block;

  /* At the non-block level, pointers to the children nodes. */
  struct spamm_node_t *child[4];

  /* At the block level, the dense matrix data. */
  double *block_dense;
};

void
spamm_log (const char *format, const char *filename, const int linenumber, ...);

int
spamm_dense_index (const int i, const int j, const int stride);

void
spamm_new (const int M, const int N, const int M_block, const int N_block, struct spamm_t *A);

void
spamm_new_node (struct spamm_node_t **node);

void
spamm_delete (struct spamm_t *A);

void
spamm_dense_to_spamm (const int M, const int N, const int M_block, const int N_block, const double *A_dense, struct spamm_t *A);

double
spamm_get (const int i, const int j, const struct spamm_t *A);

void
spamm_set (const int i, const int j, const double Aij, struct spamm_t *A);

void
spamm_print_dense (const int M, const int N, const double *A_dense);

void
spamm_print_spamm (const struct spamm_t *A);

void
spamm_print_node (const struct spamm_node_t *node);
