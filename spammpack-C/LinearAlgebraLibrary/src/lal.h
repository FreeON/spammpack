/* Interface to library.
 */

typedef struct lal_matrix_t
{
  int M, N;
  double *data;
}
lal_matrix_t;

int
lal_allocate (const int M, const int N, lal_matrix_t **A);

void
lal_free (lal_matrix_t **A);

double
lal_get (const int i, const int j, const lal_matrix_t *A);

void
lal_set (const int i, const int j, const double Aij, lal_matrix_t *A);

int
lal_equals (lal_matrix_t *A, lal_matrix_t *B);

void
lal_zero (lal_matrix_t *A);

void
lal_rand (lal_matrix_t *A);

void
lal_print (lal_matrix_t *A);

lal_matrix_t *
lal_transpose (lal_matrix_t *A);

void
lal_add (lal_matrix_t *A, lal_matrix_t *B, lal_matrix_t *C);

void
lal_dgemm (const char *transA, const char *transB, const int M, const int N,
    const int K, const double alpha, const lal_matrix_t *A, const int lda,
    const lal_matrix_t *B, const int ldb, const double beta, lal_matrix_t *C,
    const int ldc);

void
lal_pointer_to_integer (int *f90_A, struct lal_matrix_t *A);

void
lal_integer_to_pointer (int *f90_A, struct lal_matrix_t **A);
