#include <spamm.h>

#include <stdio.h>
#include <stdlib.h>

void
sgeev_ (char * jobvl, char *jobvr, int *N, float *A, int *lda, float *wr, float *wi, float *vl, int *ldvl, float *vr, int *ldvr, float *work, int *lwork, int *info);

int
main (int argc, char **argv)
{
  float *A_dense;
  struct spamm_matrix_t *A;

  unsigned int N[2] = { 513, 513 };
  unsigned int i[2];

  unsigned int chunk_tier = 3;
  short use_linear_tree = 0;

  float *wr;
  float *wi;
  float *vl;
  float *vr;
  float *work;
  int ldvl;
  int ldvr;
  int lwork;
  int info;

  float f_min_reference;
  float f_max_reference;

  float f_min;
  float f_max;

  A_dense = malloc(N[0]*N[1]*sizeof(float));

  for(i[0] = 0; i[0] < N[0]; i[0]++) {
    for(i[1] = 0; i[1] < N[1]; i[1]++)
    {
      A_dense[i[0]+i[1]*N[1]] = rand()/(float) RAND_MAX;
    }
  }

  A = spamm_convert_dense_to_spamm(2, N, chunk_tier, use_linear_tree, column_major, A_dense);

  /* Get the eigenvalues. */
  wr = malloc(N[0]*sizeof(float));
  wi = malloc(N[0]*sizeof(float));
  ldvl = 1;
  ldvr = 1;
  vl = malloc(ldvl*sizeof(float));
  vr = malloc(ldvr*sizeof(float));
  lwork = 3*N[0];
  work = malloc(lwork*sizeof(float));
  sgeev_("N", "N", &N[0], A_dense, &N[0], wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);

  for(i[0] = 0; i[0] < N[0]; i[0]++)
  {
    if(i[0] == 0)
    {
      f_min_reference = wr[i[0]];
      f_max_reference = wr[i[0]];
    }

    else
    {
      if(wr[i[0]] < f_min_reference)
      {
        f_min_reference = wr[i[0]];
      }
      if(wr[i[0]] > f_max_reference)
      {
        f_max_reference = wr[i[0]];
      }
    }
  }

  printf("f_reference = [ %f, %f ]\n", f_min_reference, f_max_reference);

  spamm_spectral_bounds(&f_min, &f_max, A);

  printf("f = [ %f, %f ]\n", f_min, f_max);

  /* Free memory. */
  free(work);
  free(vr);
  free(vl);
  free(wi);
  free(wr);
  free(A_dense);
  spamm_delete(&A);
}
