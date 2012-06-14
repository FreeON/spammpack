void
cmultiply_ (float restrict *C, float restrict *A, float restrict *B)
{
	int i, j, k;

	for (i = 0; i < 16; i++) {
		for (j = 0; j < 16; j++) {
			for (k = 0; k < 16; k++)
			{
				C[i+j*16] += A[i+k*16]*B[k+j*16];
			}
		}
	}
}
