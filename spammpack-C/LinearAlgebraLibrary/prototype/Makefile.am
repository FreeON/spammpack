if HAVE_CUDA
cudadirs = cublas cuda_gemm_kernel
else
cudadirs =
endif

if HAVE_BLAS
blasdirs = kernel_comparison
else
blasdirs =
endif

SUBDIRS = $(cudadirs) $(blasdirs) linear_quadtree Z_curve stream_multiply
