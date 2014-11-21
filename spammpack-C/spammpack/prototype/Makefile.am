#if HAVE_CUDA
#cudadirs = cublas cuda_gemm_kernel
#endif
#
#if HAVE_BLAS
#blasdirs = kernel_comparison
#endif

if HAVE_CHARM
charmdirs = charm
endif

#if PROTOTYPES
#otherdirs = linear_quadtree Z_curve stream_multiply hash_tables
#endif

SUBDIRS = $(cudadirs) $(blasdirs) $(charmdirs) $(otherdirs)
