NVCC = nvcc

noinst_PROGRAMS = kernel

kernel_SOURCES = main.c kernel.cu

SUFFIXES = .cu
.cu.$(OBJEXT) :
	$(NVCC) -c -o $@ $<
