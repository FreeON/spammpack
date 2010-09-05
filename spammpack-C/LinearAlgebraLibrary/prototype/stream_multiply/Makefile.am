noinst_PROGRAMS = \
  stream_multiply \
  memory_access \
  Drepper_Follow \
	Drepper_matrix_multiply

AM_LDFLAGS = -all-static

stream_multiply_SOURCES = \
  stream_multiply.c \
  stream_kernel_1.S \
  stream_kernel_2.S \
  stream_kernel_3.S \
  stream_kernel_4.S \
  stream_kernel_5.S

tags :
	ctags *.c *.S
