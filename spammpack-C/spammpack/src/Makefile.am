ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

lib_LTLIBRARIES = libspamm.la

libspamm_la_SOURCES = \
  $(matrix_tree_sources) \
  $(matrix_kernel_sources)

matrix_tree_sources = \
  spamm.h \
  spamm_config.h \
  spamm_add.c \
  spamm_allocate.c \
  spamm_block_index.c \
  spamm_delete.c \
  spamm_delete_node.c \
  spamm_dense_index.c \
  spamm_dense_to_spamm.c \
  spamm_free.c \
  spamm_get.c \
  spamm_int_to_binary.c \
  spamm_linear_index.c \
  spamm_linear_to_coordinates.c \
  spamm_logger.c \
  spamm_multiply.c \
  spamm_multiply_scalar.c \
  spamm_new.c \
  spamm_new_childnode.c \
  spamm_new_node.c \
  spamm_number_nonzero.c \
  spamm_print_dense.c \
  spamm_print_node.c \
  spamm_print_spamm.c \
  spamm_print_tree.c \
  spamm_read_MM.c \
  spamm_set.c \
  spamm_spamm_to_dense.c \
  spamm_tree_stats.c \
  spamm_version.c

matrix_kernel_sources = \
  spamm_kernel.h \
  spamm_stream_kernel.S \
  spamm_trivial_kernel.c

EXTRA_DIST = spamm_mainpage.h.in

tags :
	ctags --C-kinds=+l+x *.h *.c *.S ../config.h
