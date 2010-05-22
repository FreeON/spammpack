SUBDIRS = sparsekit

ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

lib_LTLIBRARIES = libspamm.la

libspamm_la_SOURCES = \
  spamm.h \
  spamm_add.c \
  spamm_delete.c \
  spamm_dense_index.c \
  spamm_dense_to_spamm.c \
  spamm_get.c \
  spamm_linked_list.c \
  spamm_logger.c \
  spamm_multiply.c \
  spamm_new.c \
  spamm_number_nonzero.c \
  spamm_print_dense.c \
  spamm_print_node.c \
  spamm_print_spamm.c \
  spamm_print_tree.c \
  spamm_read_MM.c \
  spamm_set.c \
  spamm_spamm_to_dense.c \
  spamm_tree_stats.c

tags :
	ctags *.h *.c
