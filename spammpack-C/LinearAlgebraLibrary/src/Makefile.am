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
  spamm_tree_pack.c \
  spamm_tree_stats.c \
	$(linked_list_sources) \
	$(memory_manager_sources)

linked_list_sources = \
  spamm_ll.h \
  spamm_ll_append.c \
  spamm_ll_delete.c \
  spamm_ll_get.c \
  spamm_ll_initialize.c \
  spamm_ll_initialize_node.c \
  spamm_ll_print.c \
  spamm_ll_print_node.c \
  spamm_ll_sort.c \
  spamm_ll_swap.c

memory_manager_sources = \
  spamm_mm.h \
  spamm_mm_allocate.c \
  spamm_mm_grow.c \
  spamm_mm_initialize.c \
  spamm_mm_print.c

EXTRA_DIST = spamm_mainpage.h.in

tags :
	ctags *.h *.c
