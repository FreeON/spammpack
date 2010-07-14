ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

lib_LTLIBRARIES = libspamm.la

libspamm_la_SOURCES = \
  $(matrix_tree_sources) \
  $(linked_list_sources) \
  $(memory_manager_sources) \
  $(matrix_kernel_sources)

matrix_tree_sources = \
  spamm.h \
  spamm_add.c \
  spamm_compare_int.c \
  spamm_compare_multiply_stream_element.c \
  spamm_delete.c \
  spamm_delete_node.c \
  spamm_dense_index.c \
  spamm_dense_to_spamm.c \
  spamm_get.c \
  spamm_int_to_binary.c \
  spamm_linear_to_coordinates.c \
  spamm_logger.c \
  spamm_mask.c \
  spamm_multiply.c \
  spamm_multiply_scalar.c \
  spamm_new.c \
  spamm_new_linear_quadtree_node.c \
  spamm_new_node.c \
  spamm_number_nonzero.c \
  spamm_print_dense.c \
  spamm_print_multiply_stream.c \
  spamm_print_node.c \
  spamm_print_spamm.c \
  spamm_print_tree.c \
  spamm_read_MM.c \
  spamm_set.c \
  spamm_spamm_to_dense.c \
  spamm_swap_block_dense.c \
  spamm_swap_floating_point_t.c \
  spamm_swap_floating_point_t_pointer.c \
  spamm_swap_linear_quadtree.c \
  spamm_swap_multiply_stream.c \
  spamm_swap_short_unsigned_int.c \
  spamm_swap_spamm_node_t_pointer.c \
  spamm_swap_unsigned_int.c \
  spamm_tree_pack.c \
  spamm_tree_stats.c \
  spamm_version.c

linked_list_sources = \
  spamm_ll.h \
  spamm_ll_append.c \
  spamm_ll_delete.c \
  spamm_ll_delete_node.c \
  spamm_ll_get.c \
  spamm_ll_iterator_delete.c \
  spamm_ll_iterator_first.c \
  spamm_ll_iterator_new.c \
  spamm_ll_iterator_next.c \
  spamm_ll_new.c \
  spamm_ll_new_node.c \
  spamm_ll_print.c \
  spamm_ll_print_node.c \
  spamm_ll_sort.c \
  spamm_ll_sort_data.c \
  spamm_ll_swap.c

memory_manager_sources = \
  spamm_mm.h \
  spamm_mm_allocate.c \
  spamm_mm_delete.c \
  spamm_mm_delete_chunk.c \
  spamm_mm_grow.c \
  spamm_mm_new.c \
  spamm_mm_new_chunk.c \
  spamm_mm_print.c

matrix_kernel_sources = \
  spamm_trivial_kernel.c

EXTRA_DIST = spamm_mainpage.h.in

tags :
	ctags *.h *.c
