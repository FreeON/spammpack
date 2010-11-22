ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

lib_LTLIBRARIES = libspamm.la

noinst_LTLIBRARIES = libspamm_tree.la libspamm_kernel.la

libspamm_la_SOURCES =
libspamm_la_LIBADD = $(noinst_LTLIBRARIES)

libspamm_kernel_la_SOURCES = \
  spamm_kernel.h \
  spamm_stream_kernel.S

libspamm_tree_la_SOURCES = \
  spamm.h \
  spamm_config.h \
  spamm_delete.c \
  spamm_delete_node.c \
  spamm_delete_block.c \
  spamm_get.c \
  spamm_hashtable.c \
  spamm_index_2D.c \
  spamm_multiply.c \
  spamm_new.c \
  spamm_new_block.c \
  spamm_new_node.c \
  spamm_print.c \
  spamm_set.c

tags :
	ctags --C-kinds=+l+x *.h *.c *.S ../config.h
