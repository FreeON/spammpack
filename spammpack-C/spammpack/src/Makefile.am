ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

lib_LTLIBRARIES = libspamm.la

libspamm_la_SOURCES = \
  $(matrix_tree_sources) \
  $(matrix_kernel_sources)

matrix_tree_sources = \
  spamm.h

matrix_kernel_sources = \
  spamm_kernel.h \
  spamm_stream_kernel.S

tags :
	ctags --C-kinds=+l+x *.h *.c *.S ../config.h
