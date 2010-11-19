ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

lib_LTLIBRARIES = libspamm.la

noinst_LTLIBRARIES = libspamm_double.la libspamm_single.la libspamm_kernel.la

libspamm_la_SOURCES =
libspamm_la_LIBADD = $(noinst_LTLIBRARIES)

libspamm_kernel_la_SOURCES = $(matrix_kernel_sources)

libspamm_double_la_CPPFLAGS = -DPRECISION=DOUBLE
libspamm_double_la_SOURCES = $(matrix_tree_sources)

libspamm_single_la_CPPFLAGS = -DPRECISION=SINGLE
libspamm_single_la_SOURCES = $(matrix_tree_sources)

matrix_tree_sources = \
  spamm.h

matrix_kernel_sources = \
  spamm_kernel.h \
  spamm_stream_kernel.S

tags :
	ctags --C-kinds=+l+x *.h *.c *.S ../config.h
