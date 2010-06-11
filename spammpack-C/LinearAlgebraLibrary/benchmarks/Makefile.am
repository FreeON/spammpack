ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

AM_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/sparsekit
AM_LDFLAGS = $(top_builddir)/src/libspamm.la $(top_builddir)/src/sparsekit/libsparsekit.la

EXTRA_DIST = \
  e20r0000.mtx \
  sparseHamiltonian.mtx \
  generateMatrix.py \
  README

noinst_PROGRAMS = multiply_MM SP2_MM latency

multiply_MM_SOURCES = multiply_MM.c
SP2_MM_SOURCES = SP2_MM.c
latency_SOURCES = latency.c

tags :
	ctags ../src/*.{c,h}
