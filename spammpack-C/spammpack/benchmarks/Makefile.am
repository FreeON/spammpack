ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

AM_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/sparsekit
AM_LDFLAGS = $(top_builddir)/src/libspamm.la $(top_builddir)/sparsekit/libsparsekit.la

EXTRA_DIST = \
  e20r0000.mtx \
  sparseHamiltonian.mtx \
  generateMatrix.py \
  README

#noinst_PROGRAMS = \
#  multiply_MM \
#  multiply \
#  SP2_MM \
#  latency

tags :
	ctags ../src/*.{c,h}
