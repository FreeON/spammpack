ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

AM_CPPFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS = $(top_builddir)/src/libspamm.la

noinst_PROGRAMS = \
  create_spamm \
  check_spamm \
  convert_dense \
  hashtable_remove \
  hashtable_insert \
  list \
  list_sort \
  multiply_spamm

TESTS = $(noinst_PROGRAMS)
