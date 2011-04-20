AM_CTAGSFLAGS = --C-kinds=+l+x

ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

AM_CPPFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS = $(top_builddir)/src/libspamm.la

noinst_PROGRAMS = \
  create_spamm \
  expand \
  stats \
  check_spamm \
  convert_dense \
  hashtable_remove \
  hashtable_insert \
  list \
  list_sort \
  multiply_spamm

TESTS = \
  create_spamm \
  expand \
  stats \
  check_spamm \
  convert_dense \
  hashtable_remove \
  hashtable_insert \
  list \
  list_sort \
  multiply_spamm_kernel_standard_SSE.sh \
  multiply_spamm_kernel_standard_SSE4_1.sh \
  multiply_spamm_kernel_Z_curve_SSE.sh \
  multiply_spamm_kernel_Z_curve_SSE4_1.sh \
  multiply_spamm_kernel_hierarchical_SSE.sh \
  multiply_spamm_kernel_hierarchical_SSE4_1.sh \
  multiply_spamm_kernel_hierarchical_Z_curve_SSE.sh \
  multiply_spamm_kernel_hierarchical_Z_curve_SSE4_1.sh
