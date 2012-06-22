AM_CTAGSFLAGS = --C-kinds=+l+x
TAGS_FILES = $(top_srcdir)/src/*.h $(top_srcdir)/src/*.c $(top_srcdir)/src/*.S

ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/interfaces/Fortran
AM_LDFLAGS = $(top_builddir)/src/libspamm.la
#AM_LDFLAGS = $(top_builddir)/src/libspamm.la -all-static

if HAVE_SSE_IN_CPU
multiply_spamm_SSE = multiply_spamm_kernel_standard_SSE.sh
endif

if HAVE_SSE4_1_IN_CPU
multiply_spamm_SSE4_1 = multiply_spamm_kernel_standard_SSE4_1.sh
endif

multiply_tests = $(multiply_spamm_SSE) $(multiply_spamm_SSE4_1)

check_LTLIBRARIES = libspammtest.la
libspammtest_la_SOURCES = \
  test_utilities.F90 \
  spammtests.F90

check_PROGRAMS = \
  check_spamm \
  convert_dense \
  create_spamm \
  expand \
  hashtable_insert \
  hashtable_remove \
  list \
  stats \
  list_sort \
  multiply_recursive \
  multiply_spamm \
  get_version \
  SP2

noinst_PROGRAMS = \
  convert_to_binary

get_version_SOURCES = get_version.F90

noinst_SCRIPTS = \
  create_spamm_row_major.sh \
  create_spamm_Z_curve.sh \
  convert_dense_row_major.sh \
  convert_dense_Z_curve.sh \
  multiply_spamm_kernel_standard_SSE.sh \
  multiply_spamm_kernel_standard_SSE4_1.sh

EXTRA_DIST = $(noinst_SCRIPTS)

TESTS = \
  hashtable_remove \
  hashtable_insert \
  list \
  list_sort \
  expand \
  stats \
  check_spamm \
	$(multiply_tests) \
  multiply_recursive
