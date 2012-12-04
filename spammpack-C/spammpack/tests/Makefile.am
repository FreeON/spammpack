AM_CTAGSFLAGS = --C-kinds=+l+x --Fortran-kinds=+i+L
TAGS_FILES = \
  $(top_srcdir)/src/*.h \
  $(top_srcdir)/src/*.c \
  $(top_srcdir)/src/*.S

ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

AM_CFLAGS = $(OPENMP_CFLAGS)
AM_FCFLAGS = $(OPENMP_FCFLAGS)
AM_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/src
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

spammtests.lo : test_utilities.lo

noinst_PROGRAMS = \
  convert_to_binary

08_get_version_SOURCES = 08_get_version.F90

spamm_multiply_SOURCES = spamm_multiply.F90
spamm_multiply_LDFLAGS = libspammtest.la
spamm_multiply_DEPENDENCIES = libspammtest.la

SP2_SOURCES = SP2.F90 spamm_project.F90
SP2_LDFLAGS = libspammtest.la
SP2_DEPENDENCIES = libspammtest.la

06_chunk_test_SOURCES = 06_chunk_test.F90
06_chunk_test_LDFLAGS = \
  $(top_builddir)/src/libspammpack_fortran.la \
  $(top_builddir)/src/libspamm.la

noinst_SCRIPTS = \
  create_spamm_row_major.sh \
  create_spamm_Z_curve.sh \
  convert_dense_row_major.sh \
  convert_dense_Z_curve.sh \
  multiply_spamm_kernel_standard_SSE.sh \
  multiply_spamm_kernel_standard_SSE4_1.sh

EXTRA_DIST = $(noinst_SCRIPTS)

check_PROGRAMS = \
  03_create_spamm \
  04_copy_spamm \
  05_new_chunk \
  06_chunk_test \
  07_convert_dense \
  08_get_version \
  11_ipow \
  12_index \
  12_linear_index \
  13_check_spamm \
  14_add_spamm_dense \
  15_add_spamm_sparse \
  16_stats \
  17_multiply_spamm \
  18_spectral_bounds

# spamm_multiply

TESTS = \
  $(check_PROGRAMS)

#$(multiply_tests)

CLEANFILES = *.mod
