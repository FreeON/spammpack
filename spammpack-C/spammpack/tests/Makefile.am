AM_CTAGSFLAGS = --C-kinds=+l+x --Fortran-kinds=+i+L
TAGS_FILES = \
  $(top_builddir)/config.h \
  $(top_srcdir)/src/*.h \
  $(top_srcdir)/src/*.c \
  $(top_srcdir)/src/*.S

ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

check_LTLIBRARIES = libtest.la

libtest_la_SOURCES = \
	libtest.c \
  libtest_generate_matrix.c \
  libtest_compare_spamm_to_dense.c \
	test.h

noinst_HEADERS = \
  libtest_generate_matrix_sources.c \
  libtest_compare_spamm_to_dense_sources.c

AM_CFLAGS = $(OPENMP_CFLAGS)
AM_FCFLAGS = $(OPENMP_FCFLAGS)
AM_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/src -I$(top_builddir)/interfaces/Fortran
AM_LDFLAGS = $(top_builddir)/src/libspamm.la

if HAVE_SSE_IN_CPU
multiply_spamm_SSE = multiply_spamm_kernel_standard_SSE.sh
endif

if HAVE_SSE4_1_IN_CPU
multiply_spamm_SSE4_1 = multiply_spamm_kernel_standard_SSE4_1.sh
endif

multiply_tests = $(multiply_spamm_SSE) $(multiply_spamm_SSE4_1)

03_create_spamm_LDADD = libtest.la
04_copy_spamm_LDADD = libtest.la
06_chunk_test_SOURCES = 06_chunk_test.F90
06_chunk_test_LDADD = $(top_builddir)/interfaces/Fortran/libspammpack_fortran.la
08_get_version_SOURCES = 08_get_version.F90
08_get_version_LDADD = $(top_builddir)/interfaces/Fortran/libspammpack_fortran.la
14_add_spamm_LDADD = libtest.la
17_multiply_spamm_LDADD = libtest.la
20_seteq_SOURCES = 20_seteq.F90
20_seteq_LDADD = $(top_builddir)/interfaces/Fortran/libspammpack_fortran.la
21_trace_SOURCES = 21_trace.F90
21_trace_LDADD = $(top_builddir)/interfaces/Fortran/libspammpack_fortran.la

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
	14_add_spamm \
	16_stats \
	17_multiply_spamm \
	18_spectral_bounds \
	20_seteq \
	21_trace

TESTS = \
  $(check_PROGRAMS)

CLEANFILES = *.mod
