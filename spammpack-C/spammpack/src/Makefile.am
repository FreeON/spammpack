AM_CTAGSFLAGS = --C-kinds=+l+x
TAGS_FILES = $(top_builddir)/config.h generate_SSE_assembly.py

AM_CFLAGS = $(OPENMP_CFLAGS)
AM_CPPFLAGS = $(OPENMP_CFLAGS)
AM_FCFLAGS = $(OPENMP_FCFLAGS)

ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

lib_LTLIBRARIES = libspamm.la libspammpack_fortran.la
noinst_LTLIBRARIES = libspamm_tree.la libspamm_kernel.la

if INTERNAL_LAPACK
LIBADD_INTERNAL_LAPACK = \
  $(top_builddir)/../lapack/lapack/libfreeonlapack.la \
  $(top_builddir)/../lapack/blas/libfreeonblas.la \
  $(top_builddir)/../lapack/install/libfreeonlapackextra.la
endif

libspamm_la_LDFLAGS = -version-info 0:0:0
libspamm_la_SOURCES =
libspamm_la_LIBADD  = libspamm_tree.la libspamm_kernel.la $(LIBADD_INTERNAL_LAPACK)
#libspamm_la_DEPENDENCIES = $(LIBADD_INTERNAL_LAPACK)

noinst_SCRIPTS = generate_SSE_assembly.py spammOffsets.py
noinst_PROGRAMS = print_data_sizes print_chunk_offsets

print_chunk_offsets_LDFLAGS = libspamm.la
print_chunk_offsets_DEPENDENCIES = libspamm.la

EXTRA_DIST = spamm_list_driver.c generate_SSE_assembly.py

spammOffsets.py : print_data_sizes$(EXEEXT)
	$(builddir)/print_data_sizes$(EXEEXT) > $@

spamm_stream_kernel_SSE.S : $(srcdir)/generate_SSE_assembly.py spammOffsets.py
	PYTHONPATH=$(builddir) $(PYTHON) $(srcdir)/generate_SSE_assembly.py -N 4 --SSE 1 --name spamm_stream_kernel_SSE > $@

spamm_stream_kernel_SSE4_1.S : $(srcdir)/generate_SSE_assembly.py spammOffsets.py
	PYTHONPATH=$(builddir) $(PYTHON) $(srcdir)/generate_SSE_assembly.py -N 4 --SSE 4.1 --name spamm_stream_kernel_SSE4_1 > $@

generated_kernel_files = \
  spamm_stream_kernel_SSE.S \
  spamm_stream_kernel_SSE4_1.S

BUILT_SOURCES = $(generated_kernel_files) spammpack.mod

spammpack.mod : libspammpack_fortran.la
spammpack.lo : spamm_fortran_derived.lo

CLEANFILES = \
  $(generated_kernel_files) \
  *.mod

DISTCLEANFILES = \
  $(builddir)/spammOffsets.py \
  $(builddir)/spammOffsets.pyc \
  $(builddir)/__pycache__/*

nodist_libspamm_kernel_la_SOURCES = $(generated_kernel_files)

libspammpack_fortran_la_SOURCES = \
  spamm_fortran_chunk.F90 \
  spamm_fortran_derived.F90 \
  spamm_fortran_interface.c \
  spamm_fortran_interface.h \
  spammpack.F90

libspamm_kernel_la_SOURCES = \
  $(libspamm_kernel_headers) \
  spamm_stream_external_sgemm.c \
  spamm_kernel.c

libspamm_kernel_headers = \
  spamm_kernel.h

libspamm_tree_la_SOURCES = \
  $(libspamm_tree_headers) \
  $(libspamm_tree_private_headers) \
  spamm_manual_mainpage.h \
  spamm_add.c \
  spamm_allocate.c \
  spamm_check.c \
  spamm_chunk.c \
  spamm_convert.c \
  spamm_copy.c \
  spamm_delete.c \
  spamm_error.c \
  spamm_get.c \
  spamm_hashtable.c \
  spamm_index.c \
  spamm_ipow.c \
  spamm_list.c \
  spamm_maintenance.c \
  spamm_multiply.c \
  spamm_new.c \
  spamm_print.c \
  spamm_set.c \
  spamm_sgemm.c \
  spamm_stats.c \
  spamm_timer.c \
  spamm_uint_to_bin_string.c \
  spamm_version.c

libspamm_tree_headers = \
  spamm.h \
  spamm_chunk.h \
  spamm_error.h \
  spamm_errno.h \
  spamm_general.h \
  spamm_hashtable.h \
  spamm_list.h \
  spamm_timer.h \
  spamm_types.h

libspamm_tree_private_headers = \
  spamm_types_private.h

include_HEADERS = \
  $(libspamm_kernel_headers) \
  $(libspamm_tree_headers) \
  spammpack.mod
