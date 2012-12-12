AM_CTAGSFLAGS = --C-kinds=+l+x
TAGS_FILES = $(top_builddir)/config.h *.py $(pkginclude_HEADERS)

AM_CFLAGS = $(OPENMP_CFLAGS)
AM_CPPFLAGS = $(OPENMP_CFLAGS)
AM_FCFLAGS = $(OPENMP_FCFLAGS)

ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

lib_LTLIBRARIES = libspamm.la libspammpack_fortran.la
noinst_LTLIBRARIES = libspamm_tree.la libspamm_kernel.la

libspamm_la_LDFLAGS = -version-info 0:0:0
libspamm_la_SOURCES =
libspamm_la_LIBADD  = libspamm_tree.la libspamm_kernel.la

noinst_SCRIPTS = SSERegister.py generate_SSE_assembly_for_chunks.py # spammOffsets.py
noinst_PROGRAMS = print_chunk_offsets

print_chunk_offsets_LDFLAGS = libspamm.la
print_chunk_offsets_DEPENDENCIES = libspamm.la

EXTRA_DIST = SSERegister.py generate_SSE_assembly_for_chunks.py

spamm_stream_kernel.S : $(srcdir)/generate_SSE_assembly_for_chunks.py
	PYTHONPATH=$(builddir) $(PYTHON) $(srcdir)/generate_SSE_assembly_for_chunks.py > $@

generated_kernel_files = spamm_stream_kernel.S

BUILT_SOURCES = $(generated_kernel_files) spammpack.mod

spammpack.mod : libspammpack_fortran.la

CLEANFILES = \
  $(generated_kernel_files) \
  *.mod

DISTCLEANFILES = \
  $(builddir)/spammOffsets.py \
  $(builddir)/spammOffsets.pyc \
  $(builddir)/__pycache__/*

nodist_libspamm_kernel_la_SOURCES = $(generated_kernel_files)

libspammpack_fortran_la_SOURCES = \
  spamm_fortran_interface.c \
  spammpack_algebra.F90 \
  spammpack_types.F90 \
  spammpack.F90

# Fortran dependencies.
spammpack.lo : spammpack_types.lo spammpack_algebra.lo
spammpack_algebra.lo : spammpack_types.lo

libspamm_kernel_la_SOURCES = \
  $(libspamm_kernel_headers) \
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
  spamm_index.c \
  spamm_ipow.c \
  spamm_multiply.c \
  spamm_new.c \
  spamm_print.c \
  spamm_set.c \
  spamm_sgemm.c \
  spamm_sort.c \
  spamm_spectral_bounds.c \
  spamm_stats.c \
  spamm_timer.c \
  spamm_uint_to_bin_string.c \
  spamm_version.c

pkginclude_HEADERS = \
  spamm_sort_source.c

libspamm_tree_headers = \
  spamm.h \
  spamm_chunk.h \
  spamm_error.h \
  spamm_errno.h \
  spamm_general.h \
  spamm_sort.h \
  spamm_timer.h \
  spamm_types.h

libspamm_tree_private_headers = \
  spamm_types_private.h

include_HEADERS = \
  $(libspamm_kernel_headers) \
  $(libspamm_tree_headers) \
  spammpack.mod
