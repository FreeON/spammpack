AM_CTAGSFLAGS = --C-kinds=+l+x
TAGS_FILES = $(top_builddir)/config.h *.py $(pkginclude_HEADERS)

AM_CFLAGS = $(OPENMP_CFLAGS)
AM_CPPFLAGS = $(OPENMP_CFLAGS)
AM_FCFLAGS = $(OPENMP_FCFLAGS)

ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

lib_LTLIBRARIES = libspamm.la

libspamm_la_LDFLAGS = -version-info 0:1:0

noinst_SCRIPTS = SSERegister.py generate_SSE_assembly_for_chunks.py

EXTRA_DIST = SSERegister.py generate_SSE_assembly_for_chunks.py

spamm_stream_kernel.S : $(srcdir)/generate_SSE_assembly_for_chunks.py
	PYTHONPATH=$(builddir) $(PYTHON) $(srcdir)/generate_SSE_assembly_for_chunks.py > $@

generated_kernel_files = spamm_stream_kernel.S

BUILT_SOURCES = $(generated_kernel_files)

CLEANFILES = $(generated_kernel_files)

DISTCLEANFILES = \
  $(builddir)/spammOffsets.py \
  $(builddir)/spammOffsets.pyc \
  $(builddir)/__pycache__/*

libspamm_la_SOURCES = \
  $(libspamm_kernel_headers) \
  $(libspamm_tree_headers) \
  $(generated_kernel_files) \
  $(libspamm_tree_private_headers) \
  spamm_kernel.c \
  spamm_add.c \
  spamm_allocate.c \
  spamm_blas.c \
  spamm_check.c \
  spamm_chunk.c \
  spamm_convert.c \
  spamm_copy.c \
  spamm_delete.c \
  spamm_error.c \
  spamm_get.c \
  spamm_index.c \
  spamm_introspection.c \
  spamm_ipow.c \
  spamm_multiply.c \
  spamm_new.c \
  spamm_omp_init.c \
  spamm_print.c \
  spamm_prune.c \
  spamm_set.c \
  spamm_sort.c \
  spamm_spectral_bounds.c \
  spamm_stats.c \
  spamm_timer.c \
  spamm_trace.c \
  spamm_uint_to_bin_string.c \
  spamm_version.c

pkginclude_HEADERS = \
  spamm_sort_source.c \
  spamm_blas_source.c

libspamm_kernel_headers = spamm_kernel.h

libspamm_tree_headers = \
  spamm.h \
  spamm_chunk.h \
  spamm_errno.h \
  spamm_error.h \
  spamm_general.h \
  spamm_introspection.h \
  spamm_prune.h \
  spamm_sort.h \
  spamm_timer.h \
  spamm_types.h

libspamm_tree_private_headers = spamm_types_private.h

include_HEADERS = \
  $(libspamm_kernel_headers) \
  $(libspamm_tree_headers)
