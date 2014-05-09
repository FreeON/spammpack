TAGS_FILES = $(top_builddir)/config.h

AM_CFLAGS = $(OPENMP_CFLAGS) $(MIC_CFLAGS)
AM_LDFLAGS = $(MIC_LDFLAGS)

lib_LTLIBRARIES = libchunk.la

libchunk_la_SOURCES = \
	chunk.c \
	chunk_block.c \
	chunk_tiled.c \
	chunk_tree.c \
	chunk_mic.c \
	chunk_mic.h \
	lapack_interface.h \
	matrix_4x4.S

include_HEADERS = \
	chunk.h \
	chunk_block.h \
	chunk_tiled.h \
	chunk_tree.h
