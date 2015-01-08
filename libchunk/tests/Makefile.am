dist_noinst_SCRIPTS = scaling-mustang.sh

noinst_PROGRAMS = \
	chunk_multiply \
	chunk_tree_multiply \
	chunk_SP2

check_PROGRAMS = \
	chunk_add \
	chunk_add_identity \
	chunk_scale \
	chunk_trace \
	chunk_tree_alloc

chunk_multiply_SOURCES = \
	chunk_multiply.c \
	bcsr.h \
	bcsr.c

chunk_SP2_SOURCES = \
	chunk_SP2.c \
	bcsr.h \
	bcsr.c

AM_CFLAGS = $(OPENMP_CFLAGS) $(MIC_CFLAGS)
AM_CPPFLAGS = -I $(top_srcdir)/src
AM_LDFLAGS = ../src/libchunk.la $(MIC_LDFLAGS)

TESTS = \
	$(check_PROGRAMS) \
	$(noinst_PROGRAMS)
