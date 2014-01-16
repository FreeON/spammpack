if HAVE_PYTHON
python_tests = test_sp2.sh
endif

noinst_PROGRAMS = \
	chunk_multiply \
	chunk_tree_multiply

check_PROGRAMS = \
	chunk_add \
	chunk_add_identity \
	chunk_scale \
	chunk_trace \
	chunk_tree_alloc

AM_CFLAGS = $(OPENMP_CFLAGS)
AM_CPPFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS = $(top_builddir)/src/libchunk.la

TESTS = \
	$(python_tests) \
	$(check_PROGRAMS) \
	$(noinst_PROGRAMS) \
	spamm_test_SP2_full.sh \
	spamm_test_add_full.sh \
	spamm_test_add_decay.sh \
	spamm_test_add_diagonal.sh \
	spamm_test_multiply_full.sh \
	spamm_test_multiply_decay.sh \
	spamm_test_multiply_diagonal.sh \
	spamm_test_scale_full.sh \
	spamm_test_scale_decay.sh \
	spamm_test_scale_diagonal.sh \
	spamm_test_trace_full.sh \
	spamm_test_trace_decay.sh \
	spamm_test_trace_diagonal.sh \
	spamm_test_addIdentity_full.sh \
	spamm_test_addIdentity_decay.sh \
	spamm_test_addIdentity_diagonal.sh

AM_TESTS_ENVIRONMENT = top_srcdir=${top_srcdir} top_builddir=${top_builddir}

DISTCLEANFILES = $(TESTS) random.mm

TAGS_FILES = $(top_srcdir)/src/*.c $(top_srcdir)/src/*.h
