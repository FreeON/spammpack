if HAVE_PYTHON
python_tests = test_sp2.sh
endif

TESTS = \
	$(python_tests) \
	spamm_test_add_full.sh \
	spamm_test_add_decay.sh \
	spamm_test_add_diagonal.sh \
	spamm_test_multiply_full.sh \
	spamm_test_multiply_decay.sh \
	spamm_test_multiply_diagonal.sh \
	spamm_test_trace_full.sh \
	spamm_test_trace_decay.sh \
	spamm_test_trace_diagonal.sh

AM_TESTS_ENVIRONMENT = top_srcdir=${top_srcdir} top_builddir=${top_builddir}

DISTCLEANFILES = $(TESTS) random.mm
