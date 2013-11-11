if HAVE_PYTHON
python_tests = test_sp2.sh
endif

check_SCRIPTS = \
	$(python_tests) \
	spamm_test_add.sh \
	spamm_test_multiply.sh \
	spamm_test_trace.sh

TESTS = $(check_SCRIPTS)

AM_TESTS_ENVIRONMENT = top_srcdir=${top_srcdir} top_builddir=${top_builddir}

DISTCLEANFILES = $(TESTS) random.mm
