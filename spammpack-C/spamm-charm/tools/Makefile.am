if HAVE_PYTHON
python_PYTHON = plotPEMap.py sp2.py spectral_bounds.py matrix_market.py
TESTS = $(top_builddir)/sp2_test.sh
AM_TESTS_ENVIRONMENT = export top_builddir=${top_builddir};
check_SCRIPTS = $(TESTS)
endif

DISTCLEANFILES = $(TESTS) $(top_builddir)/random.mm
