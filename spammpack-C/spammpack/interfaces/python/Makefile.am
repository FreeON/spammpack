pkgpython_PYTHON = spammpack.py
pkgpyexec_LTLIBRARIES = _spammpack.la

_spammpack_la_SOURCES = $(srcdir)/spamm_version_wrap.c spamm_version.i
_spammpack_la_CPPFLAGS = -I@PYTHON_INCLUDE@ -I$(top_srcdir)/src
_spammpack_la_LDFLAGS = -module
_spammpack_la_LIBADD = $(top_builddir)/src/libspamm.la

$(srcdir)/spamm_version_wrap.c : spamm_version.i
	swig -python -I$(top_srcdir)/src -o $@ $<
