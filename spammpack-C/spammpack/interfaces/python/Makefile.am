nodist_pkgpython_PYTHON = spammpack.py
pkgpyexec_LTLIBRARIES = _spammpack.la

spammpack.py : spammpack_wrap.c

_spammpack_la_SOURCES = spammpack.i spammpack_wrap.c
_spammpack_la_CPPFLAGS = -I$(PYTHON_INCLUDE) -I$(top_srcdir)/src -I$(top_builddir)/src
_spammpack_la_LDFLAGS = -module
_spammpack_la_LIBADD = $(top_builddir)/src/libspamm.la

%_wrap.c : $(srcdir)/%.i
	swig -python -I$(top_srcdir)/src -o $@ $<

CLEANFILES = *_wrap.c spammpack.py
