include $(top_srcdir)/Makefile.charm

charmtest_SOURCES = \
  charmtest.cc \
  charmtest.ci

charmtest$(EXEEXT) : charmtest.$(OBJEXT)
	$(CHARMC) -o $@ $^

charmtest$(OBJECT) : charmtest.cc testmodule.decl.h testmodule.def.h
	$(CHARMC) -o $@ -c $^

testmodule.decl.h testmodule.def.h : charmtest.ci
	$(CHARMC) $^

CLEANFILES = testmodule.decl.h testmodule.def.h
