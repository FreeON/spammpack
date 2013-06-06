if HAVE_CHARM

include $(top_srcdir)/Makefile.charm

bin_PROGRAMS = spammmultiply

BUILT_SOURCES = \
	spamm.decl.h \
	spamm.def.h \
	spammmultiply.def.h \
	spammmultiply.def.h

INTERFACE_FILES = \
	spamm.ci \
	spammmultiply.ci

spammmultiply_SOURCES = \
	$(INTERFACE_FILES) \
	spamm.cc \
	spammmultiply.cc

endif
