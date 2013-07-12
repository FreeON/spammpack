EXTRA_DIST = $(INTERFACE_FILES)

INTERFACE_FILES = \
	matrix.ci \
	node.ci \
	reductions.ci

bin_PROGRAMS = reductions

reductions_SOURCES = \
	matrix.cc \
	matrix.h \
	node.cc \
	node.h \
	reductions.cc

%.decl.h %.def.h : %.ci
	$(AM_V_CXX)$(CHARMC) -preprocess $<

AM_LDFLAGS = -language charm++

BUILT_SOURCES = $(foreach i, $(INTERFACE_FILES), $(basename $i).decl.h $(basename $i).def.h)

CLEANFILES = $(BUILT_SOURCES) charmrun
