EXTRA_DIST = $(INTERFACE_FILES)

INTERFACE_FILES = \
	matrix.ci \
	messages.ci \
	node.ci \
	matmul.ci

bin_PROGRAMS = matmul

matmul_SOURCES = \
	matrix.cc \
	matrix.h \
	messages.cc \
	messages.h \
	node.cc \
	node.h \
	matmul.cc

%.decl.h %.def.h : %.ci
	$(AM_V_CXX)$(CHARMC) -preprocess $<

AM_LDFLAGS = -language charm++

BUILT_SOURCES = $(foreach i, $(INTERFACE_FILES), $(basename $i).decl.h $(basename $i).def.h)

CLEANFILES = $(BUILT_SOURCES) charmrun
