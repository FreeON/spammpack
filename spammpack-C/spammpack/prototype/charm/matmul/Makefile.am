EXTRA_DIST = $(INTERFACE_FILES)

INTERFACE_FILES = \
	matmul.ci \
	matrix.ci \
	node.ci

AM_CXXFLAGS = @LOADBALANCER_FLAGS@

bin_PROGRAMS = matmul

matmul_SOURCES = \
	matmul.cc \
	matrix.cc \
	matrix.h \
	node.cc \
	node.h

%.decl.h %.def.h : %.ci
	$(AM_V_CXX)$(CHARMC) -preprocess $<

AM_LDFLAGS = -language charm++

BUILT_SOURCES = $(foreach i, $(INTERFACE_FILES), $(basename $i).decl.h $(basename $i).def.h)

CLEANFILES = $(BUILT_SOURCES) charmrun
