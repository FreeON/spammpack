EXTRA_DIST = \
	$(INTERFACE_FILES) \
	freeon-submit.sh \
	freeon.job.template \
	moonlight.sh

INTERFACE_FILES = \
	matmul.ci \
	matrix.ci \
	messages.ci \
	node.ci

AM_CXXFLAGS = @LOADBALANCER_FLAGS@

bin_PROGRAMS = matmul

matmul_SOURCES = \
	index.h \
	logger.h \
	matmul.cc \
	matrix.cc \
	matrix.h \
	messages.cc \
	messages.h \
	node.cc \
	node.h \
	timer.cc \
	timer.h \
	types.h

N = 2048
B = 64

RSH = rsh

run : matmul
	for i in 1 2 4 8 12 16 20 24 28 32 36 40 44 48; do \
		./charmrun +p$$i ++remote-shell $(RSH) matmul -N $(N) -b $(B); \
	done

run-local : matmul
	for i in 1 2 4 8 12 16 20 24 28 32 36 40 44 48; do \
		./matmul +p$$i -N $(N) -b $(B); \
	done

%.decl.h %.def.h : %.ci
	$(AM_V_CXX)$(CHARMC) -preprocess $<

AM_LDFLAGS = -language charm++

BUILT_SOURCES = $(foreach i, $(INTERFACE_FILES), $(basename $i).decl.h $(basename $i).def.h)

CLEANFILES = $(BUILT_SOURCES) charmrun
