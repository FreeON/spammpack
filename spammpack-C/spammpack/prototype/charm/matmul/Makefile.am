EXTRA_DIST = \
	$(INTERFACE_FILES) \
	freeon-submit.sh \
	freeon.job.template \
	moonlight-submit.sh \
	moonlight.sh.template \
	mustang-submit.sh \
	mustang-results.sh \
	mustang.job.template \
	configure.mustang \
	configure.LB

AM_LDFLAGS = -language charm++
AM_CXXFLAGS = @LOADBALANCER_FLAGS@

bin_PROGRAMS = matmul

INTERFACE_FILES = \
	matmul.ci \
	matrix.ci \
	messages.ci \
	multiply.ci \
	node.ci \
	reduction.ci

matmul_SOURCES = \
	index.h \
	logger.h \
	matmul.cc \
	matrix.cc \
	matrix.h \
	messages.cc \
	messages.h \
	multiply.cc \
	multiply.h \
	node.cc \
	node.h \
	reduction.cc \
	reduction.h \
	timer.cc \
	timer.h \
	types.h

$(matmul_OBJECTS) : $(BUILT_SOURCES)

N = 2048
B = 64

RSH = rsh

run : matmul$(EXEEXT)
	for i in 1 2 4 8 12 16 20 24 28 32 36 40 44 48; do \
		./charmrun +p$$i ++remote-shell $(RSH) matmul -N $(N) -b $(B); \
	done

run-local : matmul$(EXEEXT)
	for i in 1 2 4 8 12 16 20 24 28 32 36 40 44 48; do \
		./matmul +p$$i -N $(N) -b $(B); \
	done

valgrind : matmul$(EXEEXT)
	valgrind --log-file=valgrind.%p.log \
		--track-origins=yes \
		--leak-check=full \
		--show-reachable=yes \
		./$< -N 1

docs :
	$(DOXYGEN)

%.decl.h %.def.h : %.ci
	$(AM_V_CXX)$(CHARMC) -preprocess $<

BUILT_SOURCES = $(foreach i, $(INTERFACE_FILES), $(basename $i).decl.h $(basename $i).def.h)

CLEANFILES = $(BUILT_SOURCES) charmrun
