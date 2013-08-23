EXTRA_DIST = \
	$(INTERFACE_FILES) \
	configure.LB \
	configure.mustang \
	freeon-submit.sh \
	freeon.job.template \
	moonlight-submit.sh \
	moonlight.sh.template \
	mustang-results.sh \
	mustang-submit.sh \
	mustang.job.template

AM_LDFLAGS = -language charm++ @LOADBALANCER_FLAGS@ @DEFAULT_LB_FLAGS@ @CHARMDEBUG_FLAGS@

bin_PROGRAMS = matmul

INTERFACE_FILES = \
	matmul.ci \
	matrix.ci \
	messages.ci \
	multiply.ci \
	multiplyelement.ci \
	node.ci

matmul_SOURCES = \
	index.h \
	blas_interface.h \
	logger.cc \
	logger.h \
	matmul.cc \
	matmul.h \
	matrix.cc \
	matrix.h \
	messages.cc \
	messages.h \
	multiply.cc \
	multiply.h \
	multiplyelement.cc \
	multiplyelement.h \
	node.cc \
	node.h \
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

run-LB : matmul$(EXEEXT)
	./charmrun +p4 matmul -N 2048 -b 256 --verify +LBDebug

valgrind : matmul$(EXEEXT)
	valgrind --log-file=valgrind.%p.log \
		--track-origins=yes \
		--leak-check=full \
		--show-reachable=yes \
		./$< -N 1

docs :
	$(DOXYGEN)

%.decl.h %.def.h : %.ci
	$(AM_V_CXX)$(CHARMC) $(DEFAULT_INCLUDES) -preprocess $<

BUILT_SOURCES = $(foreach i, $(INTERFACE_FILES), $(basename $i).decl.h $(basename $i).def.h)

CLEANFILES = $(BUILT_SOURCES) charmrun
