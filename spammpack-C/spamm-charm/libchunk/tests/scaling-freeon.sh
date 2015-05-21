#!/bin/bash
#
# vim: tw=0

echo "The following environment variables can be used to"
echo "control this script:"
echo
echo "BUILD_TYPE = ${BUILD_TYPE:=serial} {serial,openmp}"
echo "BUILD_COMPILER = ${BUILD_COMPILER:=gcc} {gcc,intel}"
echo "MAX_TIER = ${MAX_TIER:=5}"
echo "REPEAT = ${REPEAT:=4}"
echo "FULL_REPEAT = ${FULL_REPEAT:=1}"
echo "NO_WORK = ${NO_WORK:=FALSE} {TRUE,FALSE}"
echo "MATRIX_TYPE = ${MATRIX_TYPE:=BCSR} {BCSR,full}"
echo "EXTRA_WORK = ${EXTRA_WORK:=0}"
echo "N_BASIC = ${N_BASIC:=4}"
echo "N_CHUNK = ${N_CHUNK:=128}"
echo "N = ${N:=-1}"
echo "LOCK = ${LOCK:=TRUE}"
echo "BCSR = ${BCSR:=~/data-sets/water/h2o_90.OrthoD}"
echo

if [[ ${BUILD_TYPE} = "serial" ]]; then
    echo "serial version..."
    THREADS=( 1 )
    NUMA_POLICY="--physcpubind=0 --membind=0"
    CONFIGURE_ARGS="--disable-openmp"
    CONFIGURE_ARGS+=" --disable-shared"
    CONFIGURE_ARGS+=" --disable-assert"
    CONFIGURE_ARGS+=" --enable-extra-work=${EXTRA_WORK}"
    CONFIGURE_ARGS+=" --enable-chunk-block-transpose=0"
    if [[ ${NO_WORK} = "TRUE" ]]; then
        CONFIGURE_ARGS+=" --enable-no-work"
    elif [[ ! ${NO_WORK} = "FALSE" ]]; then
        echo "unknown value for NO_WORK (${NO_WORK})"
        exit 1
    fi
elif [[ ${BUILD_TYPE} = "openmp" ]]; then
    echo "OpenMP version..."
    THREADS=( 1 2 4 8 12 16 )
    #THREADS=( 1 2 4 8 12 16 20 24 28 32 36 40 44 48 )
    #THREADS=( 1 2 )
    NUMA_POLICY="--interleave=all"
    CONFIGURE_ARGS="--enable-openmp"
    CONFIGURE_ARGS+=" --disable-shared"
    CONFIGURE_ARGS+=" --disable-assert"
    CONFIGURE_ARGS+=" --enable-extra-work=${EXTRA_WORK}"
    CONFIGURE_ARGS+=" --enable-chunk-block-transpose=0"
    CONFIGURE_ARGS+=" --enable-chunk-tree-max-tier=${MAX_TIER}"
    if [[ ${NO_WORK} = "TRUE" ]]; then
        CONFIGURE_ARGS+=" --enable-no-work"
    elif [[ ! ${NO_WORK} = "FALSE" ]]; then
        echo "unknown value for NO_WORK (${NO_WORK})"
        exit 1
    fi
    if [[ ${LOCK} = "FALSE" ]]; then
        CONFIGURE_ARGS+=" --disable-lock"
    elif [[ ! ${LOCK} = "TRUE" ]]; then
        echo "unknown value for LOCK (${LOCK})"
        exit 1
    fi
else
    echo "unknown build type, either serial or openmp"
    exit 1
fi

pushd ..
if [[ ${BUILD_COMPILER} = "gcc" ]]; then
    # Building with gcc:
    ./configure \
        ${CONFIGURE_ARGS} \
        CC=gcc \
        CFLAGS="-Ofast -g -ftree-vectorizer-verbose=2 -mfpmath=sse -msse2" || exit
elif [[ ${BUILD_COMPILER} = "intel" ]]; then
    ./configure \
        ${CONFIGURE_ARGS} \
        CC=icc \
        CFLAGS="-O3 -g -no-prec-div -static -fp-model fast=2 -qopt-report -xSSE2 -ipo" || exit
else
    echo "unknown build compiler ${BUILD_COMPILER}"
    exit 1
fi

make clean || exit
make -j V=1 || exit
popd

head ../config.log \
    | grep "$.*config" \
    | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output || exit

echo "running ${BUILD_TYPE}-${BUILD_COMPILER}" \
    | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output || exit
echo "numactl policy: ${NUMA_POLICY}" \
    | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output || exit

if [[ ${N} -gt 0 ]]; then
  N_ARG="--N ${N}"
else
  N_ARG=""
fi

for tolerance in 1e-10 1e-6; do
  COMMAND_LINE="--tolerance ${tolerance}"
  COMMAND_LINE+=" --N_basic ${N_BASIC}"
  COMMAND_LINE+=" --N_chunk ${N_CHUNK}"
  COMMAND_LINE+=" ${N_ARG}"
  COMMAND_LINE+=" --type ${MATRIX_TYPE}"
  COMMAND_LINE+=" --bcsr ${BCSR}"
  COMMAND_LINE+=" --no-verify"
  COMMAND_LINE+=" --repeat ${REPEAT}"
  echo "running ./chunk_multiply ${COMMAND_LINE}"
  for threads in ${THREADS[*]}; do
    for (( repeat = 0; repeat < ${FULL_REPEAT}; repeat++ )); do
      OMP_NUM_THREADS=${threads} \
        numactl ${NUMA_POLICY} -- \
        ./chunk_multiply ${COMMAND_LINE} \
        | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output \
        || exit
    done \
      | grep "seconds" \
      | awk "{ times[N++] = \$15; };
             END {
               avg_time = 0;
               for(i = 0; i < N; i++) {
                 avg_time = avg_time+times[i];
               }
               avg_time = avg_time/N;
               var = 0;
               for(i = 0; i < N; i++) {
                 var = var+(times[i]-avg_time)**2;
               }
               var = var/N;
               printf(\"%d %e %e\n\", ${threads}, avg_time, sqrt(var));
             }"
  done | awk 'BEGIN { printf("%3s %12s %12s %5s %8s\n", "P", "walltime", "std.", "s/u", "eff."); };
              { if(!t0) { t0 = $2 };
                printf("% 3d %e %e %5.2f %7.3f%%\n", $1, $2, $3, t0/$2, 100*t0/$2/$1);
              }'
done
