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
echo "N_BASIC = ${N_BASIC:=4}"
echo "N_CHUNK = ${N_CHUNK:=128}"
echo "BCSR = ${BCSR:=~/data-sets/water/h2o_90.OrthoD}"
echo

if [[ ${BUILD_TYPE} = "serial" ]]; then
    echo "serial version..."
    THREADS=( 1 )
    NUMA_POLICY="--physcpubind=0 --membind=0"
    CONFIGURE_ARGS="--disable-openmp"
    CONFIGURE_ARGS+=" --disable-shared"
    CONFIGURE_ARGS+=" --disable-assert"
    CONFIGURE_ARGS+=" --enable-chunk-block-transpose=0"
    if [[ ${NO_WORK} = "TRUE" ]]; then
        CONFIGURE_ARGS+=" --enable-no-work"
    elif [[ ! ${NO_WORK} = "FALSE" ]]; then
        echo "unknown value for NO_WORK (${NO_WORK})"
        exit 1
    fi
elif [[ ${BUILD_TYPE} = "openmp" ]]; then
    echo "OpenMP version..."
    THREADS=( 1 2 4 8 12 16 20 24 28 32 36 40 44 48 )
    #THREADS=( 1 2 )
    NUMA_POLICY="--interleave=all"
    CONFIGURE_ARGS="--enable-openmp"
    CONFIGURE_ARGS+=" --disable-shared"
    CONFIGURE_ARGS+=" --disable-assert"
    CONFIGURE_ARGS+=" --enable-chunk-block-transpose=0"
    CONFIGURE_ARGS+=" --enable-chunk-tree-max-tier=${MAX_TIER}"
    if [[ ${NO_WORK} = "TRUE" ]]; then
        CONFIGURE_ARGS+=" --enable-no-work"
    elif [[ ! ${NO_WORK} = "FALSE" ]]; then
        echo "unknown value for NO_WORK (${NO_WORK})"
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
        CFLAGS="-O3 -g -no-prec-div -static -fp-model fast=2 -qopt-report -xSSE2" || exit
else
    echo "unknown build compiler ${BUILD_COMPILER}"
    exit 1
fi

make clean || exit
make V=1 || exit
popd

head ../config.log \
    | grep "$.*config" \
    | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output || exit

echo "running ${BUILD_TYPE}-${BUILD_COMPILER}" \
    | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output || exit
echo "numactl policy: ${NUMA_POLICY}" \
    | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output || exit

for tolerance in 1e-10 1e-6; do
    for threads in ${THREADS[*]}; do
      for (( repeat = 0; repeat < ${FULL_REPEAT}; repeat++ )); do
        OMP_NUM_THREADS=${threads} \
          numactl ${NUMA_POLICY} -- \
          ./chunk_multiply \
          --tolerance ${tolerance} \
          --N_basic ${N_BASIC} \
          --N_chunk ${N_CHUNK} \
          --type BCSR \
          --bcsr ${BCSR} \
          --no-verify \
          --repeat ${REPEAT} \
          | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output \
          || exit
      done
    done
done
