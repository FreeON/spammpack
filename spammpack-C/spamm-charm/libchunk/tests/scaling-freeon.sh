#!/bin/bash
#
# vim: tw=0

if [[ ${BUILD_TYPE:=serial} = "serial" ]]; then
  echo "serial version..."
  THREADS=( 1 )
  NUMA_POLICY="--physcpubind=0 --membind=0"
  CONFIGURE_ARGS="--disable-openmp"
  CONFIGURE_ARGS+=" --disable-assert"
  CONFIGURE_ARGS+=" --enable-chunk-block-transpose=0"
  CONFIGURE_ARGS+=" --enable-no-work"
elif [[ ${BUILD_TYPE} = "openmp" ]]; then
  echo "OpenMP version..."
  THREADS=( 1 2 4 8 12 16 20 24 28 32 36 40 44 48 )
  NUMA_POLICY="--interleave=all"
  CONFIGURE_ARGS="--enable-openmp"
  CONFIGURE_ARGS+=" --disable-assert"
  CONFIGURE_ARGS+=" --enable-chunk-block-transpose=0"
  CONFIGURE_ARGS+=" --enable-chunk-tree-max-tier=3"
  #CONFIGURE_ARGS+=" --enable-no-work"
else
  echo "unknown build type, either serial or openmp"
  exit 1
fi

pushd ..
if [[ ${BUILD_COMPILER:=gcc} = "gcc" ]]; then
  # Building with gcc:
  ./configure \
    ${CONFIGURE_ARGS} \
    CC=gcc \
    CFLAGS="-Ofast -g -ftree-vectorizer-verbose=2 -mfpmath=sse -msse2" || exit
elif [[ ${BUILD_COMPILER} = "intel" ]]; then
  ./configure \
    ${CONFIGURE_ARGS} \
    CC=icc \
    CFLAGS="-O3 -g -no-prec-div -static -fp-model fast=2 -qopt-report -xHost -ipo" || exit
else
  echo "unknown build type ${BUILD_COMPILER}"
  exit 1
fi

make clean || exit
make V=1 || exit
popd

head ../config.log | grep "$.*config" | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output || exit

echo "running ${BUILD_TYPE}-${BUILD_COMPILER}" | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output || exit
echo "numactl policy: ${NUMA_POLICY}" | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output || exit

for tolerance in 1e-10 1e-6; do
  for threads in ${THREADS[*]}; do
    OMP_NUM_THREADS=${threads} \
      numactl ${NUMA_POLICY} -- \
      ./chunk_multiply \
      --tolerance ${tolerance} \
      --N_basic 4 \
      --type BCSR \
      --bcsr ~/data-sets/water/h2o_90.OrthoD \
      --no-verify \
      --repeat 20 | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output \
      || exit
  done
done
