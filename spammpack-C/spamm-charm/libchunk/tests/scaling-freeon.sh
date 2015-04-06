#!/bin/bash
#
# vim: tw=0

if [[ ${BUILD_TYPE:=serial} = "serial" ]]; then
  echo "serial version..."
  THREADS=( 1 )
  NUMA_POLICY="--physcpubind=0 --membind=0"
  CONFIGURE="--disable-openmp"
  CONFIGURE+=" --enable-chunk-block-transpose=0"
  #CONFIGURE+=" --enable-no-work"
elif [[ ${BUILD_TYPE} = "openmp" ]]; then
  echo "OpenMP version..."
  THREADS=( 1 2 4 8 12 16 20 24 28 32 36 40 44 48 )
  NUMA_POLICY="--interleave=all"
  CONFIGURE="--enable-openmp"
  CONFIGURE+=" --enable-chunk-block-transpose=0"
  #CONFIGURE+=" --enable-no-work"
else
  echo "unknown build type, either serial or openmp"
  exit 1
fi

pushd ..
if [[ ${BUILD_COMPILER:=gcc} = "gcc" ]]; then
  # Building with gcc:
  ./configure --disable-assert \
    --enable-chunk-tree-max-tier=4 \
    ${CONFIGURE} \
    CC=gcc \
    CFLAGS="-Ofast -g -ftree-vectorizer-verbose=2 -mfpmath=sse -msse2" || exit
elif [[ ${BUILD_COMPILER} = "intel" ]]; then
  ./configure --disable-assert \
    --enable-chunk-tree-max-tier=4 \
    ${CONFIGURE} \
    CC=icc \
    CFLAGS="-O3 -g -no-prec-div -static -fp-model fast=2 -vec-report -xSSE2" || exit
else
  echo "unknown build type ${BUILD_COMPILER}"
  exit 1
fi

make clean || exit
make V=1 || exit
popd

echo "running ${BUILD_TYPE}-${BUILD_COMPILER}"
echo "numactl policy: ${NUMA_POLICY}"

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
      --repeat 10 | tee --append ${BUILD_TYPE}-${BUILD_COMPILER}.output \
      || exit
  done
done
