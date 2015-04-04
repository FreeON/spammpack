#!/bin/bash
#
# vim: tw=0

echo "building with ${BUILD_TYPE:=gcc}..."

pushd ..
if [[ ${BUILD_TYPE} = "gcc" ]]; then
  # Building with gcc:
  ./configure --disable-assert \
    --enable-chunk-tree-max-tier=4 \
    --disable-openmp \
    CC=gcc \
    CFLAGS="-Ofast" || exit
elif [[ ${BUILD_TYPE} = "intel" ]]; then
  ./configure --disable-assert \
    --enable-chunk-tree-max-tier=4 \
    --disable-openmp \
    CC=icc \
    CFLAGS="-O3 -no-prec-div -static -fp-model fast=2" || exit
else
  echo "unknown build type ${BUILD_TYPE}"
  exit 1
fi

make clean || exit
make || exit
popd

for tolerance in 1e-10 1e-6; do
  numactl \
    --membind=0 --pyscpubind=0 -- \
    ./chunk_multiply \
    --tolerance ${tolerance} \
    --N_basic 4 \
    --type BCSR \
    --bcsr ~/data-sets/water/h2o_90.OrthoD \
    --no-verify \
    --repeat 10 | tee --append serial-${BUILD_TYPE}.output \
    || exit
done
