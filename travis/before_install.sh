#!/bin/bash

add-apt-repository ppa:ubuntu-toolchain-r/test -y
apt-get update -qq
apt-get install -y gfortran-${GCC_VERSION}
