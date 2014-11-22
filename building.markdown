---
title: Spammpack - Building
---

# Building

Download a zip or tar archive using the links on the right. Unpack and change
into the top level source directory:

    $ tar xf FreeON-spammpack-v1.0.0-104-g24641da.tar.gz
    $ cd FreeON-spammpack-24641da/

The build system uses `CMake`. Configure the sources by running

    $ cmake . -DCMAKE_INSTALL_PREFIX=/some/path

And build them using `make`

    $ make

{% include responsive_ad.html %}

# Testing

Spammpack comes with a regression test suite. If you want to make sure (and
why wouldn't you?) that everything built correctly, run

    $ make test

In addition, we let
[Travis-CI](https://travis-ci.org)
run through the test suite.  The latest result of this build/test process is
shown on the right. Current condition:

| Branch | Build Status |
| ------ | ------------ |
| master | [![Build Status](https://travis-ci.org/FreeON/spammpack.svg?branch=master)](https://travis-ci.org/FreeON/spammpack) |
| v1.0   | [![Build Status](https://travis-ci.org/FreeON/spammpack.svg?branch=v1.0)](https://travis-ci.org/FreeON/spammpack) |

# Installing

The `install` target will install the library into
`${CMAKE_INSTALL_PREFIX}/lib`.

    $ make install
