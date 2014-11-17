---
layout: default
title: Spammpack - Building
---

## Building

Download a zip or tar archive using the links on the right. Unpack and change
into the top level source directory:

    $ tar xf FreeON-spammpack-v1.0.0-104-g24641da.tar.gz
    $ cd FreeON-spammpack-24641da/

The build system uses `CMake`. Configure the sources by running

    $ cmake . -DCMAKE_INSTALL_PREFIX=/some/path

And build them using `make`

    $ make

{% include responsive_ad.html %}

## Installing

The `install` target will install the library into
`${CMAKE_INSTALL_PREFIX}/lib`.

    $ make install
