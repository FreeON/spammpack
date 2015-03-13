---
title: Spammpack
---

# Spammpack

Spammpack is an implementation of the Sparse Approximate Matrix
Multiplication (SpAMM) algorithm introduced in Ref.
{% cite Challacombe2010 %}, and further developed in Refs.
{% cite Bock2013 Bock2014 %}.  It provides a matrix data type, and an
approximate matrix product, which exhibits linear scaling
computational complexity for matrices with decay. The product error
and the performance of the multiply can be tuned by choosing an
appropriate tolerance. The library can be compiled for serial
executation or parallel execution on shared memory systems with an
OpenMP capable compiler.  Currently under heavy development, the
latest version can be obtained through git:

    $ git clone https://github.com/FreeON/spammpack.git

Alternatively we offer a snapshot of the latest version in a zip or
tar archive through the links on the right hand side. See
[Building](/spammpack/building.html) for more details on how to build
the library.

# License

Spammpack is licensed under the
[BSD 3-clause license](http://opensource.org/licenses/BSD-3-Clause).

{% include responsive_ad.html %}

# Project Statistics

The latest project statistics including a source code line count and
some other factoids can be found on
[OpenHub](https://www.openhub.net/p/spammpack):

{% include openhub.html %}

# API Documentation

The [API](/spammpack/html/) is documented using
[Doxygen](http://www.doxygen.org).

# Authors

The principal authors of spammpack are:

  - Matt Challacombe @mattchallacombe
  - Nicolas Bock @nicolasbock

# Bibliography

{% bibliography --cited_in_order %}
