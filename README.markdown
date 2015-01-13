# Build Instructions

This library uses CMake. Simply run

~~~~~
cmake .
make
~~~~~

to build it.

# Building with full debugging

Configurating the library using

~~~
cmake -DCMAKE_BUILD_TYPE=Debug .
~~~

will set compiler flag defaults that include bounds checking (we currently
support GNU and Intel compilers).

# Building with a specific compiler

When `cmake` is run, it searches for a suitable Fortran compiler. If the one
it finds is not the one you would like to use, you have to explicitly specify
the compiler as

~~~
cmake -DCMAKE_Fortran_COMPILER=FCompiler
~~~

where the `FCompiler` can be the name of the compiler executable with or
without an absolute path.

# Further Information

For more information, see the doxygen documentation for details, or visit
http://freeon.org/spammpack.
