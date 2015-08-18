SpAMM, the Sparse Approximate Matrix Multiply.

# Build Instructions

This library uses CMake. Simply run

~~~~~
cmake .
make
~~~~~

to build it. For extra convenience, we provide a shell script that
configures and builds spammpack out-of-source in a directory called
`build`.

~~~
./build.sh
~~~

will create a new directory `build` and compile the sources in that
directory. The script honors the following environment variables, and
passes them to the `cmake` command.

~~~
CMAKE_BUILD_TYPE
SPAMM_DEBUG_LEVEL
CMAKE_VERBOSE_MAKEFILE
CMAKE_C_COMPILER
CMAKE_Fortran_COMPILER
~~~

# Building with full debugging

Configurating the library using

~~~
cmake -DCMAKE_BUILD_TYPE=Debug -DSPAMM_DEBUG_LEVEL=2 .
~~~

will set compiler flag defaults that include bounds checking (we
currently support GNU and Intel compilers).

# Building with a specific compiler

When `cmake` is run, it searches for a suitable Fortran compiler. If
the one it finds is not the one you would like to use, you have to
explicitly specify the compiler as

~~~
cmake -DCMAKE_Fortran_COMPILER=FCompiler
~~~

where the `FCompiler` can be the name of the compiler executable with
or without an absolute path.

# Emacs suggested settings

For emacs users we have a suggested set of settings we use. Have a
look at
[emacs-settings.el](https://github.com/FreeON/spammpack/blob/master/emacs-settings.el). This
file can simply be copy and pasted into an existing `~/.emacs`
configuration file.

# Some other useful things to know

The `ctags` and `etags` targets generate tags files for `vim` and
`emacs`, respectively. The tags file is written to the repository root
directory.

~~~
make etags
~~~

# Further Information

For more information, see the doxygen documentation for details, or
visit http://freeon.org/spammpack.
