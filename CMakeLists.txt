# Copyright (c) 2015, Los Alamos National Laboratory
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 2.8)
project(spammpack C Fortran)

set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_PATCH 0)
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")

if(NOT CMAKE_C_COMPILER)
  message(FATAL_ERROR "Can not find a C compiler")
endif()

if(NOT CMAKE_Fortran_COMPILER)
  message(FATAL_ERROR "Can not find Fortran compiler")
endif()

# Default build type, override by explicitly running cmake with
# -DCMAKE_BUILD_TYPE={Debug,Release}.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(BUILD_SHARED_LIBS FALSE
  CACHE BOOL "Whether to build a shared library.")

# Make sure we get the right RPATH. We want the RPATH to be the build
# path and after installing, the actual installed library path.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

if(NOT CMAKE_Fortran_FLAGS)
  # Set some default compiler flags.
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(TEMP_FLAGS
      "-O0 -g -gdwarf-3"
      "-fvar-tracking-assignments -fcheck=all"
      "-Wall -Wimplicit-interface -Wuninitialized -Wimplicit-procedure")
    try_compile(HAVE_ADDRESS_SANITIZER
      ${CMAKE_BINARY_DIR}/CMakeTmp
      ${CMAKE_SOURCE_DIR}/cmake-tests/address-sanitizer
      address-sanitizer
      OUTPUT_VARIABLE ADDRESS_SANITIZER_OUTPUT)
    if(HAVE_ADDRESS_SANITIZER)
      message(STATUS "Adding address sanitizer compiler flag")
      list(APPEND TEMP_FLAGS "-fsanitize=address")
    else()
      message(STATUS "Address sanitizer not supported")
    endif()
    string(REGEX REPLACE ";" " " CMAKE_Fortran_FLAGS_DEBUG "${TEMP_FLAGS}")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -g")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS_DEBUG
      "-stand f08 -O0 -g"
      "-debug all -check all"
      "-warn unused -traceback")
    string(REGEX REPLACE ";" " " CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -no-prec-div -fp-model fast=2 -vec-report -g")
  else()
    message(STATUS "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
  endif()
endif()

unset(COMPILE_DEFINITIONS)

try_compile(COMPILE_RESULT_1
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/cmake-tests/error_stop.F90
  OUTPUT_VARIABLE COMPILE_OUTPUT_1)
if(NOT COMPILE_RESULT_1)
  message(FATAL_ERROR "The Fortran compiler does not understand the "
    "'error stop' statement")
endif()

try_run(RUN_RESULT_2 COMPILE_RESULT_2
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/cmake-tests/deferred_string_length.F90
  COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT_2)
if(NOT COMPILE_RESULT_2)
  message(FATAL_ERROR "The Fortran compiler can not compile deferred "
    "length string expressions")
endif()

if(NOT RUN_RESULT_2 EQUAL 0)
  message(WARNING
    "The Fortran compiler has trouble with deferred length string "
    "assignments, which might render screen output rather ugly. Please "
    "consider using a more recent compiler that fully implements "
    "this part of the the Fortran 2003 standard (e.g. >=gcc-4.8).")
else()
  set(HAVE_DEFERRED_STRING_LENGTH TRUE)
endif()

try_compile(COMPILE_RESULT_3
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/cmake-tests/long_lines.F90
  OUTPUT_VARIABLE COMPILE_OUTPUT_3)
if(COMPILE_RESULT_3)
  message(STATUS "The Fortran compiler supports lines > 132 characters")
  set(LOG_FILENAMES TRUE)
endif()

try_compile(COMPILE_TEST_4
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/cmake-tests/finalize.F90
  OUTPUT_VARIABLE COMPILE_OUTPUT_4)
if(COMPILE_TEST_4)
  set(HAVE_FINALIZE TRUE)
else()
  message(STATUS "The Fortran compiler does not support final type bound procedures")
endif()

try_compile(COMPILE_TEST_5
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/cmake-tests/constructor.F90
  OUTPUT_VARIABLE COMPILE_OUTPUT_5)
if(COMPILE_TEST_5)
  set(HAVE_CONSTRUCTOR TRUE)
else()
  message(STATUS "The Fortran compiler does not support derived type constructors")
endif()

set(SPAMM_DOCS TRUE
  CACHE LOGICAL "Build the spammpack API documentation.")

if(SPAMM_DOCS)
  include(FindDoxygen)
  if(DOXYGEN_FOUND)
    if(DOXYGEN_DOT_FOUND)
      set(HAVE_DOT "YES")
    endif()
    configure_file(Doxyfile.in Doxyfile)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE}
      COMMENT "Generating API documentation")
  else()
    message(STATUS "Can not build documentation")
  endif()
endif()

# Fortran support was added in 3.1.0.
if(${CMAKE_VERSION} VERSION_LESS 3.1)
  include(${CMAKE_SOURCE_DIR}/cmake-scripts/FindOpenMP.cmake)
else()
  include(FindOpenMP)
endif()

include(FindLAPACK)
if(LAPACK_FOUND)
  list(APPEND COMPILE_DEFINITIONS "LAPACK_FOUND")
endif()

include(FindPythonInterp)
if(NOT PYTHONINTERP_FOUND)
  message(FATAL "Need a working python interpreter")
endif()

unset(SPAMM_EXTRA_DEFINES)

if(NOT SPAMM_DEBUG_LEVEL)
  set(SPAMM_DEBUG_LEVEL 0)
else()
  if(SPAMM_DEBUG_LEVEL LESS 0)
    message(FATAL_ERROR "SPAMM_DEBUG_LEVEL has to be >= 0")
  endif()
  message(STATUS "Setting SPAMM_DEBUG_LEVEL to ${SPAMM_DEBUG_LEVEL}")
endif()

set(SPAMM_DEBUG_LEVEL ${SPAMM_DEBUG_LEVEL}
  CACHE STRING "The debugging output verbosity (valid range 0 <= LEVEL <= 2)")
list(APPEND SPAMM_EXTRA_DEFINES "SPAMM_DEBUG_LEVEL=${SPAMM_DEBUG_LEVEL}")

set(SPAMM_COUNTERS FALSE
  CACHE BOOL "Add non-zero and number_operations counters.")

if(SPAMM_COUNTERS)
  message(STATUS "Compiling with SpAMM counters")
  list(APPEND SPAMM_EXTRA_DEFINES SPAMM_COUNTERS)
endif()

mark_as_advanced(SPAMM_STORE_TRANSPOSE)
set(SPAMM_STORE_TRANSPOSE FALSE
  CACHE BOOL "Store the transpose sub-matrix to aid vectorization.")

if(SPAMM_STORE_TRANSPOSE)
  message(STATUS "Matrix will store the transpose")
  list(APPEND SPAMM_EXTRA_DEFINES SPAMM_STORE_TRANSPOSE)
endif()

if(HAVE_DEFERRED_STRING_LENGTH)
  list(APPEND SPAMM_EXTRA_DEFINES HAVE_DEFERRED_STRING_LENGTH)
endif()

if(LOG_FILENAMES)
  list(APPEND SPAMM_EXTRA_DEFINES LOG_FILENAMES)
endif()

if(HAVE_FINALIZE)
  list(APPEND SPAMM_EXTRA_DEFINES HAVE_FINALIZE)
endif()

if(HAVE_FINALIZE)
  list(APPEND SPAMM_EXTRA_DEFINES HAVE_CONSTRUCTOR)
endif()

set(VALGRIND_TESTS FALSE
  CACHE BOOL "Build the valgrind tests.")

if(VALGRIND_TESTS)
  find_program(VALGRIND valgrind)
  if(VALGRIND)
    message(STATUS "Setting up valgrind tests")
  endif()
endif()

set(SPAMM_PRECISION double
  CACHE STRING "The real precision of the library, single or double.")

string(TOLOWER ${SPAMM_PRECISION} SPAMM_PRECISION)

if(SPAMM_PRECISION STREQUAL "double")
  set(SPAMM_KIND_EXPRESSION "0d0")
else()
  set(SPAMM_KIND_EXPRESSION "0e0")
endif()
message(STATUS "Building ${SPAMM_PRECISION} precision version")
list(APPEND COMPILE_DEFINITIONS "SPAMM_KIND_EXPRESSION=kind(${SPAMM_KIND_EXPRESSION})")

set(SPAMM_BLOCK_SIZE 4
  CACHE STRING "The matrix size of basic submatrix blocks at the leaves.")
list(APPEND COMPILE_DEFINITIONS "SPAMM_BLOCK_SIZE=${SPAMM_BLOCK_SIZE}")

set(SPAMM_CHUNK_SIZE 32
  CACHE STRING "The matrix size of a chunk.")
list(APPEND COMPILE_DEFINITIONS "SPAMM_CHUNK_SIZE=${SPAMM_CHUNK_SIZE}")

mark_as_advanced(SPAMM_CHUNK_BLOCKS)
math(EXPR SPAMM_CHUNK_BLOCKS "${SPAMM_CHUNK_SIZE}/${SPAMM_BLOCK_SIZE}")
set(SPAMM_CHUNK_BLOCKS ${SPAMM_CHUNK_BLOCKS}
  CACHE STRING "The number of submatrix blocks in a chunk.")
list(APPEND COMPILE_DEFINITIONS "SPAMM_CHUNK_BLOCKS=${SPAMM_CHUNK_BLOCKS}")

mark_as_advanced(SPAMM_CHUNK_NODES)
math(EXPR SPAMM_CHUNK_NODES "(4*${SPAMM_CHUNK_BLOCKS}*${SPAMM_CHUNK_BLOCKS}-1)/3.")
set(SPAMM_CHUNK_NODES ${SPAMM_CHUNK_NODES}
  CACHE STRING "The number of tree nodes inside a chunk.")
list(APPEND COMPILE_DEFINITIONS "SPAMM_CHUNK_NODES=${SPAMM_CHUNK_NODES}")

find_program(CTAGS ctags)

if(CTAGS)
  add_custom_target(ctags
    ${CTAGS}
    --sort=foldcase
    --Fortran-kinds=+iL
    --recurse
    ${CMAKE_SOURCE_DIR}/src/*.F90
    ${CMAKE_SOURCE_DIR}/src/*.c
    ${CMAKE_SOURCE_DIR}/src/*.h)
endif()

find_program(ETAGS etags)

if(ETAGS)
  add_custom_target(etags
    ${ETAGS}
    ${CMAKE_SOURCE_DIR}/{src,tests,utilities,spammsand}/*.F90
    ${CMAKE_SOURCE_DIR}/{src,tests,utilities,spammsand}/*.c
    ${CMAKE_SOURCE_DIR}/{src,tests,utilities,spammsand}/*.h
    )
endif()

set(CPACK_PACKAGE_NAME "spammpack")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
  "An implementation of the Sparse Approximate Matrix Multiply (SpAMM) Algorithm")
set(CPACK_PACKAGE_CONTACT "www.freeon.org")

set(CPACK_GENERATOR DEB RPM)
set(CPACK_SOURCE_GENERATOR TBZ2)

set(CPACK_RPM_PACKAGE_LICENSE "BSD")
set(CPACK_RPM_PACKAGE_URL "http://www.freeon.org")

#add_subdirectory( src-solvers )
#add_subdirectory( tools )

add_subdirectory( src )
add_subdirectory( spammsand )

#add_subdirectory( tests )
#add_subdirectory( utilities )
