project( chunklib C ASM )
cmake_minimum_required( VERSION 2.8 )

include( FindOpenMP )
include( CheckCSourceCompiles )
include( CheckCCompilerFlag )
include( CheckFunctionExists )

check_c_compiler_flag( "-std=gnu99" C99_CONFORM )
if( NOT C99_CONFORM )
  message( WARNING "Can not set C compiler to C99 standard" )
else()
  set( STDC "-std=gnu99" )
endif()

check_function_exists( clock_gettime FOUND_GETTIME )
if( NOT FOUND_GETTIME )
  message( WARNING "Could not find gettime" )
endif()

check_function_exists( ceil FOUND_CEIL )
if( NOT FOUND_CEIL )
  # Try linking with m.
  set( CMAKE_REQUIRED_LIBRARIES "-lm" )
  check_function_exists( ceil FOUND_CEIL_WITH_M )
  if( NOT FOUND_CEIL_WITH_M )
    message( FATAL_ERROR "Could not find ceil() function" )
  endif()
  set( MATHLIB "-lm" )
endif()

# Check for supported clocks. Thanks to Jonathan Lifflander
# <jliffl2@illinois.edu> for suggesting this test.
foreach( CLOCKTYPE_
    CLOCK_MONOTONIC_RAW
    CLOCK_MONOTONIC
    CLOCK_REALTIME )
  check_c_source_compiles(
    "#include <time.h>
    int main () {
    struct timespec p;
    clock_gettime(${CLOCKTYPE_}, &p);
    return 0;
    }"
    TEST_COMPILED )
  if( TEST_COMPILED )
    message( STATUS "Setting clocktype to ${CLOCKTYPE_}" )
    set( CLOCKTYPE ${CLOCKTYPE_}
      CACHE STRING "The default clock type for clock_gettime()" )
    break()
  endif()
endforeach()

if( NOT TEST_COMPILED )
  message( FATAL_ERROR "could not find a suitable clock type" )
endif()

set( ENABLE_MIC FALSE
  CACHE BOOL "Enable Xeon Phi target" )
if( ENABLE_MIC )
  set( MIC_ALLOC TRUE )
  set( MIC_CFLAGS "-mmic" )
  set( MIC_LDFLAGS "-static" )
endif()

set( FORCE_OMP_NUM_THREADS FALSE
  CACHE STRING
  "Set the number of threads in OpenMP parallel sections to N" )

set( CHUNK_IMPLEMENTATION tree
  CACHE STRING "The chunk implementation (tiled, tree), default = tree" )

set( BLOCK_MULTIPLY block
  CACHE STRING
  "The dense multiply implementation (block, blas), default = block" )

if( BLOCK_MULTIPLY STREQUAL "block" )
  set( BLOCK_MULTIPLY TRUE )
elseif( BLOCK_MULTIPLY STREQUAL "blas" )
  set( BLOCK_BLAS TRUE )
else()
  message( FATAL_ERROR "Unknown multiply implementation" )
endif()

set( CHUNK_TREE_MAX_TIER 6
  CACHE STRING
  "The maximum tier at which new tasks are created when recursing the chunk
  tree" )

set( CHUNK_BLOCK_TRANSPOSE 32
  CACHE STRING
  "The minimum size N of the block matrix that will transpose a matrix to
  allow for loop vectorization" )

set( BLOCK_SLEEP 0
  CACHE STRING "Number of nanoseconds to sleep after block multiply" )

set( MEASURE_COMPLEXITY FALSE
  CACHE BOOL "Measure the product complexity" )

set( DEBUG_OUTPUT FALSE
  CACHE BOOL "Enable a lot of extra output" )

set( DEBUG_OUTPUT_MATRIX FALSE
  CACHE BOOL "Print the matrices for debugging" )

configure_file( config.h.in config.h )

enable_testing()

add_subdirectory( src )
add_subdirectory( tests )
