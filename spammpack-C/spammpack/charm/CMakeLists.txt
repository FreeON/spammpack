include( CheckFunctionExists )
include( CheckCSourceCompiles )

project( spamm C CXX )
cmake_minimum_required( VERSION 2.8 )

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/cmake_config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

# Check for supported clocks. Thanks to Jonathan Lifflander
# <jliffl2@illinois.edu> for suggesting this test.
check_function_exists( clock_gettime HAVE_CLOCK_GETTIME )
message( "${HAVE_CLOCK_GETTIME}" )
if( ${HAVE_CLOCK_GETTIME} )
  set( HAVE_CLOCK_GETTIME )
  check_library_exists( rt HAVE_LIBRT )
  if( not ${HAVE_LIBRT} )
    message( FATAL_ERROR "I can't find the clock_gettime function in the
    following libraries: ${LIBRARIES}" )
  endif()
  check_function_exists( clock_gettime HAVE_CLOCK_GETTIME )
endif()

set( CLOCKTYPE "unset" )
foreach( CHECK_CLOCKTYPE CLOCK_MONOTONIC_RAW CLOCK_MONOTONIC CLOCK_REALTIME )
  check_C_source_compiles( "main() { struct timespec p;
  clock_gettime(${CLOCKTYPE}, &p); }" CLOCKTYPE)
endforeach()
if( CLOCKTYPE strequal "unset" )
  message( FATAL_ERROR "no known clocktype" )
endif()
