cmake_minimum_required( VERSION 2.8 )

include( CheckLibraryExists )
include( CheckCSourceCompiles )

project( spamm-charm C CXX )
set( PACKAGE_VERSION 1.0 )

set( CHARMC "charmc" CACHE STRING "The Charm++ compiler command" )

check_library_exists( rt clock_gettime "" HAVE_CLOCK_GETTIME )
if( HAVE_CLOCK_GETTIME STREQUAL "unset" )
  message( FATAL_ERROR "I need clock_gettime()" )
endif()
list( APPEND CMAKE_REQUIRED_LIBRARIES "-lrt" )

# Check for supported clocks. Thanks to Jonathan Lifflander
# <jliffl2@illinois.edu> for suggesting this test.
foreach( TRY_CLOCKTYPE
    CLOCK_MONOTONIC_RAW
    CLOCK_MONOTONIC
    CLOCK_REALTIME
    )

  unset( HAVE_CLOCKTYPE CACHE )
  message( STATUS "checking clocktype ${TRY_CLOCKTYPE}" )
  check_c_source_compiles( "#include <time.h>
int main() {
  struct timespec p;
  clock_gettime(${TRY_CLOCKTYPE}, &p);
  return 0;
}"
    HAVE_CLOCKTYPE)

  if( HAVE_CLOCKTYPE )
    SET( CLOCKTYPE ${TRY_CLOCKTYPE} )
    break()
  endif()
endforeach()

if( NOT HAVE_CLOCKTYPE )
  message( FATAL_ERROR "no known clocktype" )
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake_config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h
  )
