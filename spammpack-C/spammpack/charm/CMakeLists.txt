INCLUDE( CheckFunctionExists )
INCLUDE( CheckLibraryExists )
INCLUDE( CheckCSourceCompiles )

PROJECT( spamm C CXX )
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake_config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

# Check for supported clocks. Thanks to Jonathan Lifflander
# <jliffl2@illinois.edu> for suggesting this test.
SET( HAVE_CLOCK_GETTIME "unset" )
CHECK_FUNCTION_EXISTS( clock_gettime HAVE_CLOCK_GETTIME )
IF( HAVE_CLOCK_GETTIME STREQUAL "unset" )
  CHECK_LIBRARY_EXISTS( rt clock_gettime "" HAVE_LIBRT )
  IF( NOT HAVE_LIBRT )
    MESSAGE( FATAL_ERROR "I can't find the clock_gettime function in the following libraries: ${LIBRARIES}" )
  ENDIF()
ENDIF()

SET( CLOCKTYPE "unset" )
FOREACH( CHECK_CLOCKTYPE CLOCK_MONOTONIC_RAW CLOCK_MONOTONIC CLOCK_REALTIME )
  SET( CLOCKTYPE "unset" )
  MESSAGE( "checking ${CHECK_CLOCKTYPE}" )
  CHECK_C_SOURCE_COMPILES( "#include <time.h>
main() { struct timespec p; clock_gettime(${CHECK_CLOCKTYPE}, &p); }" CLOCKTYPE)
  IF( NOT CLOCKTYPE STREQUAL "unset" )
    BREAK()
  ENDIF()
ENDFOREACH()

IF( CLOCKTYPE STREQUAL "unset" )
  MESSAGe( FATAL_ERROR "no known clocktype" )
ENDIF()
