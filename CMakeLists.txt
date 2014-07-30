project( spammpack Fortran C )
cmake_minimum_required( VERSION 2.8 )

if( CMAKE_C_COMPILER_ID STREQUAL "GNU" )
  set( CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall" )
elseif( CMAKE_C_COMPILER_ID STREQUAL "Intel" )
  set( CMAKE_C_FLAGS_DEBUG "-O0 -g -check=conversions,stack,uninit" )
else()
  message( STATUS "Unknown C compiler ${CMAKE_C_COMPILER_ID}" )
endif()

if( CMAKE_Fotran_COMPILER_ID STREQUAL "GNU" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fcheck=all" )
elseif( CMAKE_Fotran_COMPILER_ID STREQUAL "Intel" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-stand f08 -O0 -g -debug all -check all -traceback -assume realloc_lhs -heap-arrays" )
else()
  message( STATUS "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}" )
endif()

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug )
  message( STATUS
    "Setting build type to ${CMAKE_BUILD_TYPE}. If this is not "
    "what was intended, please specify another one using the "
    "variable CMAKE_BUILD_TYPE."
    )
endif()

include( FindOpenMP )

list( APPEND CMAKE_C_FLAGS "-DTIMER_GETTIMEOFDAY" )

add_subdirectory( src )
add_subdirectory( tests )
