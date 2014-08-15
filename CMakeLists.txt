project( spammpack Fortran C )
cmake_minimum_required( VERSION 2.8 )

set( CPACK_PACKAGE_NAME "spammpack" )
set( CPACK_PACKAGE_VERSION_MAJOR 1 )
set( CPACK_PACKAGE_VERSION_MINOR 0 )
set( CPACK_PACKAGE_VERSION_PATCH 0 )
set( CPACK_PACKAGE_DESCRIPTION_SUMMARY
  "An implementation of the Sparse Approximate Matrix Multiply (SpAMM) Algorithm"
  )
set( CPACK_PACKAGE_CONTACT "www.freeon.org" )

set( CPACK_GENERATOR DEB RPM )
set( CPACK_SOURCE_GENERATOR TBZ2 )

set( CPACK_RPM_PACKAGE_LICENSE "BSD" )
set( CPACK_RPM_PACKAGE_URL "http://www.freeon.org" )

set( CPACK_DEBIAN_PACKAGE_HOMEPAGE "http://www.freeon.org" )

list( APPEND CPACK_SOURCE_IGNORE_FILES
  "/${CPACK_PACKAGE_NAME}.*deb$"
  "/${CPACK_PACKAGE_NAME}.*rpm$"
  "/${CPACK_PACKAGE_NAME}.*tar.bz2$"
  "/.*[.]cmake$"
  "/CMakeCache[.]txt$"
  "/CMakeFiles/"
  "/Doxyfile$"
  "/Makefile$"
  "/[.]git/"
  "/_CPack_Packages/"
  "/html/"
  "/install_manifest[.]txt"
  "/latex/"
  "/lib.*[.]a$"
  "/lib.*[.]so.*$"
  "/man/"
  "[.]mod$"
  "tests/testmatrix.*$"
  )

include( CPack )

if( CMAKE_C_COMPILER_ID STREQUAL "GNU" )
  set( CMAKE_C_FLAGS_DEBUG "-O0 -g -gdwarf-3 -fvar-tracking-assignments -Wall" )
elseif( CMAKE_C_COMPILER_ID STREQUAL "Intel" )
  set( CMAKE_C_FLAGS_DEBUG "-O0 -g -check=conversions,stack,uninit" )
else()
  message( STATUS "Unknown C compiler ${CMAKE_C_COMPILER_ID}" )
endif()

if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -gdwarf-3 -fvar-tracking-assignments -fcheck=all -Wall" )
elseif( CMAKE_Fotran_COMPILER_ID STREQUAL "Intel" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-stand f08 -O0 -g -debug all -check all -traceback -assume realloc_lhs -heap-arrays" )
else()
  message( STATUS "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}" )
endif()

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug )
  message( STATUS
    "Setting build type to ${CMAKE_BUILD_TYPE}. If this was not "
    "intended, please specify another one using the variable "
    "CMAKE_BUILD_TYPE."
    )
endif()

include( CheckFortranSourceCompiles.cmake )

set( Fortran_TEST_SOURCE
"
program test
error stop
end program test
"
)

check_fortran_source_compiles( ${Fortran_TEST_SOURCE} RESULT_TEST )
if( NOT RESULT_TEST )
  message( FATAL_ERROR "The Fortran compiler does not understand the 'error stop' statement" )
endif()

include( FindDoxygen )
if( DOXYGEN_FOUND )
  add_custom_target( docs doxygen )
  configure_file( Doxyfile.in Doxyfile )
endif()

include( ${CMAKE_SOURCE_DIR}/FindOpenMP.cmake )

list( APPEND
  CMAKE_C_FLAGS "-DTIMER_GETTIMEOFDAY"
  )

if( NOT SPAMM_SINGLE )
  set( SPAMM_SINGLE FALSE
    CACHE BOOL "Build the single precision version"
    )
else()
  set( SPAMM_SINGLE TRUE
    CACHE BOOL "Build the single precision version"
    )
  list( APPEND CMAKE_Fortran_FLAGS "-DSPAMM_SINGLE" )
endif()

find_program( VALGRIND valgrind )
if( VALGRIND )
  message( STATUS "Setting up valgrind tests" )
endif()

enable_testing()

add_subdirectory( src )
add_subdirectory( tests )
