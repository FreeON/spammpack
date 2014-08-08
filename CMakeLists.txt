project( spammpack Fortran C )
cmake_minimum_required( VERSION 2.8 )

set( CPACK_PACKAGE_NAME "spammpack" )
set( CPACK_PACKAGE_VERSION_MAJOR 1 )
set( CPACK_PACKAGE_VERSION_MINOR 0 )
set( CPACK_PACKAGE_VERSION_PATCH 1 )
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
  "/.git/"
  "/CMakeFiles/"
  "/_CPack_Packages/"
  "tests/testmatrix.*$"
  "/Makefile$"
  "/.*\\\\.cmake$"
  "/CMakeCache\\\\.txt$"
  "/${CPACK_PACKAGE_NAME}.*tar.*$"
  "/${CPACK_PACKAGE_NAME}.*rpm$"
  "/${CPACK_PACKAGE_NAME}.*deb$"
  "/install_manifest\\\\.txt"
  "\\\\.mod$"
  "/lib.*\\\\.so.*$"
  "/lib.*\\\\.a$"
  "/Doxyfile$"
  "/latex/"
  "/html/"
  "/man/"
  )

message( STATUS "ignoring ${CPACK_SOURCE_IGNORE_FILES}" )

include( CPack )

if( CMAKE_C_COMPILER_ID STREQUAL "GNU" )
  set( CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall" )
elseif( CMAKE_C_COMPILER_ID STREQUAL "Intel" )
  set( CMAKE_C_FLAGS_DEBUG "-O0 -g -check=conversions,stack,uninit" )
else()
  message( STATUS "Unknown C compiler ${CMAKE_C_COMPILER_ID}" )
endif()

if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
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

include( FindDoxygen )
if( DOXYGEN_FOUND )
  add_custom_target( docs doxygen )
  configure_file( Doxyfile.in Doxyfile )
endif()

include( FindOpenMP )

list( APPEND CMAKE_C_FLAGS "-DTIMER_GETTIMEOFDAY" )

add_subdirectory( src )
add_subdirectory( tests )
