project( spammpack Fortran C )
cmake_minimum_required( VERSION 2.8 )

set( CPACK_PACKAGE_NAME "spammpack" )
set( CPACK_PACKAGE_VERSION_MAJOR 1 )
set( CPACK_PACKAGE_VERSION_MINOR 1 )
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

set( INCLUDE_INSTALL_DIR include )
set( LIB_INSTALL_DIR lib )

include( CMakePackageConfigHelpers )
configure_package_config_file( spammpackConfig.cmake.in spammpackConfig.cmake
  INSTALL_DESTINATION ${LIB_INSTALL_DIR}/spammpack/cmake
  PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR
  )
write_basic_package_version_file( spammpackConfigVersion.cmake
  VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}
  COMPATIBILITY SameMajorVersion
  )
install( FILES
  ${CMAKE_BINARY_DIR}/spammpackConfig.cmake
  ${CMAKE_BINARY_DIR}/spammpackConfigVersion.cmake
  DESTINATION ${LIB_INSTALL_DIR}/spammpack/cmake
  )

if( CMAKE_C_COMPILER_ID STREQUAL "GNU" )
  set( CMAKE_C_FLAGS_DEBUG "-O0 -g -gdwarf-3 -fvar-tracking-assignments -Wall" )
elseif( CMAKE_C_COMPILER_ID STREQUAL "Intel" )
  set( CMAKE_C_FLAGS_DEBUG "-O0 -g -check=conversions,stack,uninit" )
else()
  message( STATUS "Unknown C compiler ${CMAKE_C_COMPILER_ID}" )
endif()

if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -gdwarf-3 -fvar-tracking-assignments -fcheck=all -Wall -ffree-line-length-none" )
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-stand f08 -O0 -g -debug all -check all -warn unused -traceback -assume realloc_lhs -heap-arrays" )
else()
  message( STATUS "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}" )
endif()

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release )
  message( STATUS
    "Setting build type to ${CMAKE_BUILD_TYPE}. Known build types "
    "are Debug and Release, which can be specified using the variable "
    "CMAKE_BUILD_TYPE."
    )
endif()

if( CMAKE_BUILD_TYPE STREQUAL Debug )
  message( STATUS "C compiler flags: ${CMAKE_C_FLAGS_DEBUG}" )
  message( STATUS "Fortran compiler flags: ${CMAKE_Fortran_FLAGS_DEBUG}" )
endif()

try_compile( COMPILE_RESULT_1
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/error_stop.F90
  )
if( NOT COMPILE_RESULT_1 )
  message( FATAL_ERROR "The Fortran compiler does not understand the 'error stop' statement" )
endif()

try_run( RUN_RESULT_2 COMPILE_RESULT_2
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/deferred_string_length.F90
  )
if( NOT COMPILE_RESULT_2 )
  message( FATAL_ERROR "The Fortran compiler can not compile deferred length string expressions" )
endif()
if( NOT RUN_RESULT_2 EQUAL 0 )
  message( WARNING
    "The Fortran compiler has trouble with deferred length string "
    "assignments, which might render screen output rather ugly. Please "
    "consider using a more recent compiler that fully implements "
    "the Fortran 2003 standard.")
endif()

include( FindDoxygen )
if( DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND )
  add_custom_target( docs ${DOXYGEN_EXECUTABLE} )
  configure_file( Doxyfile.in Doxyfile )
else()
  message( STATUS "Can not build documentation" )
endif()

include( ${CMAKE_SOURCE_DIR}/FindOpenMP.cmake )

if( NOT DEBUG_LEVEL )
  set( DEBUG_LEVEL 0
    CACHE STRING "The debugging output verbosity"
    )
else()
  if( DEBUG_LEVEL LESS 0 )
    message( FATAL_ERROR "DEBUG_LEVEL has to be >= 0" )
  endif()
  message( STATUS "Setting DEBUG_LEVEL to ${DEBUG_LEVEL}" )
endif()
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DDEBUG_LEVEL=${DEBUG_LEVEL}" )

# For now hard-coded.
set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DTIMER_GETTIMEOFDAY" )

if( NOT VALGRIND_TESTS )
  set( VALGRIND_TESTS FALSE
    CACHE BOOL "Build the valgrind tests"
    )
endif()

if( VALGRIND_TESTS )
  find_program( VALGRIND valgrind )
  if( VALGRIND )
    message( STATUS "Setting up valgrind tests" )
  endif()
endif()

find_program( CTAGS ctags )

enable_testing()

add_subdirectory( src )
add_subdirectory( tests )
