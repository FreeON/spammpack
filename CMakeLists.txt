project( spammpack Fortran )
cmake_minimum_required( VERSION 2.8 )

if( NOT CMAKE_Fortran_COMPILER )
  message( FATAL_ERROR "Can not find Fortran compiler" )
endif()

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

set( CPACK_DEBIAN_PACKAGE_HOMEPAGE "http://www.freeon.org/spammpack" )

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

include( CMakePackageConfigHelpers OPTIONAL RESULT_VARIABLE PACKAGE_CONFIG_FOUND )
if( PACKAGE_CONFIG_FOUND )
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
endif()

set( CMAKE_TRY_COMPILE_CONFIGURATION Debug )

if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
  set( CMAKE_Fortran_FLAGS_DEBUG
    "-O0 -g -gdwarf-3"
    "-fvar-tracking-assignments -fcheck=all -fsanitize=address"
    "-ffree-line-length-none"
    "-Wall -Wimplicit-interface -Wimplicit-procedure"
    )
  string( REGEX REPLACE ";" " " CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-Ofast -g -ffree-line-length-none" )
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
  set( CMAKE_Fortran_FLAGS_DEBUG
    "-stand f08 -O0 -g"
    "-extend-source"
    "-debug all -check all"
    "-warn unused -traceback"
    )
  string( REGEX REPLACE ";" " " CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -no-prec-div -fp-model fast=2 -vec-report -g" )
else()
  message( STATUS "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}" )
endif()

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release )
  message( WARNING
    "Setting build type to ${CMAKE_BUILD_TYPE}. Known "
    "build types are Debug and Release, which can be specified using the "
    "variable CMAKE_BUILD_TYPE."
    )
endif()

set( CMAKE_TRY_COMPILE_CONFIGURATION ${CMAKE_BUILD_TYPE} )

try_compile( COMPILE_RESULT_1
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/cmake-tests/error_stop.F90
  )
if( NOT COMPILE_RESULT_1 )
  message( FATAL_ERROR "The Fortran compiler does not understand the 'error stop' statement" )
endif()

try_run( RUN_RESULT_2 COMPILE_RESULT_2
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/cmake-tests/deferred_string_length.F90
  )
if( NOT COMPILE_RESULT_2 )
  message( FATAL_ERROR "The Fortran compiler can not compile deferred length string expressions" )
endif()

if( NOT RUN_RESULT_2 EQUAL 0 )
  message( WARNING
    "The Fortran compiler has trouble with deferred length string "
    "assignments, which might render screen output rather ugly. Please "
    "consider using a more recent compiler that fully implements "
    "this part of the the Fortran 2003 standard (e.g. >=gcc-4.8).")
else()
  set( HAVE_DEFERRED_STRING_LENGTH TRUE )
endif()

try_compile( COMPILE_RESULT_3
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/cmake-tests/long_lines.F90
  )
if( COMPILE_RESULT_3 )
  message( STATUS "The Fortran compiler supports lines > 132 characters" )
  set( LOG_FILENAMES TRUE )
endif()

include( FindDoxygen )
if( DOXYGEN_FOUND )
  if( DOXYGEN_DOT_FOUND )
    set( HAVE_DOT "YES" )
  endif()
  configure_file( Doxyfile.in Doxyfile )
  add_custom_target( doc ${DOXYGEN_EXECUTABLE} )
else()
  message( STATUS "Can not build documentation" )
endif()

include( ${CMAKE_SOURCE_DIR}/cmake-scripts/FindOpenMP.cmake )

include( FindLAPACK )

if( LAPACK_FOUND )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DLAPACK_FOUND" )
endif()

include( FindPythonInterp )
if( NOT PYTHONINTERP_FOUND )
  message( FATAL "Need a working python interpreter" )
endif()

if( NOT SPAMM_DEBUG_LEVEL )
  set( SPAMM_DEBUG_LEVEL 0 )
else()
  if( SPAMM_DEBUG_LEVEL LESS 0 )
    message( FATAL_ERROR "SPAMM_DEBUG_LEVEL has to be >= 0" )
  endif()
  message( STATUS "Setting SPAMM_DEBUG_LEVEL to ${SPAMM_DEBUG_LEVEL}" )
endif()

set( SPAMM_DEBUG_LEVEL ${SPAMM_DEBUG_LEVEL}
  CACHE STRING "The debugging output verbosity (valid range 0 <= LEVEL <= 2)" )
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DSPAMM_DEBUG_LEVEL=${SPAMM_DEBUG_LEVEL}" )

set( SPAMM_COUNTERS FALSE
  CACHE BOOL "Add non-zero and number_operations counters" )

set( SPAMM_STORE_TRANSPOSE FALSE
  CACHE BOOL "Store the transpose sub-matrix to aid vectorization" )

unset( SPAMM_EXTRA_DEFINES )

if( SPAMM_COUNTERS )
  list( APPEND SPAMM_EXTRA_DEFINES SPAMM_COUNTERS )
endif()

if( SPAMM_STORE_TRANSPOSE )
  list( APPEND SPAMM_EXTRA_DEFINES SPAMM_STORE_TRANSPOSE )
endif()

if( HAVE_DEFERRED_STRING_LENGTH )
  list( APPEND SPAMM_EXTRA_DEFINES HAVE_DEFERRED_STRING_LENGTH )
endif()

if( LOG_FILENAMES )
  list( APPEND SPAMM_EXTRA_DEFINES LOG_FILENAMES )
endif()

if( HAVE_DEFERRED_STRING_LENGTH )
  list( APPEND SPAMM_EXTRA_DEFINES HAVE_DEFERRED_STRING_LENGTH )
endif()

list( LENGTH SPAMM_EXTRA_DEFINES NUMBER_EXTRAS )
if( NUMBER_EXTRAS GREATER 0 )
  message( STATUS "SPAMM_EXTRA_DEFINES: ${SPAMM_EXTRA_DEFINES}" )
endif()

set( VALGRIND_TESTS FALSE
  CACHE BOOL "Build the valgrind tests" )

if( VALGRIND_TESTS )
  find_program( VALGRIND valgrind )
  if( VALGRIND )
    message( STATUS "Setting up valgrind tests" )
  endif()
endif()

find_program( CTAGS ctags )

if( CTAGS )
  add_custom_target( ctags
    ${CTAGS}
    --sort=foldcase
    --Fortran-kinds=+iL
    --recurse
    ${CMAKE_SOURCE_DIR}/src/*.F90
    ${CMAKE_SOURCE_DIR}/src/*.c
    ${CMAKE_SOURCE_DIR}/src/*.h
    )
endif()

enable_testing()

add_subdirectory( src-Fortran )
add_subdirectory( tests )
