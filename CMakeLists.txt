cmake_minimum_required(VERSION 2.8)
project(spammpack Fortran)

# Set some default compiler flags.
set(SPAMM_Fortran_FLAGS_GNU_DEBUG "-O0 -g")
set(SPAMM_Fortran_FLAGS_GNU_RELEASE "-O2 -g")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG ${SPAMM_Fortran_FLAGS_GNU_DEBUG})
  set(CMAKE_Fortran_FLAGS_RELEASE ${SPAMM_Fortran_FLAGS_GNU_RELEASE})
else()
  message(STATUS "Unknown comiler ${CMAKE_Fortran_COMPILER_ID}, using CMake defaults for compiler flags")
endif()

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(src)

set(SPAMM_TESTS FALSE CACHE BOOL "Compile test suite")
if(SPAMM_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
