cmake_minimum_required( VERSION 2.8 )
project( tasking C )

message( STATUS "Found C compiler ${CMAKE_C_COMPILER_ID}" )
if( CMAKE_C_COMPILER_ID STREQUAL "GNU" )
  list( APPEND EXTRA_FLAGS "-std=c99" )
elseif( CMAKE_C_COMPILER_ID STREQUAL "Intel" )
  list( APPEND EXTRA_FLAGS "-std=c99" )
endif()

add_executable( tasking-cpu tasking-cpu.c )
add_executable( tasking-memory tasking-memory.c )

set( OPENMP TRUE CACHE BOOL "Whether to compile using OpenMP" )

if( OPENMP )
  include( FindOpenMP )
  if( OPENMP_FOUND )
    list( APPEND EXTRA_FLAGS ${OpenMP_FLAGS} )
  endif()
endif()

target_compile_definitions( tasking-cpu PUBLIC "_POSIX_C_SOURCE=199309L" )
target_compile_options( tasking-cpu PUBLIC ${EXTRA_FLAGS} )

target_compile_definitions( tasking-memory PUBLIC "_POSIX_C_SOURCE=199309L" )
target_compile_options( tasking-memory PUBLIC ${EXTRA_FLAGS} )
