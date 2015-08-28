cmake_minimum_required(VERSION 2.8)
project(productspace C)

find_program(CYTHON_EXECUTABLE cython)
if(NOT CYTHON_EXECUTABLE)
  message(FATAL_ERROR "need cython to continue")
endif()

set(Python_ADDITIONAL_VERSIONS 2.7)
include(FindPythonLibs)

include_directories(${PYTHON_INCLUDE_DIRS})

set(THREADS_PREFER_PTHREAD_FLAG TRUE)
include(FindThreads)

set(PYTHON_SOURCES src/productspace.pyx)
add_custom_command(OUTPUT productspacemodule.c
  COMMAND ${CYTHON_EXECUTABLE} -a ${PYTHON_SOURCES} -o productspacemodule.c
  MAIN_DEPENDENCY ${PYTHON_SOURCES})

add_library(productspacemodule SHARED
  productspacemodule.c)
set_target_properties(productspacemodule
  PROPERTIES COMPILE_FLAGS ${CMAKE_THREAD_LIBS_INIT}
  LINK_FLAGS ${CMAKE_THREAD_LIBS_INIT}
  PREFIX "")