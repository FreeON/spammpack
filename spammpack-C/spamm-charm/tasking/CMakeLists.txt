cmake_minimum_required( VERSION 2.8 )
project( tasking C )

add_executable( tasking-cpu tasking-cpu.c )
add_executable( tasking-memory tasking-memory.c )

if( OPENMP )
  include( FindOpenMP )
  if( OPENMP_FOUND )
    set_target_properties( tasking-cpu
      PROPERTIES
      COMPILE_FLAGS "${TEMP_C_FLAGS}"
      LINK_FLAGS "${TEMP_LINK_FLAGS}"
      )
    set_target_properties( tasking-memory
      PROPERTIES
      COMPILE_FLAGS "${TEMP_C_FLAGS}"
      LINK_FLAGS "${TEMP_LINK_FLAGS}"
      )
  endif()
endif()
