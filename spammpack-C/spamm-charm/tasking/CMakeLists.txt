cmake_minimum_required( VERSION 2.8 )
project( tasking C )

add_executable( tasking tasking.c )

if( OPENMP )
  include( FindOpenMP )
  if( OPENMP_FOUND )
    get_target_property( TEMP_C_FLAGS tasking COMPILE_FLAGS )
    if( TEMP_C_FLAGS )
      set( TEMP_C_FLAGS "${TEMP_C_FLAGS} ${OpenMP_C_FLAGS}" )
    else()
      set( TEMP_C_FLAGS "${OpenMP_C_FLAGS}" )
    endif()
    get_target_property( TEMP_LINK_FLAGS tasking LINK_FLAGS )
    if( TEMP_LINK_FLAGS )
      set( TEMP_LINK_FLAGS "${TEMP_LINK_FLAGS} ${OpenMP_C_FLAGS}" )
    else()
      set( TEMP_LINK_FLAGS "${OpenMP_C_FLAGS}" )
    endif()
    set_target_properties( tasking
      PROPERTIES
      COMPILE_FLAGS "${TEMP_C_FLAGS}"
      LINK_FLAGS "${TEMP_LINK_FLAGS}"
      )
  endif()
endif()