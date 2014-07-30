project( spammpack Fortran C )
cmake_minimum_required( VERSION 2.8 )

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug )
  message( STATUS
    "Setting build type to ${CMAKE_BUILD_TYPE}. If this is not "
    "what was intended, please specify another one using the "
    "variable CMAKE_BUILD_TYPE."
    )
endif()

include( FindOpenMP )

list( APPEND
  CMAKE_C_FLAGS
  ${OpenMP_C_FLAGS}
  )

list( APPEND
  CMAKE_Fortran_FLAGS
  ${OpenMP_C_FLAGS}
  )

add_subdirectory( src )
