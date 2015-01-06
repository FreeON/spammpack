set( spammpack-sources
  spamm_algebra.F90
  spamm_convert.F90
  spamm_globals.F90
  spamm_inverse.F90
  spamm_management.F90
  spamm_order_1.F90
  spamm_order_1_double.F90
  spamm_order_1_single.F90
  spamm_project.F90
  spamm_types.F90
  spamm_utilities.F90
  spammpack.F90
  )

set( spammpack-modules
  spamm_algebra.mod
  spamm_convert.mod
  spamm_globals.mod
  spamm_management.mod
  spamm_project.mod
  spamm_types.mod
  spammpack.mod
  spamm_utility_macros.h
  )

set( LIBRARY_BASENAME "spammpack" )

set_directory_properties(
  PROPERTIES
  COMPILE_DEFINITIONS "${SPAMM_EXTRA_DEFINES}"
  )

include( ${CMAKE_SOURCE_DIR}/cmake-scripts/BuildLibrary.cmake )

build_library( FALSE "" )

if( OPENMP_FOUND )
  build_library( TRUE "${OpenMP_Fortran_FLAGS}")
  add_dependencies( spammpack-threaded-shared spammpack-serial-static )
endif()
