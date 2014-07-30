set( spammpack-sources
  SpAMM_ALGEBRA.F90
  SpAMM_CONVERT.F90
  SpAMM_DERIVED.F90
  SpAMM_GLOBALS.F90
  SpAMM_MNGMENT.F90
  SpAMM_PROJECT.F90
  spamm_exit.c
  spamm_get_time.c
  spamm_trap.c
  spammpack.F90
  )

add_library( spammpack_serial_static STATIC ${spammpack-sources} )
add_library( spammpack_serial_shared SHARED ${spammpack-sources} )
set_target_properties( spammpack_serial_shared
  PROPERTIES
  SOVERSION 0
  VERSION 0.0.0
  )

if( OPENMP_FOUND )
  add_library( spammpack_threaded_static STATIC ${spammpack-sources} )
  add_library( spammpack_threaded_shared SHARED ${spammpack-sources} )
  set_target_properties( spammpack_threaded_static
    PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS}
    )
  set_target_properties( spammpack_threaded_shared
    PROPERTIES
    SOVERSION 0
    VERSION 0.0.0
    COMPILE_FLAGS ${OpenMP_C_FLAGS}
    )
endif()
