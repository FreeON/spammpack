set( spammpack-sources
  SpAMM_ALGEBRA.F90
  SpAMM_CONVERT.F90
  SpAMM_PROJECT.F90
  spamm_exit.c
  spamm_get_time.c
  spamm_globals.F90
  spamm_management.F90
  spamm_trap.c
  spamm_types.F90
  spammpack.F90
  )

add_library( spammpack_serial_static STATIC ${spammpack-sources} )
add_library( spammpack_serial_shared SHARED ${spammpack-sources} )
set_target_properties( spammpack_serial_static
  PROPERTIES
  OUTPUT_NAME spammpack
  )
set_target_properties( spammpack_serial_shared
  PROPERTIES
  SOVERSION 0
  VERSION 0.0.0
  OUTPUT_NAME spammpack
  )

if( OPENMP_FOUND )
  add_library( spammpack_threaded_static STATIC ${spammpack-sources} )
  add_library( spammpack_threaded_shared SHARED ${spammpack-sources} )
  set_target_properties( spammpack_threaded_static
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS}
    OUTPUT_NAME spammpack-threaded
    )
  set_target_properties( spammpack_threaded_shared
    PROPERTIES
    SOVERSION 0
    VERSION 0.0.0
    COMPILE_FLAGS ${OpenMP_C_FLAGS}
    OUTPUT_NAME spammpack-threaded
    )
  install(
    TARGETS spammpack_threaded_static spammpack_threaded_shared
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
endif()

install(
  TARGETS spammpack_serial_static spammpack_serial_shared
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
install(
  FILES 
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_algebra.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_convert.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_types.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_globals.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_mngment.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_project.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spammpack.mod
  DESTINATION include
  )

add_custom_target( ctags
  ctags *.F90 *.c
  )
