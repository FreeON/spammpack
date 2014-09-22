set( spammpack-sources-C
  spamm_exit.c
  spamm_get_time.c
  spamm_trap.c
  )

set( spammpack-sources-Fortran
  spamm_algebra.F90
  spamm_c_bindings.F90
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

if( SPAMM_SINGLE )
  set( LIBRARY_BASENAME "spammpack_single" )
else()
  set( LIBRARY_BASENAME "spammpack" )
endif()

add_library( spammpack-serial-C SHARED OBJECT ${spammpack-sources-C} )
set_target_properties( spammpack-serial-C
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  )

add_library( spammpack-serial-Fortran SHARED OBJECT ${spammpack-sources-Fortran} )
set_target_properties( spammpack-serial-Fortran
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  )

add_library( spammpack_serial_static
  STATIC
  $<TARGET_OBJECTS:spammpack-serial-C>
  $<TARGET_OBJECTS:spammpack-serial-Fortran>
  )
set_target_properties( spammpack_serial_static
  PROPERTIES
  OUTPUT_NAME ${LIBRARY_BASENAME}
  )

add_library( spammpack_serial_shared
  SHARED
  $<TARGET_OBJECTS:spammpack-serial-C>
  $<TARGET_OBJECTS:spammpack-serial-Fortran>
  )
set_target_properties( spammpack_serial_shared
  PROPERTIES
  SOVERSION 0
  VERSION 0.0.0
  OUTPUT_NAME ${LIBRARY_BASENAME}
  )

install(
  TARGETS spammpack_serial_static spammpack_serial_shared
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

if( OPENMP_FOUND )
  add_library( spammpack_threaded-C OBJECT ${spammpack-sources-C} )
  set_target_properties( spammpack_threaded-C
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_C_FLAGS}
    POSITION_INDEPENDENT_CODE TRUE
    )

  add_library( spammpack_threaded-Fortran OBJECT ${spammpack-sources-Fortran} )
  set_target_properties( spammpack_threaded-Fortran
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    POSITION_INDEPENDENT_CODE TRUE
    )

  add_library( spammpack_threaded_static
    STATIC
    $<TARGET_OBJECTS:spammpack-serial-C>
    $<TARGET_OBJECTS:spammpack-serial-Fortran>
    )
  set_target_properties( spammpack_threaded_static
    PROPERTIES
    OUTPUT_NAME ${LIBRARY_BASENAME}_threaded
    )

  add_library( spammpack_threaded_shared
    SHARED
    $<TARGET_OBJECTS:spammpack-serial-C>
    $<TARGET_OBJECTS:spammpack-serial-Fortran>
    )
  set_target_properties( spammpack_threaded_shared
    PROPERTIES
    SOVERSION 0
    VERSION 0.0.0
    COMPILE_FLAGS ${OpenMP_C_FLAGS}
    OUTPUT_NAME ${LIBRARY_BASENAME}_threaded
    )

  install(
    TARGETS spammpack_threaded_static spammpack_threaded_shared
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
endif()

install(
  FILES
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_algebra.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_c_bindings.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_convert.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_globals.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_management.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_project.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_types.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spammpack.mod
  DESTINATION include
  )

if( CTAGS )
  add_custom_target( ctags
    ${CTAGS} --sort=foldcase --Fortran-kinds=+iL ${spammpack-sources-C} ${spammpack-sources-Fortran} *.h
    )
endif()
