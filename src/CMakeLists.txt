set( spammpack-sources
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

set( LIBRARY_BASENAME "spammpack" )

add_library( spammpack-serial-shared SHARED ${spammpack-sources} )
set_target_properties( spammpack-serial-shared
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  COMPILE_DEFINITIONS "${SPAMM_EXTRA_DEFINES}"
  SOVERSION 0
  VERSION 0.0.0
  OUTPUT_NAME ${LIBRARY_BASENAME}
 )

add_library( spammpack-serial-static STATIC ${spammpack-sources} )
set_target_properties( spammpack-serial-static
  PROPERTIES
  OUTPUT_NAME ${LIBRARY_BASENAME}
  )

install(
  TARGETS spammpack-serial-static spammpack-serial-shared
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

if( OPENMP_FOUND )
  add_library( spammpack-threaded-shared SHARED ${spammpack-sources} )
  set_target_properties( spammpack-threaded-shared
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    POSITION_INDEPENDENT_CODE TRUE
    SOVERSION 0
    VERSION 0.0.0
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    OUTPUT_NAME ${LIBRARY_BASENAME}_threaded
    )

  add_library( spammpack-threaded-static STATIC ${spammpack-sources} )
  set_target_properties( spammpack-threaded-static
    PROPERTIES
    OUTPUT_NAME ${LIBRARY_BASENAME}_threaded
    )

  install(
    TARGETS spammpack-threaded-static spammpack-threaded-shared
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
    ${CTAGS}
    --sort=foldcase
    --Fortran-kinds=+iL
    ${spammpack-sources-Fortran}
    *.h
    )
endif()
