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
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_algebra.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_c_bindings.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_convert.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_globals.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_management.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_project.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spamm_types.mod
  ${CMAKE_CURRENT_BINARY_DIR}/spammpack.mod
  )

set( LIBRARY_BASENAME "spammpack" )

include_directories( ${CMAKE_CURRENT_BINARY_DIR}/include-serial )

file( MAKE_DIRECTORY "include-serial" )
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
# Without these additional (chaining) dependencies, a parallel build gets all
# confused because it tries to reuse the .mod files in some weird way.
# Parallel builds are still possible within a library, but not across
# unfortunately.
add_dependencies( spammpack-serial-static spammpack-serial-shared )
set_target_properties( spammpack-serial-static
  PROPERTIES
  COMPILE_DEFINITIONS "${SPAMM_EXTRA_DEFINES}"
  OUTPUT_NAME ${LIBRARY_BASENAME}
  )

install(
  TARGETS spammpack-serial-static spammpack-serial-shared
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

file( COPY
  ${spammpack-modules}
  DESTINATION include-serial
  )

install(
  FILES
  ${CMAKE_CURRENT_BINARY_DIR}/include-serial/spamm_algebra.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include-serial/spamm_c_bindings.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include-serial/spamm_convert.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include-serial/spamm_globals.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include-serial/spamm_management.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include-serial/spamm_project.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include-serial/spamm_types.mod
  ${CMAKE_CURRENT_BINARY_DIR}/include-serial/spammpack.mod
  DESTINATION include-serial
  )

if( OPENMP_FOUND )
  include_directories( ${CMAKE_CURRENT_BINARY_DIR}/include-serial )

  file( MAKE_DIRECTORY "include-threaded" )
  add_library( spammpack-threaded-shared SHARED ${spammpack-sources} )
  add_dependencies( spammpack-threaded-shared spammpack-serial-static )
  set_target_properties( spammpack-threaded-shared
    PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    COMPILE_DEFINITIONS "${SPAMM_EXTRA_DEFINES}"
    COMPILE_FLAGS "${OpenMP_Fortran_FLAGS} ${MODULE_OUTPUT_ARGUMENT} ${CMAKE_CURRENT_BINARY_DIR}/include-threaded"
    SOVERSION 0
    VERSION 0.0.0
    OUTPUT_NAME ${LIBRARY_BASENAME}_threaded
    )

  add_library( spammpack-threaded-static STATIC ${spammpack-sources} )
  add_dependencies( spammpack-threaded-static spammpack-threaded-shared )
  set_target_properties( spammpack-threaded-static
    PROPERTIES
    COMPILE_DEFINITIONS "${SPAMM_EXTRA_DEFINES}"
    COMPILE_FLAGS "${OpenMP_Fortran_FLAGS} ${MODULE_OUTPUT_ARGUMENT} ${CMAKE_CURRENT_BINARY_DIR}/include-threaded"
    OUTPUT_NAME ${LIBRARY_BASENAME}_threaded
    )

  install(
    TARGETS spammpack-threaded-static spammpack-threaded-shared
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

  install(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/include-threaded/spamm_algebra.mod
    ${CMAKE_CURRENT_BINARY_DIR}/include-threaded/spamm_c_bindings.mod
    ${CMAKE_CURRENT_BINARY_DIR}/include-threaded/spamm_convert.mod
    ${CMAKE_CURRENT_BINARY_DIR}/include-threaded/spamm_globals.mod
    ${CMAKE_CURRENT_BINARY_DIR}/include-threaded/spamm_management.mod
    ${CMAKE_CURRENT_BINARY_DIR}/include-threaded/spamm_project.mod
    ${CMAKE_CURRENT_BINARY_DIR}/include-threaded/spamm_types.mod
    ${CMAKE_CURRENT_BINARY_DIR}/include-threaded/spammpack.mod
    DESTINATION include-threaded
    )

endif()

if( CTAGS )
  add_custom_target( ctags
    ${CTAGS}
    --sort=foldcase
    --Fortran-kinds=+iL
    ${spammpack-sources-Fortran}
    *.h
    )
endif()
