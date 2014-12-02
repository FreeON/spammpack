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
  )

set( LIBRARY_BASENAME "spammpack" )

set_directory_properties(
  PROPERTIES
  COMPILE_DEFINITIONS "${SPAMM_EXTRA_DEFINES}"
  )

add_library( spammpack-serial-shared SHARED ${spammpack-sources} )
set_target_properties( spammpack-serial-shared
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  SOVERSION 0
  VERSION 0.0.0
  OUTPUT_NAME ${LIBRARY_BASENAME}
 )

add_library( spammpack-serial-static STATIC ${spammpack-sources} )
# Without these additional (chaining) dependencies, a parallel build gets all
# confused because of the way CMake builds Fortran sources: The compilation is
# done in the src dir, but uses -o to write the .o file into a subdirectory
# named after the target. CMake does not use -J or something similar to place
# the generated .mod file in that directory as well. Instead it copies the
# .mod file from the src directory into the target directory.  Parallel builds
# are still possible within a library, but not across unfortunately.
add_dependencies( spammpack-serial-static spammpack-serial-shared )
set_target_properties( spammpack-serial-static
  PROPERTIES
  OUTPUT_NAME ${LIBRARY_BASENAME}
  )

install(
  TARGETS spammpack-serial-static spammpack-serial-shared
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

file( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include-serial" )
add_custom_command( TARGET spammpack-serial-static
  POST_BUILD
  COMMAND cp ${spammpack-modules} "include-serial"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

unset( serial-modfiles )
foreach( mod ${spammpack-modules} )
  list( APPEND serial-modfiles
    ${CMAKE_CURRENT_BINARY_DIR}/include-serial/${mod}
    )
endforeach()

install(
  FILES ${serial-modfiles}
  DESTINATION include-serial
  )

if( OPENMP_FOUND )
  add_library( spammpack-threaded-shared SHARED ${spammpack-sources} )
  add_dependencies( spammpack-threaded-shared spammpack-serial-static )
  set_target_properties( spammpack-threaded-shared
    PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
    SOVERSION 0
    VERSION 0.0.0
    OUTPUT_NAME ${LIBRARY_BASENAME}_threaded
    )

  add_library( spammpack-threaded-static STATIC ${spammpack-sources} )
  add_dependencies( spammpack-threaded-static spammpack-threaded-shared )
  set_target_properties( spammpack-threaded-static
    PROPERTIES
    COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
    OUTPUT_NAME ${LIBRARY_BASENAME}_threaded
    )

  install(
    TARGETS spammpack-threaded-static spammpack-threaded-shared
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

  file( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include-threaded" )
  add_custom_command( TARGET spammpack-threaded-static
    POST_BUILD
    COMMAND cp ${spammpack-modules} "include-threaded"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

  unset( threaded-modfiles )
  foreach( mod ${spammpack-modules} )
    list( APPEND threaded-modfiles
      ${CMAKE_CURRENT_BINARY_DIR}/include-threaded/${mod}
      )
  endforeach()

  install(
    FILES ${threaded-modfiles}
    DESTINATION include-threaded
    )

endif()
