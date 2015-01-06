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

function( build_library use_OpenMP COMPILE_FLAGS )

  if( use_OpenMP )
    set( extension "threaded")
  else()
    set( extension "serial" )
  endif()

  add_library( spammpack-${extension}-shared SHARED ${spammpack-sources} )
  set_target_properties( spammpack-${extension}-shared
    PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    COMPILE_FLAGS "${COMPILE_FLAGS}"
    SOVERSION 0
    VERSION 0.0.0
    OUTPUT_NAME ${LIBRARY_BASENAME}_${extension}
    )

  add_library( spammpack-${extension}-static STATIC ${spammpack-sources} )
  # Without these additional (chaining) dependencies, a parallel build
  # gets all confused because of the way CMake builds Fortran sources:
  # The compilation is done in the src dir, but uses -o to write the
  # .o file into a subdirectory named after the target. CMake does not
  # use -J or something similar to place the generated .mod file in
  # that directory as well. Instead it copies the .mod file from the
  # src directory into the target directory.  Parallel builds are
  # still possible within a library, but not across unfortunately.
  add_dependencies( spammpack-${extension}-static spammpack-${extension}-shared )
  set_target_properties( spammpack-${extension}-static
    PROPERTIES
    COMPILE_FLAGS "${COMPILE_FLAGS}"
    OUTPUT_NAME ${LIBRARY_BASENAME}_${extension}
    )

  install(
    TARGETS spammpack-${extension}-static spammpack-${extension}-shared
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

  file( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include-${extension}" )
  add_custom_command( TARGET spammpack-${extension}-static
    POST_BUILD
    COMMAND cp ${spammpack-modules} "include-${extension}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

  unset( ${extension}-modfiles )
  foreach( mod ${spammpack-modules} )
    list( APPEND ${extension}-modfiles
      ${CMAKE_CURRENT_BINARY_DIR}/include-${extension}/${mod}
      )
  endforeach()

  install(
    FILES ${${extension}-modfiles}
    DESTINATION include-${extension}
    )

endfunction( build_library )

build_library( FALSE "" )

if( OPENMP_FOUND )
  build_library( TRUE "${OpenMP_Fortran_FLAGS}")
  add_dependencies( spammpack-threaded-shared spammpack-serial-static )
endif()
