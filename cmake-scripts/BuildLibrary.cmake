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
