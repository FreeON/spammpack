add_library( spammtests
  test_utilities.F90
  spammtests.F90
  )
target_include_directories( spammtests PRIVATE ${CMAKE_BINARY_DIR}/src )

add_executable( spamm_multiply
  spamm_multiply.F90
  )
target_link_libraries( spamm_multiply spammpack_serial_static spammtests )
target_include_directories( spamm_multiply PRIVATE ${CMAKE_BINARY_DIR}/src )
