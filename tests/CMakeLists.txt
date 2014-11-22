add_library( spammtests
  test_utilities.F90
  spammtests.F90
  )

add_executable( spamm_multiply
  spamm_multiply.F90
  )
target_link_libraries( spamm_multiply spammpack_serial_static spammtests )
include_directories( ${CMAKE_CURRENT_BINARY_DIR}/../src )
